#!/usr/bin/env Rscript
source("oligo_designer.R")

check <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
expect_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  check(inherits(error, "error") && grepl(pattern, conditionMessage(error)),
        paste("Expected error:", pattern))
}

base_args <- c("--genome", "genome.fa", "--genome-annotation", "annotation.tsv",
  "--target-plasmid", "target.fa", "--output-dir", "output", "--gene", "a,b", "c")
cli <- parse_designer_args(base_args)
check(identical(cli$gene, c("a", "b", "c")) && cli$ncrna_neighbor_buffer == 20L &&
  cli$ncrna_max_n20_attempts == 3L, "Automatic target CLI/defaults incorrect")
for (argument in c("--ncrna-neighbor-buffer", "--ncrna-max-n20-attempts")) {
  for (value in c("-1", "0.5", "Inf", "NaN")) {
    expect_error(parse_designer_args(c(base_args, argument, value)), "integer")
  }
}
expect_error(parse_designer_args(c(base_args, "--ncrna-max-n20-attempts", "0")), "positive integer")

annotation <- data.frame(seqid = "chr", type = c(" CDS ", " ncRNA ", "tRNA", "unknown", NA),
  start = c(100L, 201L, 321L, 500L, 700L), stop = c(180L, 300L, 400L, 600L, 800L),
  strand = c("-", "+", "-", "+", "+"), gene = c("a", "b", "c", "CDS_word", "e"),
  locus_tag = paste0("LOC", 1:5))
auto <- list(cds = character(), ncrna = character(), gene = annotation$gene)
resolved <- resolve_targets(auto, annotation)
check(identical(resolved$class, c("cds", rep("ncrna", 4))) &&
  grepl("empty type", resolved$strategy_reason[[5]]), "Structured types did not select strategies")
untyped <- annotation; untyped$type <- NULL
check(all(resolve_targets(auto, untyped)$class == "cds"), "Untyped annotation must default to CDS")
untyped$type <- NA_character_
check(all(resolve_targets(auto, untyped)$class == "cds"), "Entirely missing types must default to CDS")
auto$cds <- c("LOC2", "b")
resolved <- resolve_targets(auto, annotation)
check(nrow(resolved) == 5L && resolved$class[resolved$feature_index == 2L] == "cds" &&
  resolved$strategy_reason[resolved$feature_index == 2L] == "explicit --cds",
  "Explicit override/alias deduplication failed")
auto$ncrna <- "b"
expect_error(resolve_targets(auto, annotation), "Conflicting --cds and --ncrna")
overridden <- resolve_targets(list(cds = character(), ncrna = "LOC1", gene = "a"), annotation)
check(nrow(overridden) == 1L && overridden$class == "ncrna", "Explicit ncRNA override failed")
missing <- resolve_targets(list(cds = c("missing", "a"), ncrna = character(), gene = character()), annotation)
check(nrow(missing) == 2L && is.na(missing$feature_index[[1]]) && missing$class[[2]] == "cds",
  "Missing CDS target prevented later targets from running")

# Read an actual TSV without a type column, not just an in-memory annotation.
local({
  path <- tempfile(fileext = ".tsv"); on.exit(unlink(path))
  write_tsv(annotation[, setdiff(names(annotation), "type")], path)
  check(all(is.na(read_genome_annotation(path, "bakta")$type)), "Untyped TSV rejected")
})

input <- list(annotation = annotation, genome = DNAString(strrep("ACGT", 250)),
  genome_contig = "chr", parameters = list(ncrna_neighbor_buffer = 20L))
feature <- feature_record(input, "b")
bounds <- ncrna_deletion_bounds(input, feature)
check(bounds$start == 201L && bounds$end == 300L && bounds$reason == "",
  "Inclusive 20-nt buffers or opposite-strand neighbors incorrect")
input$annotation$stop[[1]] <- 181L
check(grepl("cannot preserve neighbor buffer", ncrna_deletion_bounds(input, feature)$reason),
  "A 19-nt gap silently reduced the buffer")
input$annotation$start[[3]] <- 320L
check(grepl("allowed 202-299", ncrna_deletion_bounds(input, feature)$reason),
  "Tight neighbors on both sides were not reported")
input$annotation$stop[[1]] <- 202L
check(grepl("overlaps.*LOC1", ncrna_deletion_bounds(input, feature)$reason),
  "Overlapping annotation was ignored")
input$annotation <- rbind(annotation, transform(annotation[2L, ], locus_tag = "other", gene = "other"))
check(grepl("overlaps.*other", ncrna_deletion_bounds(input, feature)$reason),
  "Fully coincident different feature was ignored")
input$annotation <- annotation[2L, ]; input$annotation$seqid <- "chr"
feature <- feature_record(input, "b")
check(identical(unname(unlist(ncrna_deletion_bounds(input, feature)[c("start", "end")])),
  c(1L, 1000L)), "Missing neighbors must use contig edges without buffer")
input$annotation <- rbind(input$annotation, transform(annotation[1L, ], seqid = "other_contig"))
check(ncrna_deletion_bounds(input, feature)$start == 1L, "Other contig constrained deletion")
check(identical(cut_interval(list(start = 100L, end = 2000L, length = 1901L), "ncrna", 3000L),
  c(100L, 2000L)), "Long ncRNA was reduced to its middle")

# Forced Primer3 endpoints must not admit a primer flagged as violating constraints.
local({
  wrapper <- callPrimer3
  environment(wrapper) <- new.env(parent = environment(callPrimer3))
  input_file <- tempfile(); output_file <- tempfile()
  captured <- character()
  assign("system2", function(command, args, stdout, stderr) {
    if (any(grepl("--about", args, fixed = TRUE))) return("primer3 release 2.6.1")
    captured <<- readLines(input_file)
    writeLines(c("PRIMER_PAIR_NUM_RETURNED=1", "PRIMER_LEFT_0_PROBLEMS=GC content failed", "="), output_file)
    character()
  }, environment(wrapper))
  rejected <- wrapper(strrep("ACGT", 100), name = "forced", force_left_start = 1L,
    force_right_start = 400L, p3_input_file = input_file, p3_output_file = output_file,
    primer3 = "synthetic-tool", settings = "settings", thermo.param = "config")
  check(nrow(rejected) == 0L && "SEQUENCE_FORCE_LEFT_START=1" %in% captured &&
    "SEQUENCE_FORCE_RIGHT_START=400" %in% captured,
    "Forced endpoint tags or rejection of Primer3 constraint violations failed")
  check(!file.exists(input_file) && !file.exists(output_file), "Primer3 temporary files leaked")
})

local({
  directory <- tempfile(); dir.create(directory); on.exit(unlink(directory, recursive = TRUE))
  input <- list(annotation = annotation, genome = DNAString(strrep("ACGT", 250)),
    genome_contig = "chr", output_dir = directory, parameters = list(ncrna_neighbor_buffer = 20L))
  input$annotation$stop[[1]] <- 210L
  expect_error(design_target(input, "genome", "b", "ncrna", "annotation: ncrna"), "overlaps")
  report <- readLines(file.path(directory, "TechReport", "b_results", "report.tsv"))
  check(all(c("strategy\tncrna", "strategy_reason\tannotation: ncrna", "target_start\t201",
    "target_end\t300", "ncrna_neighbor_buffer\t20", "n20_sets_checked\t0",
    "failure_stage\tdeletion_bounds") %in% report) && any(grepl("failure_reason.*LOC1", report)),
    "Neighbor-conflict failure lost strategy, geometry or diagnostic metadata")
})

make_input <- function() list(
  genome = DNAString(strrep("ACGT", 500)),
  target_plasmid_sequence = DNAString(paste0("ACTAGT", strrep("C", 40), "CTGCAG", strrep("G", 80))),
  primer_qc_cache = new.env(parent = emptyenv()),
  tools = list(primer3 = file.path(.dhole_project_dir, "oligo_designer.R"), primer3_config = "unused"),
  parameters = list(left_arm = c(min = 20L, opt = 30L, max = 30L),
    right_arm = c(min = 20L, opt = 30L, max = 30L), n20_arm_min_distance = 0L,
    cds_fs = FALSE, ncrna_fs = FALSE, bridge_sequence = DEFAULT_BRIDGE_SEQUENCE,
    sgrna_product_overlap = SGRNA_PRODUCT_OVERLAP, sgrna_scaffold = SGRNA_SCAFFOLD,
    ptarget_cassette_arc = "forward", threads = 1L, filtering_level = 3L,
    primer3_generation = primer3_generation_defaults(), primer3_buffer = primer3_buffer_parameters(),
    primer_qc = primer_qc_defaults(), site1 = "ACTAGT", site2 = "CTGCAG",
    n20_mn = 1L, n20_strands = "random", ncrna_max_n20_attempts = 3L))
make_feature <- function(strand = "+", lower = 1001L, upper = 1103L) list(
  index = 1L, start = 1001L, end = 1100L, length = 100L, strand = strand,
  deletion_bounds = list(start = lower, end = upper, reason = ""))
guides <- data.frame(n20_start = 1040:1043, n20_end = 1059:1062,
  target_sequence = paste0(c("A", "C", "G", "T"), strrep("ACG", 6), "A", "TGG"),
  genomic_location = paste0("chr:", 1040:1043), strand = "+")
selected <- list(n20_range = c(1040L, 1059L), table = guides[1L, ])

local({
  directory <- tempfile(); dir.create(directory)
  old_tool <- run_tool
  on.exit({ assign("run_tool", old_tool, .GlobalEnv); unlink(directory, recursive = TRUE) })
  requested <- character()
  assign("run_tool", function(command, args, stdout, ...) {
    requested <<- c(requested, args[[which(args == "-Target") + 1L]])
    writeLines("synthetic CHOPCHOP result", stdout)
  }, .GlobalEnv)
  feature <- make_feature(lower = 980L, upper = 1120L); feature$contig <- "chr"
  result <- run_chopchop(make_input(), "genome", feature, "ncrna", directory)
  check(identical(requested, "chr:980-1120") && identical(result, c(980L, 1120L)),
    "CHOPCHOP did not use the maximum permitted ncRNA region")
})

local({
  directory <- tempfile(); dir.create(directory)
  old_primer3 <- callPrimer3; old_qc <- evaluate_candidate_reaction
  old_outputs <- write_design_outputs
  on.exit({
    assign("callPrimer3", old_primer3, .GlobalEnv)
    assign("evaluate_candidate_reaction", old_qc, .GlobalEnv)
    assign("write_design_outputs", old_outputs, .GlobalEnv)
    unlink(directory, recursive = TRUE)
  })
  generation_calls <- 0L
  assign("callPrimer3", function(seq, size_range, Tm, Tm_diff, name, ...) {
    generation_calls <<- generation_calls + 1L
    data.frame(PRIMER_LEFT_pos = 1L, PRIMER_RIGHT_pos = nchar(seq),
      PRIMER_LEFT_SEQUENCE = substr(seq, 1L, 20L),
      PRIMER_RIGHT_SEQUENCE = reverse_complement_string(substr(seq, nchar(seq) - 19L, nchar(seq))))
  }, .GlobalEnv)
  defect <- "none"
  assign("evaluate_candidate_reaction", function(input, ...) {
    specificity <- list(passed = !defect %in% c("offtarget", "missing"), n_expected_products = if (defect == "missing") 0L else 1L,
      n_high_risk_offtarget_products = as.integer(defect == "offtarget"),
      n_all_offtarget_products = as.integer(defect == "offtarget"), n_perfect_3p_offtarget_sites = 0L)
    openprimer <- list(passed = defect == "none", failed_soft_constraints = 0L,
      penalty = as.integer(defect != "none"), max_dimer_risk = 0, abs_tm_diff = if (defect == "tm") 10 else 0,
      rejection_reason = if (defect == "none") "" else defect,
      tm_diff_evaluated = TRUE)
    policy <- evaluate_filtering_policy(specificity, openprimer, input$parameters$filtering_level, input$parameters$primer_qc)
    list(passed = policy$blocking_passed, rejection_reason = policy$blocking_reasons,
      risk_warnings = policy$warnings, specificity = specificity, openprimer = openprimer, policy = policy)
  }, .GlobalEnv)

  # Initial full deletion and asymmetric expansion work in both orientations.
  for (strand in c("+", "-")) {
    feature <- make_feature(strand)
    input <- make_input()
    exact <- design_homology_arms(input, feature, selected, "ncrna", directory)
    check(identical(as.integer(exact$ticks[2:3]), c(1000L, 1101L)), "Initial ncRNA deletion was not exact")
    near_right <- selected; near_right$n20_range <- c(1085L, 1103L)
    expanded <- design_homology_arms(input, feature, near_right, "ncrna", directory)
    check(identical(as.integer(expanded$ticks[2:3]), c(1000L, 1104L)),
      "Right expansion stopped at the left limit or discarded the original left boundary")
    before <- generation_calls
    design_homology_arms(input, feature, near_right, "ncrna", directory)
    check(generation_calls == before, "Primer3 candidates were not reused across N20 attempts")
    near_left <- selected; near_left$n20_range <- c(998L, 1017L)
    mirrored <- design_homology_arms(input, make_feature(strand, 998L, 1100L), near_left, "ncrna", directory)
    check(identical(as.integer(mirrored$ticks[2:3]), c(997L, 1101L)), "Independent left expansion failed")
  }

  # Output/screening failures must try other geometries of the same N20 set.
  output_calls <- 0L
  assign("write_design_outputs", function(input, feature, selected, arms, ...) {
    output_calls <<- output_calls + 1L
    if (arms$ticks[[3]] < 1103L) stop("screening primer QC defect")
    list(accepted = TRUE)
  }, .GlobalEnv)
  result <- design_from_grna_pool(make_input(), make_feature(), guides, "ncrna", directory,
    file.path(directory, "design.log"))
  check(result$attempts == 1L && result$arms$ticks[[3]] == 1103L && output_calls == 3L,
    "Output failure skipped an available minimal expansion of the same guide")

  first_set_boundaries <- integer()
  generation_calls <- 0L
  assign("write_design_outputs", function(input, feature, selected, arms, ...) {
    if (selected$table$n20_start[[1]] == 1040L) {
      first_set_boundaries <<- c(first_set_boundaries, arms$ticks[[3]])
      stop("screening primer QC defect for first N20")
    }
    list(accepted = TRUE)
  }, .GlobalEnv)
  retried <- design_from_grna_pool(make_input(), make_feature(), guides, "ncrna", directory,
    file.path(directory, "design.log"))
  check(retried$attempts == 2L && retried$selected$table$n20_start[[1]] == 1041L &&
    identical(as.integer(first_set_boundaries), 1101:1104) && generation_calls == 5L,
    "Next N20 was not tried after exhaustive expansion, or regenerated cached arms")

  minimal <- design_homology_arms(make_input(), make_feature(lower = 998L), selected,
    "ncrna", directory, validate_design = function(arms) {
      left_extension <- 1000L - arms$ticks[[2]]
      right_extension <- arms$ticks[[3]] - 1101L
      if (right_extension >= 3L || (left_extension >= 2L && right_extension >= 2L)) list(ok = TRUE) else NULL
    })
  check(identical(as.integer(minimal$ticks[2:3]), c(1000L, 1104L)),
    "Search preferred a symmetric 4-nt expansion over an asymmetric 3-nt expansion")
  distances <- selected; distances$n20_range <- c(1041L, 1061L)
  input <- make_input(); input$parameters$n20_arm_min_distance <- 40L
  offset <- design_homology_arms(input, make_feature(upper = 1101L), distances, "ncrna", directory)
  check(identical(as.integer(offset$ticks[2:3]), c(1000L, 1102L)),
    "Actual N20-to-arm offset was not enforced")

  # Same QC defect, same filtering-level outcome in both geometry strategies.
  for (level in 1:3) {
    for (problem in c("none", "missing", "offtarget", "tm")) {
      defect <- problem
      input <- make_input(); input$parameters$filtering_level <- level
      cds <- design_homology_arms(input, make_feature(), selected, "cds", directory)
      nc <- design_homology_arms(input, make_feature(upper = 1100L), selected, "ncrna", directory)
      check(is.null(cds) == is.null(nc), paste("QC parity failure", level, problem))
    }
  }
  defect <- "none"
  for (level in 1:3) {
    input <- make_input(); input$genome <- DNAString(strrep("CTGCAG", 334))
    input$parameters$filtering_level <- level
    for (strategy in c("cds", "ncrna")) {
      trace <- new_primer_qc_trace()
      rejected <- design_homology_arms(input, make_feature(upper = 1100L), selected, strategy, directory, trace)
      check(is.null(rejected) && any(grepl("CTGCAG", bind_rows(trace$ranking)$rejection_reason)),
        paste("Restriction site bypassed", strategy, level))
    }
  }

  # Count total sets, stop at the cap, and retain diagnostic reports on failure.
  assign("callPrimer3", function(...) data.frame(), .GlobalEnv)
  input <- make_input()
  expect_error(design_from_grna_pool(input, make_feature(upper = 1100L), guides,
    "ncrna", directory, file.path(directory, "design.log")), "checked: 3; attempt_limit_reached")
  report <- readLines(file.path(directory, "report.tsv"))
  check(any(report == "n20_sets_checked\t3") && any(report == "n20_search_termination\tattempt_limit_reached"),
    "Attempt limit missing from failure report")
  input <- make_input()
  expect_error(design_from_grna_pool(input, make_feature(upper = 1100L), guides[1:2, ],
    "ncrna", directory, file.path(directory, "design.log")), "checked: 2; candidates_exhausted")
  input <- make_input(); input$parameters$ncrna_max_n20_attempts <- 1L
  expect_error(design_from_grna_pool(input, make_feature(upper = 1100L), guides,
    "ncrna", directory, file.path(directory, "design.log")), "checked: 1; attempt_limit_reached")
  edge <- make_feature(); edge$start <- 1L
  expect_error(design_homology_arms(make_input(), edge, selected, "ncrna", directory), "Genome boundary")
  edge <- make_feature(); edge$end <- 2000L
  expect_error(design_homology_arms(make_input(), edge, selected, "ncrna", directory), "Genome boundary")
})

cat("ncRNA strategy and automatic target tests passed\n")
