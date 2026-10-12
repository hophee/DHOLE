#!/usr/bin/env Rscript
source("oligo_designer.R")

check <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
expect_restriction <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, restriction_site_error = identity)
  check(inherits(error, "restriction_site_error") && grepl(pattern, conditionMessage(error)),
        paste("Expected restriction rejection:", pattern))
}

sites <- c("ACTAGT", "CTGCAG")
bridge <- DEFAULT_BRIDGE_SEQUENCE
plasmid <- paste0(sites[[1]], strrep("C", 40), sites[[2]], strrep("G", 80))
sgrna <- paste0("ACG", sites[[1]], strrep("A", 20), SGRNA_SCAFFOLD, SGRNA_PRODUCT_OVERLAP)
left <- paste0(SGRNA_PRODUCT_OVERLAP, "AACCGG", bridge)
right <- paste0(bridge, "TTGGCC", sites[[2]], "CGT")
assemble <- function(sg = sgrna, lha = left, rha = right, b = bridge) {
  model_edited_ptargets(plasmid, sg, lha, b, rha, cassette_arc = "forward")
}

# Intended sites must remain, one of each, in a clean assembled plasmid.
clean <- assemble()
check(length(clean$sequences) == 1L &&
        length(restriction_site_positions(clean$sequences[[1]], sites[[1]], TRUE)) == 1L &&
        length(restriction_site_positions(clean$sequences[[1]], sites[[2]], TRUE)) == 1L,
      "Intended terminal restriction sites were rejected or lost")

# 6s: the internal PstI site is already in the left arm, not only in its primers.
expect_restriction(assemble(lha = paste0(SGRNA_PRODUCT_OVERLAP,
  "CAGCGTATTCGGATGCTGCAGCAGACCATAGAACAA", bridge)), "left_homology_arm:CTGCAG")
expect_restriction(assemble(rha = paste0(bridge, "GGACTAGTCC", sites[[2]], "CGT")),
                   "right_homology_arm:ACTAGT")
expect_restriction(assemble(sg = sub(strrep("A", 20), "AAAACTGCAGAAAAAAAAAA", sgrna)),
                   "sgRNA_N20_1:CTGCAG")
# Only the second guide is invalid; every reported variant must be checked.
expect_restriction(assemble(sg = c(sgrna, sub(strrep("A", 20), "AAAACTGCAGAAAAAAAAAA", sgrna))),
                   "sgRNA_N20_2:CTGCAG")
# A new site can straddle N20/scaffold while neither component contains it.
expect_restriction(assemble(sg = paste0("ACGACTAGT", strrep("A", 15), "CTGCA",
                                      SGRNA_SCAFFOLD, SGRNA_PRODUCT_OVERLAP)),
                   "sgRNA_N20_1:CTGCAG")
# A short custom overlap can hide a site until sgRNA and the left arm are joined.
expect_restriction(model_edited_ptargets(plasmid, "ACGACTAGTAAACTGC", paste0("GCAG", bridge),
  bridge, right, cassette_arc = "forward", sgrna_product_overlap = "GC"), "cassette_N20_1:CTGCAG")
# A short bridge can also hide a site until both arms are assembled.
expect_restriction(assemble(lha = paste0(SGRNA_PRODUCT_OVERLAP, "AAACTGC"),
  rha = paste0("GCAGAA", sites[[2]], "CGT"), b = "GC"), "cassette_N20_1:CTGCAG")

# Custom motifs: defaults are ordinary sequence; reverse complements count too.
custom_model <- model_edited_ptargets(
  paste0("GGATCC", strrep("C", 40), "AAGCTT", strrep("G", 80)),
  sub("ACTAGT", "GGATCC", sgrna),
  paste0(SGRNA_PRODUCT_OVERLAP, "ACTAGTCTGCAG", bridge), bridge,
  paste0(bridge, "TTGGCC", "AAGCTT", "CGT"),
  "GGATCC", "AAGCTT", cassette_arc = "forward")
check(length(custom_model$sequences) == 1L,
      "Defaults were incorrectly forbidden when custom sites were selected")
check_restriction_sites("ACTAGTCTGCAG", c("GGATCC", "AAGCTT"), "custom")
expect_restriction(check_restriction_sites("AAGAGACCAA", "GGTCTC", "custom"),
                   "custom:GGTCTC:positions=3")
check(identical(restriction_site_positions("AAAA", "AAA"), 1:2),
      "Overlapping linear sites were missed")
check(identical(restriction_site_positions("AGTGGACT", "ACTAGT", TRUE), 6L),
      "Site crossing the circular FASTA origin was missed")
expect_restriction(check_restriction_sites("AGTGGACT", sites, "pTarget",
  circular = TRUE), "pTarget:ACTAGT:positions=6")

# Existing genomic sites survive an unrelated edit, including coordinate shifts.
check_genome_edit_restriction_sites(DNAString("ACTAGTAAAAAACCCCCCTGCAG"),
                                    8L, 13L, "GG", sites)
# Replacing sequence with identical bases must not call an existing site new.
check_genome_edit_restriction_sites(DNAString("AACTGCAGAA"), 3L, 7L, "TGC", sites)
# Sites formed at either bridge edge and at a direct deletion junction are new.
expect_restriction(check_genome_edit_restriction_sites(DNAString("AAACTAAACCCC"),
  5L, 9L, "AGT", sites), "edited_genome:ACTAGT:positions=3")
expect_restriction(check_genome_edit_restriction_sites(DNAString("CCCCAAAAGTAA"),
  4L, 8L, "ACT", sites), "edited_genome:ACTAGT:positions=5")
expect_restriction(check_genome_edit_restriction_sites(DNAString("AAACTAAAGTAA"),
  5L, 8L, "", sites), "edited_genome:ACTAGT:positions=3")
# A lost site elsewhere must not cancel out a newly created site in a count comparison.
expect_restriction(check_genome_edit_restriction_sites(DNAString("AAAACTGCAGAAAA"),
  3L, 12L, "GGCTGCAG", sites), "edited_genome:CTGCAG:positions=6")

# Candidate search must skip a contaminated arm at every filtering level.
local({
  directory <- tempfile(); dir.create(directory)
  old_primer3 <- callPrimer3
  old_qc <- evaluate_candidate_reaction
  on.exit({
    assign("callPrimer3", old_primer3, .GlobalEnv)
    assign("evaluate_candidate_reaction", old_qc, .GlobalEnv)
    unlink(directory, recursive = TRUE)
  })
  assign("callPrimer3", function(seq, size_range, Tm, Tm_diff, name, ...) {
    starts <- if (name == "left_arm") c(1L, 8L) else 1L
    data.frame(PRIMER_LEFT_pos = starts, PRIMER_RIGHT_pos = nchar(seq),
      PRIMER_LEFT_SEQUENCE = substring(seq, starts, starts + 19L),
      PRIMER_RIGHT_SEQUENCE = reverse_complement_string(substr(seq, nchar(seq) - 19L, nchar(seq))))
  }, .GlobalEnv)
  assign("evaluate_candidate_reaction", function(input, ...) list(
    passed = TRUE, rejection_reason = "", risk_warnings = "",
    specificity = list(passed = TRUE, n_expected_products = 1,
      n_high_risk_offtarget_products = 0, n_all_offtarget_products = 0,
      n_perfect_3p_offtarget_sites = 0),
    openprimer = list(passed = TRUE, failed_soft_constraints = 0L,
      penalty = 0, max_dimer_risk = 0, abs_tm_diff = 0),
    policy = list(strict_passed = input$parameters$filtering_level != 1L,
                  blocking_passed = TRUE)
  ), .GlobalEnv)
  genome <- strrep("ACGT", 500)
  substr(genome, 974L, 979L) <- "CTGCAG"
  input <- list(genome = DNAString(genome), target_plasmid_sequence = DNAString(plasmid),
    tools = list(primer3 = file.path(.dhole_project_dir, "oligo_designer.R"), primer3_config = "unused"),
    parameters = list(left_arm = c(min = 20L, opt = 30L, max = 30L),
      right_arm = c(min = 40L, opt = 40L, max = 40L), n20_arm_min_distance = 40L,
      cds_fs = FALSE, bridge_sequence = bridge, sgrna_product_overlap = SGRNA_PRODUCT_OVERLAP,
      sgrna_scaffold = SGRNA_SCAFFOLD, ptarget_cassette_arc = "forward",
      threads = 1L, primer3_generation = primer3_generation_defaults(),
      primer3_buffer = primer3_buffer_parameters(), site1 = sites[[1]], site2 = sites[[2]]))
  feature <- list(start = 1001L, end = 1300L, length = 300L, strand = "+")
  selected <- list(n20_range = c(1130L, 1149L),
                   table = data.frame(target_sequence = "ACGTACGTACGTACGTACGTTGG"))
  for (level in 1:3) {
    input$parameters$filtering_level <- level
    result <- design_homology_arms(input, feature, selected, "cds", directory)
    check(!is.null(result) && result$pair$genome_start[[1]] == 979L,
          paste("Failed to select the clean arm at filtering level", level,
                "actual start:", result$pair$genome_start[[1]]))
    ranking <- bind_rows(result$primer_qc_trace$ranking)
    rejected <- grepl("left_homology_arm:CTGCAG", ranking$rejection_reason)
    check(any(rejected) && !any(ranking$selected[rejected]) &&
            !any(ranking$fallback_selected[rejected]),
          paste("Restriction rejection was bypassed at filtering level", level))
  }
  # Reverse-strand arms must be checked in their assembled orientation.
  bad_pair <- data.frame(genome_start = c(972L, 1300L), genome_end = c(1001L, 1339L))
  feature$strand <- "-"
  expect_restriction(check_candidate_restriction_sites(input, feature, selected, bad_pair, bridge),
                     "left_homology_arm:CTGCAG")
})

# Two guide sets with the same genomic span can differ in restriction safety.
local({
  directory <- tempfile(); dir.create(directory)
  old_arms <- design_homology_arms
  old_outputs <- write_design_outputs
  on.exit({
    assign("design_homology_arms", old_arms, .GlobalEnv)
    assign("write_design_outputs", old_outputs, .GlobalEnv)
    unlink(directory, recursive = TRUE)
  })
  assign("design_homology_arms", function(input, feature, selected, ...) {
    if (endsWith(selected$table$target_sequence[[1]], "CTGCA")) return(NULL)
    list(selected_pair_id = "clean")
  }, .GlobalEnv)
  assign("write_design_outputs", function(...) list(), .GlobalEnv)
  guides <- data.frame(n20_start = c(100L, 100L), n20_end = c(119L, 119L),
    target_sequence = c(paste0(strrep("A", 15), "CTGCA"), strrep("A", 20)),
    genomic_location = c("chr:100", "chr:97"), strand = c("+", "-"))
  result <- design_from_grna_pool(list(parameters = list(n20_mn = 1L, n20_strands = "random")),
    list(), guides, "cds", directory, file.path(directory, "design.log"))
  check(result$attempts == 2L && result$selected$table$target_sequence[[1]] == strrep("A", 20),
        "A rejected guide set hid a clean alternative with the same genomic span")
})

cat("Restriction-site tests passed\n")
