#!/usr/bin/env Rscript
source("oligo_designer.R")

check <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
expect_error <- function(expression, pattern, class = "error") {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  check(inherits(error, class) && grepl(pattern, conditionMessage(error)),
        paste("Expected", class, pattern))
}

# 3, 12: the same circular molecule must retain the same backbone after RC/rotation.
cassette <- paste0("GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT", "CCCC",
                   reverse_complement_string("AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC"))
plasmid <- paste0("ACTAGT", cassette, "CTGCAG", "AGTCCG")
for (sequence in c(plasmid, reverse_complement_string(plasmid),
                   paste0(substr(plasmid, 20L, nchar(plasmid)), substr(plasmid, 1L, 19L)))) {
  pair <- find_oriented_restriction_pair(sequence, "ACTAGT", "CTGCAG")
  check(identical(pair$backbone, "AGTCCG"), "pTarget retained the wrong arc")
}
expect_error(find_oriented_restriction_pair("ACTAGTCCCCCTGCAGGGG", "ACTAGT", "CTGCAG"),
             "sgRNA PCR")
expect_error(find_oriented_restriction_pair(paste0("ACTAGTGGGCTGCAG", cassette),
                                           "ACTAGT", "CTGCAG"), "не лежит между")
check(identical(circular_match_positions("AAAA", "AAA"), 1:4),
      "Overlapping/circular sites were missed")

# 4, 13: separately specified genomic N20 and PAM on both strands.
local({
  path <- tempfile(fileext = ".tsv")
  on.exit(unlink(path))
  target <- "ACGTACGTACGTACGTACGTTGG"
  sequence <- paste0(strrep("A", 49), target, strrep("A", 27),
                     reverse_complement_string(target), strrep("A", 77))
  table <- data.frame(genomic_location = c("chr:50", "chr:100"),
                      target_sequence = target, strand = c("+", "-"),
                      self_complementarity = 0L, mm0 = 0L)
  write_tsv(table, path)
  feature <- list(start = 1L, end = 200L, length = 200L)
  pool <- filter_grnas(path, feature, "ncrna", 0L, DNAString(sequence))
  plus <- pool[pool$strand == "+", ]
  minus <- pool[pool$strand == "-", ]
  check(identical(as.integer(plus[1, c("n20_start", "n20_end", "pam_start", "pam_end")]),
                  c(50L, 69L, 70L, 72L)), "Incorrect plus N20/PAM")
  check(identical(as.integer(minus[1, c("n20_start", "n20_end", "pam_start", "pam_end")]),
                  c(103L, 122L, 100L, 102L)), "Incorrect minus N20/PAM")
  check(all(prepare_grna_pool(pool, 1L, "minus")$strand == "-") &&
          all(prepare_grna_pool(pool, 1L, "plus")$strand == "+") &&
          nrow(prepare_grna_pool(pool, 1L, "both")) == 2L,
        "Single-N20 strand policy was ignored")
  expect_error(filter_grnas(path, feature, "ncrna", 0L, DNAString(strrep("A", 200))),
               "не совпадает")
})

# 8, 22: missing products remain a QC rejection; oversized products keep coordinates.
local({
  path <- tempfile(fileext = ".fasta")
  on.exit(unlink(path))
  forward <- "ACGTTGCAACGTTCGATCGA"
  reverse <- "TGCACCGATGTTACGTCAGT"
  sequence <- paste0(forward, strrep("A", 100), reverse_complement_string(reverse))
  writeXStringSet(DNAStringSet(c(chr = sequence)), path)
  refs <- make_specificity_references(path, path)
  config <- primer_qc_defaults()
  config$max_mismatches <- 0L
  config$max_product_size <- 80L
  expected <- list(reference_id = "genome::chr", start = 1L, end = 140L, size = 140L)
  cache <- new.env(parent = emptyenv())
  result <- evaluate_pair_specificity(forward, reverse, refs, expected, config, cache)
  check(!result$passed && nrow(result$amplicons) > 0L &&
          all(result$amplicons$invalid_size) && all(is.na(result$amplicons$sequence)),
        "Oversized amplicons lost their trace or materialized sequences")
  cached <- ls(cache)
  again <- evaluate_pair_specificity(forward, reverse, refs, expected, config, cache)
  check(identical(result, again) && identical(ls(cache), cached) && length(cached) == 2L,
        "Binding-site cache changed specificity results")
  missing <- evaluate_pair_specificity(strrep("G", 20), strrep("C", 20), refs,
                                      expected, config)
  trace <- new_primer_qc_trace()
  append_primer_qc_trace(trace, "missing", "pair1", specificity = missing)
  check(!missing$passed && grepl("expected_product_not_found", missing$rejection_reason) &&
          nrow(trace$binding_sites[[1]]) == 0L && nrow(trace$amplicons[[1]]) == 0L,
        "Empty specificity trace is not a normal rejection")
})

# 9: opt == max must still visit an interior deletion for a short CDS.
local({
  directory <- tempfile(); dir.create(directory)
  old_primer3 <- callPrimer3
  old_qc <- evaluate_candidate_reaction
  on.exit({
    assign("callPrimer3", old_primer3, .GlobalEnv)
    assign("evaluate_candidate_reaction", old_qc, .GlobalEnv)
    unlink(directory, recursive = TRUE)
  })
  calls <- 0L
  assign("callPrimer3", function(seq, ...) {
    calls <<- calls + 1L
    data.frame(PRIMER_LEFT_pos = 1L, PRIMER_RIGHT_pos = nchar(seq),
               PRIMER_LEFT_SEQUENCE = substr(seq, 1L, 20L),
               PRIMER_RIGHT_SEQUENCE = reverse_complement_string(substr(seq, nchar(seq)-19L, nchar(seq))))
  }, .GlobalEnv)
  assign("evaluate_candidate_reaction", function(...) stop("reached_valid_geometry"), .GlobalEnv)
  input <- list(genome = DNAString(strrep("ACGT", 500)),
                tools = list(primer3 = file.path(.dhole_project_dir, "oligo_designer.R"),
                             primer3_config = "unused"),
                parameters = list(left_arm = c(min = 30L, opt = 30L, max = 30L),
                                  right_arm = c(min = 40L, opt = 40L, max = 40L),
                                  n20_arm_min_distance = 40L, cds_fs = FALSE,
                                  primer3_buffer = primer3_buffer_parameters(), site2 = "CTGCAG"))
  expect_error(design_homology_arms(input, list(start = 1001L, end = 1300L,
               length = 300L, strand = "+"), list(n20_range = c(1130L, 1149L)),
               "cds", directory), "reached_valid_geometry")
  check(calls == 4L, "Short CDS did not reach the second geometry at fixed arm lengths")
  assign("callPrimer3", function(...) data.frame(), .GlobalEnv)
  check(is.null(design_homology_arms(input, list(start = 1001L, end = 1300L,
             length = 300L, strand = "+"), list(n20_range = c(1130L, 1149L)),
             "cds", directory)), "Zero Primer3 candidates broke the homology search")
})

# 11: an expected product plus a shorter product is not unique.
local({
  forward <- "ACGTTGCAACGTTCGATCGA"
  reverse <- "TGCACCGATGTTACGTCAGT"
  template <- paste0(forward, strrep("A", 40), reverse_complement_string(reverse),
                      strrep("A", 40), reverse_complement_string(reverse))
  expect_error(simulate_full_primer_pcr("multiple", "test", "synthetic", template,
               "F", "R", forward, reverse, forward, reverse, 45), "однозначный PCR")
})

# 16: indexes follow FASTA content, completeness, and successful completion.
local({
  directory <- tempfile(); dir.create(directory)
  fasta <- file.path(directory, "genome.fasta")
  writeLines(c(">chr", "ACGT"), fasta)
  old_run <- run_tool
  calls <- 0L
  fail_build <- FALSE
  assign("run_tool", function(command, args, ...) {
    calls <<- calls + 1L
    if (command == "faToTwoBit") writeLines("index", args[[2]]) else {
      if (fail_build) stop("interrupted index")
      for (suffix in c(".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt")) {
        writeLines("index", paste0(args[[2]], suffix))
      }
    }
  }, .GlobalEnv)
  on.exit({ assign("run_tool", old_run, .GlobalEnv); unlink(directory, recursive = TRUE) })
  input <- list(genome_path = fasta, output_dir = directory)
  assets <- prepare_chopchop_assets(input)
  prepare_chopchop_assets(input)
  check(calls == 2L, "Matching indexes were rebuilt")
  writeLines(c(">chr", "TGCA"), fasta)
  prepare_chopchop_assets(input)
  check(calls == 4L, "Changed FASTA reused stale indexes")
  unlink(file.path(assets$directory, "genome.rev.2.ebwt"))
  fail_build <- TRUE
  expect_error(prepare_chopchop_assets(input), "interrupted index")
  check(!file.exists(file.path(assets$directory, "genome.fasta.md5")),
        "Failed indexing left a success marker")
  fail_build <- FALSE
  prepare_chopchop_assets(input)
  check(calls == 8L, "Incomplete indexing was not retried")
})

# 17: importing in another working directory changes no global cwd.
local({
  previous <- getwd()
  on.exit(setwd(previous))
  setwd(tempdir())
  imported <- new.env(parent = .GlobalEnv)
  source(file.path(.dhole_project_dir, "oligo_designer.R"), local = imported)
  check(identical(getwd(), tempdir()) && is.function(imported$callPrimer3) &&
          identical(imported$.dhole_project_dir, .dhole_project_dir),
        "Absolute import depends on cwd or changes it")
})

# 19: distinguish no candidates, PRIMER_ERROR and a nonzero process exit.
local({
  directory <- tempfile("primer3 paths with spaces "); dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  executable <- file.path(directory, "fake primer3")
  input_file <- file.path(directory, "input file")
  output_file <- file.path(directory, "output file")
  report <- file.path(directory, "report file")
  make_tool <- function(body) {
    writeLines(c("#!/bin/sh",
      'if [ "$1" = "--about" ]; then echo "primer3 release 2.6.1"; exit 0; fi',
      'for arg do case "$arg" in --output=*) output=${arg#--output=} ;; esac; done',
      body), executable)
    Sys.chmod(executable, "0755")
  }
  run <- function() callPrimer3(strrep("ACGT", 100), name = "test", primer3 = executable,
    settings = file.path(directory, "settings file"), thermo.param = directory,
    p3_input_file = input_file, p3_output_file = output_file, report = report)
  make_tool('printf "PRIMER_PAIR_NUM_RETURNED=0\\n=\\n" > "$output"')
  check(is.data.frame(run()), "Zero Primer3 candidates are not a normal result")
  check(!file.exists(input_file) && !file.exists(output_file), "No-candidate temp files leaked")
  make_tool('printf "PRIMER_ERROR=broken setting\\n=\\n" > "$output"')
  expect_error(run(), "broken setting", "primer3_error")
  check(!file.exists(input_file) && !file.exists(output_file), "PRIMER_ERROR temp files leaked")
  make_tool('echo "failed deliberately" >&2; exit 7')
  expect_error(run(), "exit status 7.*failed deliberately", "primer3_error")
  check(!file.exists(input_file) && !file.exists(output_file) &&
          grepl("failed deliberately", paste(readLines(paste0(report, ".stderr.log")), collapse = "")),
        "Process failure lost stderr or leaked temp files")
})

# 22: preserve IUPAC meaning and the agreed optional frame policy.
check(identical(reverse_complement_string("ARYKMBDHVN"), "NBDHVKMRYT"), "Incorrect IUPAC RC")
check(identical(vapply(198:200, function(n) design_bridge("cds", n), character(1)),
                c("ATGACTGCCCGCAAG", "ATGACTGCCCGCAA", "ATGACTGCCCGCA")),
      "Existing optional-frame bridge rule changed")
message("Targeted regressions passed")
