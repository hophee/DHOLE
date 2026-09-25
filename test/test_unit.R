#!/usr/bin/env Rscript

source("oligo_designer.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message, call. = FALSE)
  }
}

assert_error <- function(expression, pattern) {
  error <- tryCatch(
    {
      force(expression)
      NULL
    },
    error = identity
  )
  assert_true(!is.null(error), paste("Expected error:", pattern))
  assert_true(
    grepl(pattern, conditionMessage(error)),
    sprintf(
      "Error '%s' does not match '%s'",
      conditionMessage(error),
      pattern
    )
  )
}

base_args <- c(
  "--genome", "genome.fasta",
  "--genome-annotation", "genome.gff",
  "--target-plasmid", "target.fasta",
  "--output-dir", "output",
  "--cds", "geneA"
)

defaults <- parse_designer_args(base_args)
assert_true(
  defaults$annotation_format == "gff",
  "annotation default is not gff"
)
assert_true(
  defaults$chopchop_script == "chopchop/chopchop.py",
  "CHOPCHOP path default is incorrect"
)
assert_true(
  defaults$chopchop_python == "chopchop-python",
  "CHOPCHOP Python default is incorrect"
)
assert_true(
  defaults$primer3 == "primer3/src/primer3_core",
  "Primer3 path default is incorrect"
)
assert_true(defaults$filtering_level == 2L, "filtering_level default is not 2")
assert_true(defaults$n20_mn == 1L, "n20_mn default is not 1")
assert_true(defaults$n20_strands == "random", "strand default is not random")
assert_true(identical(defaults$n20_offtarget, 0L), "MM default is not 0")
assert_true(defaults$site1 == "ACTAGT", "site1 default is not SpeI")
assert_true(defaults$site2 == "CTGCAG", "site2 default is not PstI")
assert_true(
  defaults$ptarget_cassette_arc == "shortest",
  "pTarget cassette defaults are incorrect"
)
assert_true(
  identical(defaults$sgrna_annealing_temp_c, 60),
  "sgRNA annealing temperature default is incorrect"
)
assert_true(!defaults$cds_fs && !defaults$ncrna_fs, "FS defaults are not FALSE")
assert_true(
  identical(unname(defaults$left_arm), c(300L, 350L, 400L)),
  "left arm defaults are incorrect"
)
assert_true(
  identical(unname(defaults$right_arm), c(400L, 450L, 500L)),
  "right arm defaults are incorrect"
)

minus_positions <- add_genome_positions(
  data.frame(PRIMER_LEFT_pos = 1L, PRIMER_RIGHT_pos = 20L),
  data.frame(PRIMER_LEFT_pos = 1L, PRIMER_RIGHT_pos = 20L),
  list(left_end = 100L, right_end = 50L),
  FALSE
)
assert_true(
  identical(
    unname(as.integer(minus_positions$left[1, c("genome_start", "genome_end")])),
    c(81L, 100L)
  ),
  "minus-strand left-arm coordinates are not one-based"
)
assert_true(
  defaults$n20_arm_min_distance == 40L,
  "N20-to-arm distance default is not 40"
)

configured <- parse_designer_args(c(
  base_args,
  "--filtering-level", "3",
  "--n20-mn", "3",
  "--n20-strands", "both",
  "--n20-offtarget", "0,2,4",
  "--site1", "ggatcc",
  "--site2", "aagctt",
  "--ptarget-cassette-arc", "forward",
  "--sgrna-annealing-temp-c", "62.5",
  "--cds-fs"
))
assert_true(configured$filtering_level == 3L, "filtering_level was not parsed")
assert_true(configured$n20_mn == 3L, "n20_mn was not parsed")
assert_true(configured$n20_strands == "both", "strand mode was not parsed")
assert_true(
  identical(configured$n20_offtarget, c(0L, 2L, 4L)),
  "MM thresholds were not parsed"
)
assert_true(configured$cds_fs, "cds_fs flag was not parsed")
assert_true(
  configured$site1 == "GGATCC" && configured$site2 == "AAGCTT",
  "Restriction sites were not parsed and normalized"
)
assert_true(
  configured$ptarget_cassette_arc == "forward",
  "Custom pTarget cassette parameters were not parsed"
)
assert_true(
  identical(configured$sgrna_annealing_temp_c, 62.5),
  "Custom sgRNA annealing temperature was not parsed"
)

duplicate_targets <- base_args
duplicate_targets[[length(duplicate_targets)]] <- "geneA,geneA"
assert_true(
  identical(parse_designer_args(duplicate_targets)$cds, "geneA"),
  "Duplicate targets were not removed"
)

assert_error(
  parse_designer_args(c(base_args, "--n20-offtarget", "0,-1")),
  "неотрицательные"
)
assert_error(
  parse_designer_args(c(base_args, "--left-arm-min", "401")),
  "min <= opt <= max"
)
assert_error(
  parse_designer_args(c(base_args, "--filtering-level", "4")),
  "1, 2 или 3"
)
assert_error(
  parse_designer_args(c(base_args, "--sgrna-annealing-temp-c", "100")),
  "между 0 и 100"
)
assert_error(
  parse_designer_args(c(base_args, "--ptarget-cassette-arc", "guess")),
  "shortest, forward|reverse"
)
assert_error(
  parse_designer_args(c(base_args, "--site1", "ACTNGT")),
  "A/C/G/T"
)
assert_error(
  parse_designer_args(c(base_args, "--site2", "ACTAGT")),
  "должны быть разными"
)

feature <- list(start = 100L, end = 299L, length = 200L)
grna_path <- tempfile(fileext = ".tsv")
write_tsv(
  data.frame(
    Rank = 1:4,
    `Target sequence` = rep("AAAAAAAAAAAAAAAAAAAATGG", 4),
    `Genomic location` = paste0("contig:", c(150, 160, 170, 180)),
    Strand = c("+", "-", "+", "-"),
    `Self-complementarity` = c(0, 0, 0, 1),
    MM0 = c(0, 0, 1, 0),
    MM1 = c(1, 3, 0, 0),
    check.names = FALSE
  ),
  grna_path
)
filtered <- filter_grnas(grna_path, feature, "ncrna", c(0L, 2L))
assert_true(nrow(filtered) == 1L, "MM/self-complementarity filtering failed")
assert_true(filtered$strand[[1]] == "+", "Unexpected N20 survived filtering")
assert_error(
  filter_grnas(grna_path, feature, "ncrna", c(0L, 1L, 2L)),
  "только 2 колонок MM"
)

pool <- data.frame(
  rank = 1:4,
  strand = c("+", "-", "-", "+"),
  genomic_location = paste0("contig:", 1:4),
  n20_start = 1:4,
  n20_end = 20:23
)
original_pool <- pool
visited <- 0L
fallback <- visit_grna_sets(pool, 2L, "both", function(candidate) {
  visited <<- visited + 1L
  if (identical(candidate$rank, c(1L, 3L))) candidate else NULL
})
assert_true(visited == 2L, "N20 fallback did not visit the next valid set")
assert_true(
  identical(fallback$rank, c(1L, 3L)),
  "N20 fallback returned the wrong set"
)
assert_true(identical(pool, original_pool), "N20 pool was mutated")

left <- data.frame(
  genome_start = 50L,
  genome_end = 100L,
  arm = "left"
)
right <- data.frame(
  genome_start = 201L,
  genome_end = 250L,
  arm = "right"
)
selected <- list(n20_range = c(141L, 160L))
pair <- choose_pair(
  list(left = left, right = right),
  TRUE,
  list(length = 500L),
  selected,
  minimum_n20_distance = 40L,
  restrict_frame_shift = FALSE
)
assert_true(!is.null(pair), "Exact 40 nt N20 distance was rejected")
pair_with_frame_limit <- choose_pair(
  list(left = left, right = right),
  TRUE,
  list(length = 500L),
  selected,
  minimum_n20_distance = 40L,
  restrict_frame_shift = TRUE
)
assert_true(
  is.null(pair_with_frame_limit),
  "Non-divisible deletion passed frame restriction"
)
minus_pair <- choose_pair(
  list(
    left = data.frame(genome_start = 201L, genome_end = 250L),
    right = data.frame(genome_start = 50L, genome_end = 100L)
  ),
  FALSE,
  list(length = 500L),
  selected,
  minimum_n20_distance = 40L,
  restrict_frame_shift = FALSE
)
assert_true(!is.null(minus_pair), "Minus-strand N20 distance was rejected")

write_reference_fasta <- function(sequence, name) {
  path <- tempfile("2pac-reference-", fileext = ".fasta")
  writeXStringSet(DNAStringSet(setNames(sequence, name)), path)
  path
}

make_test_references <- function(genome, target_plasmid, cas_plasmid = NULL) {
  make_specificity_references(
    write_reference_fasta(genome, "chromosome"),
    write_reference_fasta(target_plasmid, "pTarget"),
    if (is.null(cas_plasmid)) NULL else {
      write_reference_fasta(cas_plasmid, "pCas")
    }
  )
}

replace_base <- function(sequence, position) {
  bases <- strsplit(sequence, "", fixed = TRUE)[[1]]
  bases[[position]] <- if (bases[[position]] == "A") "C" else "A"
  paste(bases, collapse = "")
}

specificity_config <- primer_qc_defaults()
specificity_config$min_product_size <- 1L
specificity_config$max_product_size <- 500L
specificity_config$max_mismatches <- 1L
forward_probe <- "ACGTTGCAACGTTCGATCGA"
reverse_probe <- "TGCACCGATGTTACGTCAGT"
reverse_site <- as.character(reverseComplement(DNAString(reverse_probe)))
expected_sequence <- paste0(forward_probe, strrep("A", 20L), reverse_site)
references <- make_test_references(
  expected_sequence,
  strrep("G", 120L),
  strrep("C", 120L)
)
expected_product <- list(
  reference_id = "genome::chromosome",
  start = 1L,
  end = nchar(expected_sequence),
  size = nchar(expected_sequence),
  allowed_products = 1L
)
specificity <- evaluate_pair_specificity(
  forward_probe,
  reverse_probe,
  references,
  expected_product,
  specificity_config
)
assert_true(specificity$passed, "The exact expected amplicon was rejected")
assert_true(
  specificity$n_expected_products == 1L,
  "The exact expected amplicon was not identified uniquely"
)

multi_sequence <- paste0(
  forward_probe,
  strrep("A", 20L),
  reverse_site,
  strrep("C", 20L),
  reverse_site
)
multi_references <- make_test_references(
  multi_sequence,
  strrep("G", 160L),
  strrep("C", 160L)
)
multi_specificity <- evaluate_pair_specificity(
  forward_probe,
  reverse_probe,
  multi_references,
  list(
    reference_id = "genome::chromosome",
    start = 1L,
    end = nchar(expected_sequence),
    size = nchar(expected_sequence),
    allowed_products = 1L
  ),
  specificity_config
)
assert_true(
  nrow(multi_specificity$amplicons) == 2L,
  "Exhaustive pairing did not find both amplicons"
)
assert_true(
  multi_specificity$n_all_offtarget_products == 1L,
  "An off-target product on the intended genome was not rejected"
)
assert_true(
  multi_specificity$n_reduced_matches_missed == 1L,
  "The reduced matchProbePair result was treated as exhaustive"
)

plasmid_references <- make_test_references(
  expected_sequence,
  expected_sequence,
  expected_sequence
)
plasmid_specificity <- evaluate_pair_specificity(
  forward_probe,
  reverse_probe,
  plasmid_references,
  expected_product,
  specificity_config
)
offtarget_types <- unique(
  plasmid_specificity$amplicons$reference_type[
    plasmid_specificity$amplicons$off_target &
      !plasmid_specificity$amplicons$invalid_size
  ]
)
assert_true(
  all(c("target_plasmid", "cas_plasmid") %in% offtarget_types),
  "pTarget or pCas off-target amplicon was missed"
)

internal_mismatch_sequence <- replace_base(expected_sequence, 3L)
mismatch_specificity <- evaluate_pair_specificity(
  forward_probe,
  reverse_probe,
  make_test_references(
    internal_mismatch_sequence,
    strrep("G", 120L),
    strrep("C", 120L)
  ),
  expected_product,
  specificity_config
)
forward_site_with_mismatch <- mismatch_specificity$binding_sites[
  mismatch_specificity$binding_sites$primer_id == "forward" &
    mismatch_specificity$binding_sites$reference_type == "genome",
  ,
  drop = FALSE
]
assert_true(
  nrow(forward_site_with_mismatch) == 1L &&
    forward_site_with_mismatch$mismatches[[1]] == 1L &&
    forward_site_with_mismatch$mismatches_3p[[1]] == 0L,
  "An internal primer mismatch was not reported correctly"
)

critical_mismatch_sequence <- replace_base(
  expected_sequence,
  nchar(forward_probe)
)
critical_specificity <- evaluate_pair_specificity(
  forward_probe,
  reverse_probe,
  make_test_references(
    critical_mismatch_sequence,
    strrep("G", 120L),
    strrep("C", 120L)
  ),
  expected_product,
  specificity_config
)
assert_true(
  !critical_specificity$passed &&
    any(critical_specificity$binding_sites$mismatches_3p == 1L),
  "A critical 3-prime mismatch was allowed to form the intended product"
)

circular_bases <- rep("C", 120L)
circular_bases[101:120] <- strsplit(forward_probe, "", fixed = TRUE)[[1]]
circular_bases[21:40] <- strsplit(reverse_site, "", fixed = TRUE)[[1]]
circular_sequence <- paste(circular_bases, collapse = "")
circular_specificity <- evaluate_pair_specificity(
  forward_probe,
  reverse_probe,
  make_test_references(
    strrep("G", 160L),
    circular_sequence,
    strrep("A", 160L)
  ),
  list(
    reference_id = "target_plasmid::pTarget",
    start = 101L,
    end = 40L,
    size = 60L,
    allowed_products = 1L
  ),
  specificity_config
)
assert_true(
  circular_specificity$passed &&
    nrow(circular_specificity$amplicons) == 1L &&
    circular_specificity$amplicons$circular_wrap[[1]],
  "A circular wrap-around amplicon was not identified uniquely"
)

ranking_candidates <- data.frame(
  pair_id = c("first", "second", "third"),
  structure_passed = c(TRUE, TRUE, TRUE),
  specificity_passed = c(FALSE, TRUE, TRUE),
  openprimer_passed = c(TRUE, TRUE, TRUE),
  strict_qc_passed = c(FALSE, TRUE, TRUE),
  blocking_passed = c(FALSE, TRUE, TRUE),
  n_expected_product_deviations = c(0, 0, 0),
  n_high_risk_offtarget_products = c(0, 0, 0),
  n_all_offtarget_products = c(0, 0, 0),
  n_perfect_3p_offtarget_sites = c(0, 0, 0),
  openprimer_failed_soft_constraints = c(0, 0, 0),
  openprimer_penalty = c(0, 1, 1),
  max_dimer_risk = c(0, 1, 1),
  abs_tm_diff = c(0, 1, 1),
  primer3_pair_penalty = c(0, 1, 1),
  deleted_nt = c(0, 10, 10),
  primer3_index = c(1, 2, 3),
  stringsAsFactors = FALSE
)
ranking <- select_best_primer_pair(ranking_candidates)
assert_true(
  ranking$pair$pair_id[[1]] == "second",
  "The first rejected Primer3 pair prevented selection of the second"
)
assert_true(
  !ranking$ranking$selected[[1]],
  "A hard specificity failure was compensated by soft scores"
)
assert_true(
  select_best_primer_pair(ranking_candidates[2:3, ])$pair$pair_id[[1]] ==
    "second",
  "Primer ranking is not deterministic"
)

risky_specificity <- list(
  passed = FALSE,
  rejection_reason = "off_target_products=1",
  n_expected_products = 1L,
  n_high_risk_offtarget_products = 1L,
  n_all_offtarget_products = 1L,
  n_perfect_3p_offtarget_sites = 2L
)
optional_openprimer_failure <- list(
  passed = FALSE,
  rejection_reason = "openprimer_failed:EVAL_gc_clamp",
  abs_tm_diff = 2
)
lite_policy <- evaluate_filtering_policy(
  risky_specificity,
  optional_openprimer_failure,
  1L
)
default_policy <- evaluate_filtering_policy(
  risky_specificity,
  optional_openprimer_failure,
  2L
)
hard_policy <- evaluate_filtering_policy(
  risky_specificity,
  optional_openprimer_failure,
  3L
)
assert_true(
  lite_policy$blocking_passed && default_policy$blocking_passed &&
    !hard_policy$blocking_passed &&
    grepl("high_risk_off_target_products", hard_policy$blocking_reasons),
  "Filtering levels do not distinguish report, warning, and hard risks"
)

missing_specificity <- risky_specificity
missing_specificity$n_expected_products <- 0L
missing_openprimer <- unavailable_openprimer_result(
  "test",
  "expected product is unavailable"
)
assert_true(
  evaluate_filtering_policy(
    missing_specificity,
    missing_openprimer,
    1L
  )$blocking_passed &&
    !evaluate_filtering_policy(
      missing_specificity,
      missing_openprimer,
      2L
    )$blocking_passed,
  "Lite and default filtering do not differ on a missing intended product"
)

optional_specificity <- risky_specificity
optional_specificity$passed <- TRUE
optional_specificity$rejection_reason <- ""
optional_specificity$n_high_risk_offtarget_products <- 0L
optional_specificity$n_all_offtarget_products <- 0L
optional_specificity$n_perfect_3p_offtarget_sites <- 0L
hard_optional_policy <- evaluate_filtering_policy(
  optional_specificity,
  optional_openprimer_failure,
  3L
)
assert_true(
  hard_optional_policy$blocking_passed &&
    grepl("EVAL_gc_clamp", hard_optional_policy$warnings),
  "Hard filtering treated an optional openPrimeR criterion as blocking"
)
openprimer_unavailable_with_tm <- unavailable_openprimer_result(
  "test",
  "secondary-structure executable missing",
  3,
  "Primer3 fallback"
)
assert_true(
  evaluate_filtering_policy(
    optional_specificity,
    openprimer_unavailable_with_tm,
    3L
  )$blocking_passed,
  "Hard filtering rejected an available safe Primer3 Tm fallback"
)

temperature_risk <- optional_specificity
large_tm <- optional_openprimer_failure
large_tm$abs_tm_diff <- 6
assert_true(
  !evaluate_filtering_policy(temperature_risk, large_tm, 3L)$blocking_passed,
  "Hard filtering accepted a large Tm difference"
)

fallback_candidates <- ranking_candidates[2:3, ]
fallback_candidates$strict_qc_passed <- FALSE
fallback_candidates$blocking_passed <- TRUE
fallback_candidates$risk_warnings <- c("gc_clamp", "secondary_structure")
fallback_selection <- select_best_primer_pair(fallback_candidates)
assert_true(
  fallback_selection$pair$pair_id[[1]] == "second" &&
    fallback_selection$pair$fallback_selected[[1]],
  "Best non-blocking fallback candidate was not selected"
)

strict_preference <- fallback_candidates
strict_preference$strict_qc_passed[[2]] <- TRUE
strict_preference$openprimer_penalty <- c(0, 100)
strict_selection <- select_best_primer_pair(strict_preference)
assert_true(
  strict_selection$pair$pair_id[[1]] == "third" &&
    !strict_selection$pair$fallback_selected[[1]],
  "A better-ranked fallback displaced an available strict candidate"
)

intended_preference <- fallback_candidates
intended_preference$n_expected_product_deviations <- c(1, 0)
intended_selection <- select_best_primer_pair(intended_preference)
assert_true(
  intended_selection$pair$pair_id[[1]] == "third",
  "Lite fallback did not prefer a candidate with one intended product"
)

assert_error(
  assert_openprimer_constraints(
    c("primer_length", "gc_ratio"),
    c("primer_length", "secondary_structure")
  ),
  "secondary_structure"
)
sequence_forms <- openprimer_sequence_forms(
  "LF_LR",
  forward_probe,
  reverse_probe,
  paste0("SERVICE", forward_probe),
  paste0("TAIL", reverse_probe)
)
assert_true(
  sequence_forms$annealing_sequence[[1]] == forward_probe &&
    sequence_forms$full_oligo_sequence[[1]] != forward_probe,
  "Annealing and full oligo sequences were not kept separate"
)
assert_true(
  identical(unique(sequence_forms$reaction), "LF_LR") &&
    nrow(sequence_forms) == 2L,
  "Cross-dimerization input contains primers from another reaction"
)
physical_reactions <- bind_rows(lapply(
  c("LF_LR", "RF_RR", "scrF_scrR"),
  function(reaction) openprimer_sequence_forms(
    reaction,
    paste0(forward_probe, substr(reaction, 1L, 1L)),
    paste0(reverse_probe, substr(reaction, 2L, 2L)),
    paste0("FULLF", reaction),
    paste0("FULLR", reaction)
  )
))
assert_true(
  all(table(physical_reactions$reaction) == 2L) &&
    all(vapply(
      split(physical_reactions, physical_reactions$reaction),
      function(group) identical(group$primer_id, c("forward", "reverse")),
      logical(1)
    )),
  "Cross-dimerization groups mix primers from different PCR reactions"
)

formatted_failures <- format_openprimer_failures(
  data.frame(
    primer_length_fw = 24,
    primer_length_rev = 20,
    no_runs_fw = 5,
    no_runs_rev = 4,
    Structure_deltaG = -2.25,
    stringsAsFactors = FALSE
  ),
  c("EVAL_primer_length", "EVAL_no_runs", "EVAL_secondary_structure"),
  list(
    primer_length = c(min = 18, max = 22),
    no_runs = c(max = 4),
    secondary_structure = c(min = -1)
  )
)
assert_true(
  identical(formatted_failures, c(
    "EVAL_primer_length[primer_length_fw=24 нт (18–22 нт)]",
    "EVAL_no_runs[no_runs_fw=5 нт (≤ 4 нт)]",
    paste0(
      "EVAL_secondary_structure[Structure_deltaG=-2.25 ккал/моль ",
      "(≥ -1 ккал/моль)]"
    )
  )),
  "QC failures do not report measured and threshold values"
)
assert_true(
  identical(
    format_openprimer_failures(
      data.frame(unrelated_metric = NA_real_),
      "EVAL_unknown_constraint",
      list()
    ),
    "EVAL_unknown_constraint"
  ),
  "Unknown QC failures must preserve their machine-readable flag"
)

layout <- output_layout("results")
assert_true(
  identical(layout$wet_lab, file.path("results", "WetLab")),
  "WetLab root is incorrect"
)
assert_true(
  identical(layout$tech_report, file.path("results", "TechReport")),
  "TechReport root is incorrect"
)

buffer <- primer3_buffer_parameters()
assert_true(
  identical(
    unname(buffer),
    c(50, 1.5, 0.6, 50)
  ),
  "Primer3 buffer parameters are incorrect"
)
settings_path <- tempfile("primer3-settings-", fileext = ".txt")
write_primer3_settings(settings_path, 400L, 500L, buffer)
settings_text <- readLines(settings_path)
assert_true(
  any(settings_text == "PRIMER_SALT_MONOVALENT=50"),
  "Monovalent salt concentration is absent from Primer3 settings"
)
assert_true(
  any(settings_text == "PRIMER_SALT_DIVALENT=1.5"),
  "Divalent salt concentration is absent from Primer3 settings"
)
assert_true(
  any(settings_text == "PRIMER_DNTP_CONC=0.6"),
  "dNTP concentration is absent from Primer3 settings"
)
assert_true(
  any(settings_text == "PRIMER_DNA_CONC=50"),
  "DNA concentration is absent from Primer3 settings"
)
assert_true(
  all(c(
    "PRIMER_MIN_SIZE=18",
    "PRIMER_OPT_SIZE=21",
    "PRIMER_MAX_SIZE=22",
    "PRIMER_MAX_POLY_X=4",
    "PRIMER_PAIR_MAX_DIFF_TM=5.0",
    "PRIMER_GC_CLAMP=1"
  ) %in% settings_text),
  "Primer3 generation limits do not match high-stringency QC"
)

screening_sizes <- calculate_screening_product_sizes(
  data.frame(PRIMER_PAIR_PRODUCT_SIZE = 500L),
  DNAString(paste(rep("A", 1000L), collapse = "")),
  DNAStringSet(c(
    edited_genome = paste(rep("A", 850L), collapse = "")
  ))
)
assert_true(
  identical(
    screening_sizes,
    c(unsuccessful_insertion_bp = 500L, successful_insertion_bp = 350L)
  ),
  "Screening PCR sizes do not reflect the edited-genome length change"
)

n20_distances <- calculate_n20_arm_distances(
  data.frame(
    target_sequence = c("ACGTACGTACGTACGTACGTAGG", "TGCATGCATGCATGCATGCAAGG"),
    strand = c("+", "-"),
    n20_start = c(141L, 161L),
    n20_end = c(160L, 180L),
    stringsAsFactors = FALSE
  ),
  c(50L, 100L, 221L, 270L),
  "+"
)
assert_true(
  identical(n20_distances$left_arm_distance_bp, c(40L, 60L)) &&
    identical(n20_distances$right_arm_distance_bp, c(60L, 40L)),
  "Per-N20 distances to homology arms are incorrect"
)
minus_n20_distances <- calculate_n20_arm_distances(
  data.frame(
    target_sequence = "ACGTACGTACGTACGTACGTAGG",
    strand = "-",
    n20_start = 141L,
    n20_end = 160L,
    stringsAsFactors = FALSE
  ),
  c(50L, 100L, 221L, 270L),
  "-"
)
assert_true(
  minus_n20_distances$left_arm_distance_bp[[1]] == 60L &&
    minus_n20_distances$right_arm_distance_bp[[1]] == 40L,
  "Minus-strand arm labels were not oriented to the target"
)

cassette_sequence <- paste0(
  "ACGTACGTACGTACGTACGT",
  SGRNA_SCAFFOLD,
  "GAATTCTCTAGAGTCGAC"
)
derived_annealing <- derive_sgrna_annealing(cassette_sequence)
assert_true(
  derived_annealing$original_n20 == "ACGTACGTACGTACGTACGT" &&
    derived_annealing$scaffold == SGRNA_SCAFFOLD &&
    derived_annealing$forward == SGRNA_FORWARD_ANNEALING &&
    derived_annealing$reverse == "AAAAAAAGCACCGACTCGGTGCC",
  "The synthetic pTarget architecture or annealing sites were not recovered"
)
assert_true(
  paste0(SGRNA_REVERSE_OVERHANG, derived_annealing$reverse) ==
    "AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC",
  "The generated reverse sgRNA primer differs from the expected sequence"
)
assert_error(
  derive_sgrna_annealing("ACGT"),
  "несовместима со схемой sgRNA"
)
synthetic_ptarget <- DNAString(paste0(
  "GCAGGGGACTAGT",
  cassette_sequence,
  "CTGCAG",
  strrep("A", 100L)
))
synthetic_sgrna_template <- locate_circular_pcr_template(
  synthetic_ptarget,
  derived_annealing$forward,
  derived_annealing$reverse,
  length(synthetic_ptarget)
)
assert_true(
  synthetic_sgrna_template$strand == "+" &&
    as.character(synthetic_sgrna_template$sequence) == SGRNA_SCAFFOLD,
  "The synthetic primers do not define one exact scaffold PCR product"
)
circular_ptarget <- DNAString(paste0("GCAGGGGACTAGT", cassette_sequence, "CT"))
assembly_bridge <- "ATGACTGCCCGCAAG"
assembly_n20 <- c(
  "ACGTACGTACGTACGTACGT",
  "TGCATGCATGCATGCATGCA"
)
assembly_sgrna_products <- paste0(
  "ACGACTAGT",
  assembly_n20,
  "GGGAGCGTCAACT"
)
assembly_left_product <- paste0(
  "AGCGTCAACTAACCGG",
  assembly_bridge
)
assembly_right_product <- paste0(
  assembly_bridge,
  "TTGGCCCTGCAGCGT"
)
ptarget_model <- model_edited_ptargets(
  circular_ptarget,
  assembly_sgrna_products,
  assembly_left_product,
  assembly_bridge,
  assembly_right_product,
  name_prefix = "geneA",
  cassette_arc = "forward"
)
assert_true(
  ptarget_model$restriction_pair$orientation == "+" &&
    ptarget_model$restriction_pair$site2_start == 14L + nchar(cassette_sequence),
  "Circular restriction site crossing the FASTA origin was not found"
)
assert_true(
  length(ptarget_model$sequences) == 2L &&
    identical(
      names(ptarget_model$sequences),
      c("geneA_pTarget_N20_1", "geneA_pTarget_N20_2")
    ),
  "A separate edited pTarget was not made for every N20"
)
assert_true(
  all(startsWith(as.character(ptarget_model$sequences), "ACTAGT")) &&
    !any(startsWith(as.character(ptarget_model$sequences), "AACTAGT")),
  "Edited pTarget retains the legacy duplicated-A error"
)
minus_oriented <- find_oriented_restriction_pair(
  DNAString(reverse_complement_string("AAAACTTTTGCGTAGGG")),
  "AAAAC",
  "GCGTA"
)
assert_true(
  minus_oriented$orientation == "-",
  "Reverse-oriented restriction pair was not recognized"
)
assert_error(
  find_oriented_restriction_pair(
    DNAString("ACTAGTAAAAACTAGTCCCCCTGCAGGGG"),
    "ACTAGT",
    "CTGCAG"
  ),
  "ровно один раз"
)

pcr_forward_annealing <- "ACGTTGCAAGTCGATCGTAC"
pcr_reverse_annealing <- "TGACCTGACTAGGCTACGTA"
pcr_full_forward <- paste0("AGCGTCAACT", pcr_forward_annealing)
pcr_full_reverse <- paste0("ATGACTGCCCGCAAG", pcr_reverse_annealing)
pcr_template <- paste0(
  pcr_forward_annealing,
  "AACCGGTTAACC",
  reverse_complement_string(pcr_reverse_annealing)
)
pcr_products <- simulate_full_primer_pcr(
  "test_full_primer_product",
  "Full-primer PCR test",
  "fixture:1-52 (+)",
  pcr_template,
  "test_F",
  "test_R",
  pcr_full_forward,
  pcr_full_reverse,
  pcr_forward_annealing,
  pcr_reverse_annealing,
  57,
  primer3_buffer_parameters(),
  200L
)
assert_true(
  startsWith(pcr_products$sequence, pcr_full_forward) &&
    endsWith(
      pcr_products$sequence,
      reverse_complement_string(pcr_full_reverse)
    ),
  "DECIPHER PCR product does not include full service-tailed primers"
)

openprimer_report_metrics <- format_openprimer_report_metrics(data.frame(
  constraints_passed = TRUE,
  gc_ratio_fw = 0.5,
  Tm_C_fw = 62.3456,
  Cross_Dimer_DeltaG = -1.25,
  EVAL_gc_ratio = TRUE,
  penalty = 1.5,
  stringsAsFactors = FALSE
))
assert_true(
  any(
    openprimer_report_metrics$metric == "GC-состав forward-праймера, %" &
      openprimer_report_metrics$value == "50"
  ) &&
    any(
      openprimer_report_metrics$metric == "Проверка GC-состава" &
        openprimer_report_metrics$value == "пройдено"
    ),
  "openPrimeR metrics were not formatted with readable names and values"
)

wet_lab_dir <- tempfile("2pac-wet-lab-")
wet_lab_sequences <- DNAStringSet(c(
  geneA_LF = "ACGTACGT",
  geneA_scrF = "TGCATGCA"
))
write_wet_lab_outputs(
  wet_lab_dir,
  list(query_name = "geneA"),
  "cds",
  wet_lab_sequences,
  c("left_arm_forward_primer", "screening_forward_primer"),
  data.frame(
    name = names(wet_lab_sequences),
    purpose = c("left_arm_forward_primer", "screening_forward_primer"),
    annealing_sequence = as.character(wet_lab_sequences),
    tm_c = c(61.2, 63.4),
    stringsAsFactors = FALSE
  ),
  screening_sizes,
  n20_distances,
  list(
    pair_id = "screening_02",
    filtering_level = 2L,
    filtering_mode = "default",
    selection_status = "selected_with_warnings",
    fallback_used = TRUE,
    warnings = "secondary_structure",
    offtarget_products = 0L,
    high_risk_offtarget_products = 0L,
    perfect_3p_offtarget_sites = 0L,
    openprimer_metrics = openprimer_report_metrics
  ),
  DNAStringSet(c(edited_genome = "AACCGGTTAACC")),
  ptarget_model$sequences,
  ptarget_model$restriction_pair,
  pcr_products
)
assert_true(
  all(c(
    "final_sequences.fasta",
    "final_sequences.txt",
    "wet_lab_report.txt",
    "edited_genome.fasta",
    "edited_pTargets.fasta",
    "pcr_products.fasta",
    "pcr_products.tsv"
  ) %in% list.files(wet_lab_dir)),
  "WetLab lacks a required file"
)
assert_true(
  !any(c("design.log", "n20_table.tsv", "report.tsv") %in%
    list.files(wet_lab_dir)),
  "WetLab contains technical pipeline output"
)
assert_true(
  length(readDNAStringSet(file.path(wet_lab_dir, "final_sequences.fasta"))) == 2L,
  "WetLab FASTA does not contain the complete final sequence set"
)
wet_lab_table <- read_tsv(
  file.path(wet_lab_dir, "final_sequences.txt"),
  show_col_types = FALSE
)
assert_true(
  identical(wet_lab_table$name, names(wet_lab_sequences)),
  "WetLab TXT does not contain the complete final sequence set"
)
assert_true(
  identical(wet_lab_table$sequence, unname(as.character(wet_lab_sequences))),
  "WetLab FASTA and TXT contain different sequences"
)
assert_true(
  length(readDNAStringSet(file.path(wet_lab_dir, "edited_pTargets.fasta"))) ==
    nrow(n20_distances),
  "WetLab does not contain one edited pTarget per N20"
)
wet_lab_pcr <- read_tsv(
  file.path(wet_lab_dir, "pcr_products.tsv"),
  show_col_types = FALSE
)
assert_true(
  identical(wet_lab_pcr$name, pcr_products$name) &&
    identical(wet_lab_pcr$sequence, pcr_products$sequence),
  "WetLab PCR table changed the modelled product identity or sequence"
)
wet_lab_report <- readLines(
  file.path(wet_lab_dir, "wet_lab_report.txt"),
  encoding = "UTF-8"
)
assert_true(
  any(grepl("Без успешного нокаута.*500", wet_lab_report)),
  "WetLab report lacks the unsuccessful-knockout PCR size"
)
assert_true(
  any(grepl("С успешным нокаутом.*350", wet_lab_report)),
  "WetLab report lacks the successful-knockout PCR size"
)
assert_true(
  any(grepl("geneA_LF.*61.2", wet_lab_report)),
  "WetLab report lacks a primer annealing temperature"
)
assert_true(
  any(grepl("N20_1.*40.*60", wet_lab_report)),
  "WetLab report lacks per-N20 distances to both homology arms"
)
assert_true(
  any(grepl("Оффтаргетные ПЦР-продукты, всего.*0", wet_lab_report)) &&
    any(grepl("GC-состав forward-праймера.*50", wet_lab_report)),
  "WetLab report lacks readable screening primer QC"
)
assert_true(
  any(grepl("Режим фильтрации праймеров.*2 \\(default\\)", wet_lab_report)) &&
    any(grepl("Предупреждения primer QC.*secondary_structure", wet_lab_report)),
  "WetLab report lacks filtering mode or fallback warnings"
)
assert_true(
  any(grepl("edited_genome.fasta.*edited_genome", wet_lab_report)) &&
    sum(grepl("edited_pTargets.fasta.*edited_pTarget", wet_lab_report)) == 2L,
  "WetLab report lacks the modelled genome or per-N20 pTargets"
)
assert_true(
  any(grepl("DECIPHER::AmplifyDNA", wet_lab_report, fixed = TRUE)) &&
    any(grepl("test_full_primer_product.*fixture:1-52.*57 C", wet_lab_report)),
  "WetLab report lacks PCR product location, conditions, or sequence"
)

error_root <- tempfile("2pac-error-")
dir.create(error_root)
stale_wet_lab_dir <- file.path(
  error_root,
  "WetLab",
  "missing_gene_results"
)
dir.create(stale_wet_lab_dir, recursive = TRUE)
stale_wet_lab_files <- file.path(
  stale_wet_lab_dir,
  c(
    "final_sequences.fasta",
    "final_sequences.txt",
    "wet_lab_report.txt",
    "edited_genome.fasta",
    "edited_pTargets.fasta",
    "pcr_products.fasta",
    "pcr_products.tsv"
  )
)
file.create(stale_wet_lab_files)
error_input <- list(
  output_dir = error_root,
  annotation = data.frame(
    locus_tag = character(),
    gene = character()
  )
)
target_error <- tryCatch(
  {
    design_target(error_input, "genome", "missing_gene", "cds")
    NULL
  },
  error = identity
)
assert_true(
  inherits(target_error, "target_design_error"),
  "Target error was not classified"
)
target_error_dir <- file.path(
  error_root,
  "TechReport",
  "missing_gene_results"
)
assert_true(
  file.exists(file.path(target_error_dir, "design.log")),
  "design.log was not created"
)
assert_true(
  file.exists(file.path(target_error_dir, "error.txt")),
  "error.txt was not created"
)
error_text <- readLines(file.path(target_error_dir, "error.txt"))
assert_true(any(grepl("stage\\tfeature_lookup", error_text)), "Stage is absent")
assert_true(any(grepl("gene\\tmissing_gene", error_text)), "Gene is absent")
assert_true(
  !any(file.exists(stale_wet_lab_files)),
  "A failed target left stale WetLab output"
)

trace_dir <- tempfile("2pac-primer-trace-")
dir.create(trace_dir)
trace <- new_primer_qc_trace()
trace$openprimer[[1]] <- data.frame(
  reaction = "LF_LR",
  pair_id = "pair_1",
  Structure = "line one\nline two",
  EVAL_secondary_structure = TRUE,
  constraints_passed = TRUE,
  penalty = 0,
  stringsAsFactors = FALSE
)
write_primer_qc_trace(trace, trace_dir)
trace_lines <- readLines(file.path(trace_dir, "primer_openprimer_qc.tsv"))
assert_true(
  length(trace_lines) == 2L && grepl("line one line two", trace_lines[[2]]),
  "Multiline openPrimeR metrics break the TSV row contract"
)

message("Unit tests passed")
