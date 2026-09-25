#!/usr/bin/env Rscript

source("oligo_designer.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) {
    stop(message, call. = FALSE)
  }
}

test_screening_fixture <- function(strand, retry = FALSE, fallback = FALSE) {
  target_dir <- tempfile("2pac-screening-fixture-")
  dir.create(target_dir)
  on.exit(unlink(target_dir, recursive = TRUE), add = TRUE)

  original_call_primer3 <- callPrimer3
  original_evaluate_candidate <- evaluate_candidate_reaction
  on.exit(
    assign("callPrimer3", original_call_primer3, envir = .GlobalEnv),
    add = TRUE
  )
  on.exit(
    assign(
      "evaluate_candidate_reaction",
      original_evaluate_candidate,
      envir = .GlobalEnv
    ),
    add = TRUE
  )

  assign(
    "callPrimer3",
    function(...) {
      data.frame(
        PRIMER_LEFT_SEQUENCE = c("ACGTCGATCGTAGCTACGTA", "GCTAGTCGATGCTACGTAGC"),
        PRIMER_RIGHT_SEQUENCE = c("TGCATCGATGCTAGTCGTAC", "CGTACGATCGTAGCATCGAC"),
        PRIMER_LEFT_pos = c(1L, 10L),
        PRIMER_RIGHT_pos = c(1200L, 1250L),
        PRIMER_LEFT_TM = c(62.5, 62.7),
        PRIMER_RIGHT_TM = c(62.6, 62.8),
        PRIMER_PAIR_PRODUCT_SIZE = c(1200L, 1241L),
        PRIMER_PAIR_PENALTY = c(0.1, 0.2),
        stringsAsFactors = FALSE
      )
    },
    envir = .GlobalEnv
  )

  assign(
    "evaluate_candidate_reaction",
    function(input, primer_row, full_forward, full_reverse, reaction, pair_id,
             trace) {
      passed <- endsWith(pair_id, "_screening_02")
      specificity <- list(
        passed = passed,
        binding_sites = data.frame(
          primer_id = c("forward", "reverse"),
          reference_id = input$genome_reference_id,
          reference_type = "genome",
          start = c(primer_row$genome_start[[1]], primer_row$genome_end[[1]] - 19L),
          end = c(primer_row$genome_start[[1]] + 19L, primer_row$genome_end[[1]]),
          strand = c("+", "-"),
          mismatches = 0L,
          mismatch_positions = "",
          mismatches_3p = 0L,
          stringsAsFactors = FALSE
        ),
        amplicons = data.frame(
          reference_id = input$genome_reference_id,
          reference_type = "genome",
          start = primer_row$genome_start[[1]],
          end = primer_row$genome_end[[1]],
          product_size = primer_row$PRIMER_PAIR_PRODUCT_SIZE[[1]],
          sequence = paste(
            rep("A", primer_row$PRIMER_PAIR_PRODUCT_SIZE[[1]]),
            collapse = ""
          ),
          intended = passed,
          off_target = !passed,
          invalid_size = FALSE,
          circular_wrap = FALSE,
          rejection_reason = if (passed) "" else "off_target",
          stringsAsFactors = FALSE
        ),
        n_high_risk_offtarget_products = as.integer(!passed),
        n_all_offtarget_products = as.integer(!passed),
        n_perfect_3p_offtarget_sites = as.integer(!passed),
        n_expected_products = as.integer(passed),
        rejection_reason = if (passed) "" else "off_target_products=1"
      )
      openprimer <- if (passed) {
        list(
          passed = !fallback,
          metrics = data.frame(
            reaction = reaction,
            constraints_passed = !fallback,
            penalty = 1,
            stringsAsFactors = FALSE
          ),
          failed_soft_constraints = 0L,
          penalty = 1,
          max_dimer_risk = 0,
          abs_tm_diff = 0.1,
          rejection_reason = if (fallback) {
            paste0(
              "openprimer_failed:EVAL_primer_length[",
              "primer_length_fw=24 нт (18–22 нт)]"
            )
          } else ""
        )
      } else NULL
      if (is.null(openprimer)) {
        openprimer <- unavailable_openprimer_result(
          reaction,
          "expected product is unavailable"
        )
      }
      policy <- evaluate_filtering_policy(
        specificity,
        openprimer,
        input$parameters$filtering_level,
        input$parameters$primer_qc
      )
      append_primer_qc_trace(
        trace,
        reaction,
        pair_id,
        specificity = specificity,
        openprimer = openprimer
      )
      list(
        passed = policy$blocking_passed,
        specificity = specificity,
        openprimer = openprimer,
        policy = policy,
        rejection_reason = policy$blocking_reasons,
        risk_warnings = policy$warnings
      )
    },
    envir = .GlobalEnv
  )

  trace <- new_primer_qc_trace()
  trace$ranking[[1]] <- data.frame(
    pair_id = "homology_fixture",
    reaction = "homology_combination",
    primer3_index = 1L,
    structure_passed = TRUE,
    specificity_passed = TRUE,
    openprimer_passed = TRUE,
    strict_qc_passed = TRUE,
    blocking_passed = TRUE,
    n_expected_product_deviations = 0L,
    n_high_risk_offtarget_products = 0L,
    n_all_offtarget_products = 0L,
    n_perfect_3p_offtarget_sites = 0L,
    openprimer_failed_soft_constraints = 0L,
    openprimer_penalty = 1,
    max_dimer_risk = 0,
    abs_tm_diff = 0.1,
    primer3_pair_penalty = 0.1,
    deleted_nt = 200L,
    selected = TRUE,
    fallback_selected = FALSE,
    risk_warnings = "",
    rejection_reason = "",
    stringsAsFactors = FALSE
  )

  pair <- data.frame(
    PRIMER_LEFT_SEQUENCE = c("ATGCGTACGATCGTACGTAG", "GATCGTAGCTAGTCGATGCA"),
    PRIMER_RIGHT_SEQUENCE = c("CTACGTACGATCGTACGCAT", "TGCATCGACTAGCTACGATC"),
    PRIMER_LEFT_pos = c(1L, 1L),
    PRIMER_RIGHT_pos = c(300L, 400L),
    PRIMER_LEFT_len = c(20L, 20L),
    PRIMER_RIGHT_len = c(20L, 20L),
    PRIMER_LEFT_TM = c(60, 60),
    PRIMER_RIGHT_TM = c(60, 60),
    genome_start = c(201L, 701L),
    genome_end = c(500L, 1100L),
    stringsAsFactors = FALSE
  )
  make_amplicon <- function(forward, reverse, width, fill) {
    paste0(
      forward,
      paste(rep(fill, width - nchar(forward) - nchar(reverse)), collapse = ""),
      reverse_complement_string(reverse)
    )
  }
  left_template <- make_amplicon(
    pair$PRIMER_LEFT_SEQUENCE[[1]],
    pair$PRIMER_RIGHT_SEQUENCE[[1]],
    300L,
    "A"
  )
  right_template <- make_amplicon(
    pair$PRIMER_LEFT_SEQUENCE[[2]],
    pair$PRIMER_RIGHT_SEQUENCE[[2]],
    400L,
    "C"
  )
  genome_sequence <- paste(rep("T", 2000L), collapse = "")
  substr(genome_sequence, 201L, 500L) <- left_template
  substr(genome_sequence, 701L, 1100L) <- right_template
  substr(genome_sequence, 10L, 30L) <- "GCTAGTCGATGCTACGTAGC"
  substr(genome_sequence, 1231L, 1250L) <- reverse_complement_string(
    "CGTACGATCGTAGCATCGAC"
  )
  ptarget_cassette <- paste0(
    "ACGTACGTACGTACGTACGT",
    SGRNA_SCAFFOLD,
    "GAATTCTCTAGAGTCGAC"
  )
  ptarget_annealing <- derive_sgrna_annealing(ptarget_cassette)
  input <- list(
    genome = DNAString(genome_sequence),
    genome_reference_id = "fixture_genome",
    genome_contig = "fixture_genome",
    target_plasmid_sequence = DNAString(paste0(
      "GGGACTAGT",
      ptarget_cassette,
      "CTGCAG",
      strrep("A", 200L)
    )),
    target_plasmid_name = "fixture_pTarget",
    parameters = list(
      filtering_level = 2L,
      cds_fs = FALSE,
      ncrna_fs = FALSE,
      n20_offtarget = 0L,
      n20_arm_min_distance = 40L,
      site1 = "ACTAGT",
      site2 = "CTGCAG",
      ptarget_cassette_arc = "shortest",
      ptarget_cassette_length = nchar(ptarget_cassette),
      ptarget_original_n20 = ptarget_annealing$original_n20,
      sgrna_scaffold = ptarget_annealing$scaffold,
      sgrna_forward_annealing = ptarget_annealing$forward,
      sgrna_reverse_annealing = ptarget_annealing$reverse,
      sgrna_annealing_temp_c = 60,
      primer3_buffer = primer3_buffer_parameters(),
      primer_qc = list(max_product_size = 2000L)
    ),
    tools = list(primer3 = "unused", primer3_config = "unused")
  )
  if (strand == "-") {
    pair <- pair[2:1, , drop = FALSE]
    forward <- pair$PRIMER_LEFT_SEQUENCE
    pair$PRIMER_LEFT_SEQUENCE <- pair$PRIMER_RIGHT_SEQUENCE
    pair$PRIMER_RIGHT_SEQUENCE <- forward
    left_template <- reverse_complement_string(right_template)
    right_template <- reverse_complement_string(substr(genome_sequence, 201L, 500L))
  }
  arms <- list(
    pair = pair,
    ticks = c(201L, 500L, 701L, 1100L),
    left = DNAString(paste0(left_template, paste(rep("G", 200L), collapse = ""))),
    right = DNAString(paste0(right_template, paste(rep("G", 100L), collapse = ""))),
    selected_pair_id = "homology_fixture",
    primer_qc_trace = trace
  )
  selected <- list(table = data.frame(
    target_sequence = "ACGTACGTACGTACGTACGTAGG",
    strand = "+",
    n20_start = 560L,
    n20_end = 579L,
    stringsAsFactors = FALSE
  ))
  feature <- list(
    display_name = "fixture",
    query_name = "fixture",
    strand = strand
  )
  log_path <- file.path(target_dir, "design.log")

  if (retry) {
    original_plasmid <- input$target_plasmid_sequence
    input$target_plasmid_sequence <- DNAString("ACTAGTCCCCCTGCAGGGG")
    failed <- tryCatch(
      write_design_outputs(input, feature, selected, arms, "cds", target_dir, log_path),
      error = identity
    )
    assert_true(inherits(failed, "error") && grepl("sgRNA|pTarget", conditionMessage(failed)),
                "The first attempt must fail after screening QC")
    assert_true(!any(bind_rows(trace$ranking)$selected),
                "A rejected output attempt still marks its pairs as final")
    input$target_plasmid_sequence <- original_plasmid
    # The next homology search supplies its selected row, retaining prior trace.
    next_homology <- trace$ranking[[1]]
    next_homology$pair_id <- "homology_retry"
    next_homology$selected <- TRUE
    trace$ranking[[length(trace$ranking) + 1L]] <- next_homology
    arms$selected_pair_id <- "homology_retry"
  }
  output_id <- paste0(arms$selected_pair_id, "_outputs_", ifelse(retry, 2L, 1L))
  expected_pair_id <- paste0(output_id, "_screening_02")
  emitted_warnings <- character()
  result <- withCallingHandlers(
    write_design_outputs(
      input,
      feature,
      selected,
      arms,
      "cds",
      target_dir,
      log_path
    ),
    warning = function(condition) {
      emitted_warnings <<- c(emitted_warnings, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  write_primer_qc_trace(result$primer_qc_trace, target_dir)
  wet_lab_dir <- file.path(target_dir, "WetLab")
  write_wet_lab_outputs(
    wet_lab_dir,
    feature,
    "cds",
    result$wet_lab$sequences,
    result$wet_lab$sequence_purposes,
    result$wet_lab$primer_metrics,
    result$wet_lab$screening_product_sizes,
    result$wet_lab$n20_distances,
    result$wet_lab$screening_qc,
    result$wet_lab$edited_genome,
    result$wet_lab$edited_ptargets,
    result$wet_lab$ptarget_site_pair,
    result$wet_lab$pcr_products
  )

  # Independently specified deletion: 501..700. Keep the current CDS bridge
  # (13 nt); its frame rule is a separate issue from genome assembly.
  genomic_bridge <- if (strand == "+") "ATGACTGCCCGCA" else "TGCGGGCAGTCAT"
  expected_genome <- paste0(
    substr(genome_sequence, 1L, 500L),
    genomic_bridge,
    substr(genome_sequence, 701L, 2000L)
  )
  for (output_dir in c(target_dir, wet_lab_dir)) {
    edited <- readDNAStringSet(file.path(output_dir, "edited_genome.fasta"))
    assert_true(
      identical(names(edited), "edited_genome") &&
        identical(unname(as.character(edited)), expected_genome) &&
        width(edited)[[1]] == 1813L,
      sprintf("Edited genome (%s) must preserve all bases outside 501..700", strand)
    )
  }
  assert_true(
    identical(unname(as.integer(result$wet_lab$screening_product_sizes)),
              c(1241L, 1054L)),
    "Screening sizes must reflect only the 200-nt deletion and 13-nt insertion"
  )
  expected_screening <- substr(expected_genome, 10L, 1063L)
  assert_true(
    identical(result$wet_lab$pcr_products$sequence[
      result$wet_lab$pcr_products$name == "screening_edited_genome"
    ], expected_screening),
    "Edited screening PCR differs from the independently specified allele"
  )
  donor <- paste0(left_template, "ATGACTGCCCGCA", right_template)
  genomic_donor <- if (strand == "+") donor else reverse_complement_string(donor)
  assert_true(
    identical(substr(expected_genome, 201L, 913L), genomic_donor) &&
      all(grepl(donor, as.character(result$wet_lab$edited_ptargets), fixed = TRUE)),
    "Edited allele and pTarget must contain the same complete donor junction"
  )

  assert_true(
    identical(result$screening_pair_id, expected_pair_id),
    "The second screening pair was not selected"
  )
  screening <- as.character(readDNAStringSet(result$screening_path))
  assert_true(
    identical(
      unname(screening),
      c("GCTAGTCGATGCTACGTAGC", "CGTACGATCGTAGCATCGAC")
    ),
    "Screening output does not contain the selected Primer3 row"
  )
  ranking <- read.delim(
    file.path(target_dir, "primer_pair_ranking.tsv"),
    check.names = FALSE
  )
  screening_ranking <- ranking[ranking$reaction == "scrF_scrR" & startsWith(ranking$pair_id, output_id), , drop = FALSE]
  assert_true(
    nrow(screening_ranking) == 2L &&
      !screening_ranking$selected[[1]] &&
      screening_ranking$selected[[2]],
    "Screening ranking did not reject row 1 and select row 2"
  )
  amplicons <- read.delim(
    file.path(target_dir, "primer_amplicons.tsv"),
    check.names = FALSE
  )
  assert_true(
    any(
      amplicons$pair_id == expected_pair_id &
        amplicons$reaction == "scrF_scrR" &
        amplicons$intended &
        !amplicons$off_target
    ),
    "Selected screening pair lacks its intended amplicon"
  )
  log_lines <- readLines(log_path)
  assert_true(
    any(grepl(paste0("primer_qc\\tTRY\\tpair_id=", output_id, "_screening_01"), log_lines)) &&
      any(grepl(paste0("primer_qc\\tREJECTED\\tpair_id=", output_id, "_screening_01"), log_lines)) &&
      any(grepl(paste0("primer_qc\\tOK\\tpair_id=", expected_pair_id), log_lines)),
    "Screening TRY/REJECTED/OK events are incomplete"
  )
  report <- readLines(file.path(target_dir, "report.tsv"))
  assert_true(
    any(report == paste0("screening_pair_id\t", expected_pair_id)),
    "report.tsv lacks the selected screening pair ID"
  )
  assert_true(
    all(c("screening_unsuccessful_insertion_bp\t1241",
          "screening_successful_insertion_bp\t1054") %in% report),
    "TechReport must contain the independently expected screening sizes"
  )
  assert_true(
    any(report == paste0(
      "primer_qc_fallback_used\t",
      ifelse(fallback, "TRUE", "FALSE")
    )) &&
      (!fallback || any(grepl(
        "primer_qc_warnings.*primer_length_fw=24 нт.*18–22 нт",
        report
      ))),
    "TechReport lacks the primer-QC fallback status or warnings"
  )
  assert_true(
    result$wet_lab$n20_distances$left_arm_distance_bp[[1]] ==
      ifelse(strand == "+", 59L, 121L) &&
      result$wet_lab$n20_distances$right_arm_distance_bp[[1]] ==
        ifelse(strand == "+", 121L, 59L),
    "WetLab data lacks the selected N20-to-arm distances"
  )
  assert_true(
    result$wet_lab$screening_qc$offtarget_products == 0L &&
      any(result$wet_lab$screening_qc$openprimer_metrics$metric ==
        "Все обязательные ограничения пройдены"),
    "WetLab data lacks screening QC"
  )
  assert_true(
    nrow(result$wet_lab$pcr_products) == 5L &&
      all(c(
        "sgRNA_N20_1",
        "left_homology_arm",
        "right_homology_arm",
        "screening_original_genome",
        "screening_edited_genome"
      ) %in% result$wet_lab$pcr_products$name) &&
      identical(
        result$wet_lab$pcr_products$length_bp,
        nchar(result$wet_lab$pcr_products$sequence)
      ),
    "The complete PCR-product set was not modelled"
  )
  sgrna_product <- result$wet_lab$pcr_products$sequence[
    result$wet_lab$pcr_products$name == "sgRNA_N20_1"
  ]
  assert_true(
    length(sgrna_product) == 1L &&
      endsWith(sgrna_product, SGRNA_PRODUCT_OVERLAP),
    "sgRNA PCR product lacks the reverse-primer assembly tail"
  )
  screening_products <- result$wet_lab$pcr_products[
    startsWith(result$wet_lab$pcr_products$name, "screening_"),
    ,
    drop = FALSE
  ]
  assert_true(
    identical(
      screening_products$length_bp,
      unname(as.integer(result$wet_lab$screening_product_sizes))
    ),
    "DECIPHER screening products disagree with reported product sizes"
  )
  left_product <- result$wet_lab$pcr_products$sequence[
    result$wet_lab$pcr_products$name == "left_homology_arm"
  ]
  full_left_primer <- as.character(result$wet_lab$sequences[
    grepl("_LF$", names(result$wet_lab$sequences))
  ])
  assert_true(
    length(full_left_primer) == 1L &&
      startsWith(left_product, full_left_primer),
    "Homology-arm PCR product lacks its full service-tailed primer"
  )
  wet_lab_report <- readLines(
    file.path(wet_lab_dir, "wet_lab_report.txt"),
    encoding = "UTF-8"
  )
  assert_true(
    all(c("Без успешного нокаута (исходный аллель), п.н.\t1241",
          "С успешным нокаутом (редактированный аллель), п.н.\t1054") %in%
        wet_lab_report),
    "WetLab report must contain the independently expected screening sizes"
  )
  assert_true(
    any(grepl(ifelse(strand == "+", "N20_1.*59.*121", "N20_1.*121.*59"),
              wet_lab_report)) &&
      any(grepl("Оффтаргетные ПЦР-продукты, всего.*0", wet_lab_report)) &&
      any(grepl("Все обязательные ограничения пройдены.*пройдено", wet_lab_report)) &&
      any(grepl("DECIPHER::AmplifyDNA", wet_lab_report, fixed = TRUE)),
    "Selected screening QC was not written to the WetLab report"
  )
  assert_true(
    any(grepl(
      paste0("QC fallback.*", ifelse(fallback, "TRUE", "FALSE")),
      wet_lab_report
    )) &&
      (!fallback || any(grepl(
        "Предупреждения primer QC.*primer_length_fw=24 нт.*18–22 нт",
        wet_lab_report
      ))),
    "WetLab report lacks the primer-QC fallback status or warnings"
  )
  if (fallback) {
    assert_true(
      any(grepl(
        "primer_qc\tWARNING\t.*primer_length_fw=24 нт.*18–22 нт",
        log_lines
      )) &&
        any(grepl(
          "QC-рисками.*primer_length_fw=24 нт.*18–22 нт",
          emitted_warnings
        )),
      "Fallback warning is missing from design.log or command output"
    )
  }
}

test_screening_fixture("+")
test_screening_fixture("-")
test_screening_fixture("+", retry = TRUE)
test_screening_fixture("+", fallback = TRUE)

message("Screening integration fixture passed")
