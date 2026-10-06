#!/usr/bin/env Rscript
source("oligo_designer.R")

local({
  # Use real settings and policy; replace only external tool execution.
  scope <- new.env(parent = globalenv())
  scope$configure_openprimer_environment <- function() invisible(NULL)
  scope$validate_openprimer_tools <- function(...) invisible(NULL)
  # read_settings may disable external-tool constraints on this test host.
  scope$assert_openprimer_constraints <- function(...) invisible(NULL)
  loader <- load_openprimer_settings
  environment(loader) <- scope
  loaded <- loader()
  stopifnot(identical(unname(openPrimeR::constraints(loaded$settings)$gc_clamp), c(0, 3)))
  scope$load_openprimer_settings <- function(...) loaded

  clamp <- 0
  tm_passed <- FALSE
  full_passed <- FALSE
  calls <- character()
  scope$evaluate_openprimer_form <- function(forward, reverse, template_sequence,
                                            settings, active_constraints, identifier) {
    full <- "self_dimerization" %in% active_constraints
    calls <<- c(calls, if (full) "full" else "annealing")
    stopifnot(identical(forward, if (full) "GGACGT" else "ACGT"),
              identical(reverse, if (full) "CCAACG" else "AACG"))
    metrics <- if (full) {
      data.frame(EVAL_self_dimerization = full_passed,
                 Self_Dimer_DeltaG = if (full_passed) 0 else -20)
    } else {
      data.frame(EVAL_gc_clamp = clamp <= 3, gc_clamp_fw = clamp, gc_clamp_rev = 2,
                 EVAL_melting_temp_range = tm_passed,
                 Tm_C_fw = if (tm_passed) 60 else 54, Tm_C_rev = 60,
                 melting_temp_diff = if (tm_passed) 0 else 6)
    }
    list(metrics = metrics, passed = if (full) full_passed else tm_passed && clamp <= 3,
         penalty = if (full && !full_passed) 20 else 0, active_constraints = active_constraints)
  }
  direct <- evaluate_openprimer_pair
  cached <- cached_openprimer_pair
  environment(direct) <- environment(cached) <- scope
  args <- list(annealing_forward = "ACGT", annealing_reverse = "AACG",
               full_forward = "GGACGT", full_reverse = "CCAACG",
               template_sequence = "ACGTAACG", reaction = "test")
  result <- do.call(direct, args)
  stopifnot(identical(calls, c("annealing", "full")),
            !result$passed, result$max_dimer_risk == 20,
            grepl("EVAL_self_dimerization", result$rejection_reason),
            grepl("low_gc_clamp", result$warnings), result$low_gc_clamp_count == 1)

  calls <- character()
  input <- list(primer_qc_cache = new.env(), parameters = list(
    primer_qc = primer_qc_defaults(), primer3_buffer = primer3_buffer_parameters()))
  cached_result <- do.call(cached, c(list(input = input), args))
  again <- do.call(cached, c(list(input = input), args))
  stopifnot(identical(calls, c("annealing", "full")),
            identical(cached_result, again), cached_result$max_dimer_risk == 20,
            grepl("EVAL_self_dimerization", cached_result$rejection_reason))

  # A zero clamp is advisory; 1 and 3 pass; 4 is a distinct high-clamp failure.
  tm_passed <- TRUE
  full_passed <- TRUE
  for (value in c(0, 1, 3, 4)) {
    clamp <- value
    result <- do.call(direct, args)
    stopifnot(result$passed == (value <= 3),
              grepl("low_gc_clamp", result$warnings) == (value == 0),
              grepl("high_gc_clamp", result$rejection_reason) == (value > 3))
  }
  specificity <- list(passed = TRUE, n_expected_products = 1L,
    n_high_risk_offtarget_products = 0, n_perfect_3p_offtarget_sites = 0)
  advisory <- list(passed = TRUE, abs_tm_diff = 0, warnings = "low_gc_clamp[gc_clamp_fw=0 (<1)]")
  policy <- evaluate_filtering_policy(specificity, advisory, 3L)
  stopifnot(policy$strict_passed, policy$blocking_passed,
            grepl("low_gc_clamp", policy$warnings))

  candidates <- data.frame(pair_id = c("zero", "short"), structure_passed = TRUE,
    blocking_passed = TRUE, strict_qc_passed = TRUE,
    n_high_risk_offtarget_products = 0, openprimer_penalty = 0,
    low_gc_clamp_count = c(1, 0), primer3_index = 1:2)
  stopifnot(select_best_primer_pair(candidates)$pair$pair_id == "short")
  candidates$openprimer_penalty <- c(0, 1)
  stopifnot(select_best_primer_pair(candidates)$pair$pair_id == "zero")
})
message("GC-clamp policy and full-sequence fallback QC regressions passed")
