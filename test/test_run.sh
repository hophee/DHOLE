#!/usr/bin/env bash
set -euo pipefail

readonly TEST_SCRIPT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
readonly TEST_DIR="$(dirname -- "$TEST_SCRIPT")"
readonly PROJECT_DIR="$(cd -- "$TEST_DIR/.." && pwd)"
readonly CONDA_ENV_NAME="oligo_design"
readonly R_RUNNER="$PROJECT_DIR/tools/run-r"
readonly OUTPUT_DIR="$TEST_DIR/test_output"
readonly WET_LAB_DIR="$OUTPUT_DIR/WetLab"
readonly TECH_REPORT_DIR="$OUTPUT_DIR/TechReport"
readonly SUMMARY="$TECH_REPORT_DIR/design_summary.tsv"
readonly FIXTURE_DIR="$OUTPUT_DIR/fixtures"
readonly PLASMID="$FIXTURE_DIR/test_target_plasmid.fasta"
readonly CAS_PLASMID="$FIXTURE_DIR/test_cas_plasmid.fasta"
readonly TEST_GENES="recA,pta,hupB"
readonly LOG_DIR="$TEST_DIR/test_logs"
readonly STAGE_TOTAL=6

verbose=0
for arg in "$@"; do
  case "$arg" in
    --verbose) verbose=1 ;;
    *)
      printf 'TEST FAILED: unknown argument: %s\n' "$arg" >&2
      exit 2
      ;;
  esac
done

mkdir -p -- "$LOG_DIR" || {
  printf 'TEST FAILED: cannot create log directory: %s\n' "$LOG_DIR" >&2
  exit 1
}

if [[ -n "${DHOLE_TEST_LOG:-}" ]]; then
  readonly TEST_LOG="$DHOLE_TEST_LOG"
  new_log=0
else
  readonly TEST_LOG="$LOG_DIR/test-$(date -u +%Y%m%dT%H%M%SZ)-$$.log"
  : > "$TEST_LOG" || {
    printf 'TEST FAILED: cannot create log: %s\n' "$TEST_LOG" >&2
    exit 1
  }
  export DHOLE_TEST_LOG="$TEST_LOG"
  new_log=1
fi

current_stage="bootstrap"
failure_reported=0
had_warnings=0

console() {
  printf '%s\n' "$*"
}

status_line() {
  console "$*"
  if ! printf '%s\n' "$*" >> "$TEST_LOG"; then
    printf 'TEST FAILED: cannot write log: %s\n' "$TEST_LOG" >&2
    failure_reported=1
    exit 74
  fi
}

log_command() {
  printf 'COMMAND' >> "$TEST_LOG"
  printf '\t%s' "$@" >> "$TEST_LOG"
  printf '\n' >> "$TEST_LOG"
}

report_failure() {
  local status="$1"
  local detail="$2"
  failure_reported=1
  status_line "OVERALL RESULT: FAIL ($detail; exit=$status)"
  status_line "Log: $TEST_LOG"
}

on_exit() {
  local status=$?
  if [[ "$status" -ne 0 && "$failure_reported" -eq 0 ]]; then
    report_failure "$status" "$current_stage"
  fi
}

trap on_exit EXIT
trap 'exit 130' INT
trap 'exit 143' TERM

if [[ "$new_log" -eq 1 ]]; then
  {
    printf 'DHOLE test run\n'
    printf 'started_at=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'project=%s\n' "$PROJECT_DIR"
    printf 'git_commit=%s\n' "$(git -C "$PROJECT_DIR" rev-parse HEAD 2>/dev/null || printf unknown)"
    log_command "$TEST_SCRIPT" "$@"
  } >> "$TEST_LOG"
  console 'DHOLE test run'
  console "Log: $TEST_LOG"
fi

if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV_NAME" ]]; then
  command -v conda >/dev/null 2>&1 || {
    report_failure 1 "conda is not available in PATH"
    exit 1
  }
  current_stage="enter Conda environment $CONDA_ENV_NAME"
  status_line "[bootstrap] Entering Conda environment: $CONDA_ENV_NAME"
  if conda run --no-capture-output --name "$CONDA_ENV_NAME" \
    bash "$TEST_SCRIPT" "$@" 2>> "$TEST_LOG"; then
    exit 0
  else
    status=$?
    if grep -q '^OVERALL RESULT: FAIL' "$TEST_LOG"; then
      failure_reported=1
    else
      report_failure "$status" "$current_stage"
      tail -n 20 "$TEST_LOG" >&2 || true
    fi
    exit "$status"
  fi
fi

run_stage() {
  local number="$1"
  local label="$2"
  local warning_exit="$3"
  shift 3
  local started=$SECONDS
  local status
  local -a pipe_status

  current_stage="$label"
  status_line "[$number/$STAGE_TOTAL] $label ..."
  {
    printf '\n=== [%s/%s] %s: START ===\n' "$number" "$STAGE_TOTAL" "$label"
    log_command "$@"
  } >> "$TEST_LOG"

  if [[ "$verbose" -eq 1 ]]; then
    set +e
    ( "$@" ) 2>&1 | tee -a "$TEST_LOG"
    pipe_status=("${PIPESTATUS[@]}")
    set -e
    status="${pipe_status[0]}"
    if [[ "${pipe_status[1]}" -ne 0 ]]; then
      status=74
    fi
  elif ( "$@" ) >> "$TEST_LOG" 2>&1; then
    status=0
  else
    status=$?
  fi

  local elapsed=$((SECONDS - started))
  case "$status" in
    0)
      status_line "[$number/$STAGE_TOTAL] $label: OK (${elapsed}s)"
      ;;
    "$warning_exit")
      had_warnings=1
      status_line "[$number/$STAGE_TOTAL] $label: WARN (${elapsed}s; details in log)"
      ;;
    *)
      status_line "[$number/$STAGE_TOTAL] $label: FAIL (${elapsed}s; exit=$status)"
      status_line "Last 20 log lines:"
      tail -n 20 "$TEST_LOG" >&2 || true
      report_failure "$status" "$label"
      exit "$status"
      ;;
  esac
  printf '=== [%s/%s] %s: END status=%s elapsed=%ss ===\n' \
    "$number" "$STAGE_TOTAL" "$label" "$status" "$elapsed" >> "$TEST_LOG"
}

fail() {
  printf 'TEST FAILED: %s\n' "$1" >&2
  exit 1
}

prepare_workspace() {
  rm -rf -- "$OUTPUT_DIR" || fail "cannot remove $OUTPUT_DIR"
  mkdir -p -- "$FIXTURE_DIR" || fail "cannot create $FIXTURE_DIR"
  [[ -s "$TEST_DIR/MG1655.fna" ]] || fail "MG1655.fna is missing or empty"
  [[ -s "$TEST_DIR/MG1655.gff" ]] || fail "MG1655.gff is missing or empty"

  printf '>synthetic_test_target_plasmid\n%s\n' \
    "GCAGGGGACTAGTACGTACGTACGTACGTACGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGAATTCTCTAGAGTCGACCTGCAG$(printf 'A%.0s' {1..200})" \
    > "$PLASMID"

  cat > "$CAS_PLASMID" <<'EOF'
>synthetic_test_cas_plasmid
TTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGT
EOF
}

run_designer_with_immediate_warnings() {
  log_command bash "$R_RUNNER" - "$PROJECT_DIR/oligo_designer.R" \
    --genome "$TEST_DIR/MG1655.fna" \
    --genome-annotation "$TEST_DIR/MG1655.gff" \
    --annotation-format gff \
    --target-plasmid "$PLASMID" \
    --cas-plasmid "$CAS_PLASMID" \
    --output-dir "$OUTPUT_DIR" \
    --left-arm-max 350 \
    --right-arm-max 450 \
    --cds "$TEST_GENES"

  bash "$R_RUNNER" - "$PROJECT_DIR/oligo_designer.R" \
    --genome "$TEST_DIR/MG1655.fna" \
    --genome-annotation "$TEST_DIR/MG1655.gff" \
    --annotation-format gff \
    --target-plasmid "$PLASMID" \
    --cas-plasmid "$CAS_PLASMID" \
    --output-dir "$OUTPUT_DIR" \
    --left-arm-max 350 \
    --right-arm-max 450 \
    --cds "$TEST_GENES" <<'RSCRIPT'
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
script <- args[[1L]]
source(script, local = globalenv())
main(args[-1L])
RSCRIPT
}

run_integration() {
  cd "$PROJECT_DIR" || fail "cannot enter project directory"
  run_designer_with_immediate_warnings ||
    fail "oligo_designer.R returned a non-zero exit code"

  [[ -s "$SUMMARY" ]] || fail "design_summary.tsv was not created"
  [[ ! -e "$OUTPUT_DIR/design_summary.tsv" ]] ||
    fail "design_summary.tsv was written outside TechReport"
  [[ -s "$TECH_REPORT_DIR/run_parameters.tsv" ]] ||
    fail "run_parameters.tsv was not created"
  [[ -s "$TECH_REPORT_DIR/chopchop_config.json" ]] ||
    fail "CHOPCHOP configuration snapshot was not created"
  [[ -s "$TECH_REPORT_DIR/genome_indexes/MG1655.2bit" ]] ||
    fail "2bit genome index was not written to TechReport"
  [[ -s "$TECH_REPORT_DIR/genome_indexes/MG1655.1.ebwt" ]] ||
    fail "Bowtie genome index was not written to TechReport"
  grep -q $'^genome_file\t' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the genome input path"
  grep -q $'^n20_arm_min_distance_nt\t40$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the N20-to-arm distance"
  grep -q $'^ptarget_site1\tACTAGT$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the default site1"
  grep -q $'^ptarget_site2\tCTGCAG$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the default site2"
  grep -q $'^ptarget_cassette_arc\tshortest$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the default pTarget cassette arc"
  grep -q $'^ptarget_sgrna_scaffold\tGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the validated sgRNA scaffold"
  grep -q $'^ptarget_sgrna_annealing_temp_c\t60$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the default sgRNA annealing temperature"
  grep -q $'^primer3_buffer_divalent_salt_mm\t1.5$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks Primer3 buffer data"
  grep -q $'^primer_qc_critical_3p_bases\t5$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks primer QC defaults"
  grep -q $'^filtering_level\t2$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the default filtering level"
  grep -q $'^openprimer_active_constraints\t' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks active openPrimeR constraints"
  grep -q $'^cas_plasmid_file\t' "$TECH_REPORT_DIR/run_parameters.tsv" ||
    fail "run_parameters.tsv lacks the pCas input path"

  local success_count=0
  local rejected_count=0
  local fallback_count=0
  local gene target_dir wet_target_dir result status stage try_count
  local fasta_records txt_records ptarget_records
  for gene in recA pta hupB; do
    target_dir="$TECH_REPORT_DIR/${gene,,}_results"
    wet_target_dir="$WET_LAB_DIR/${gene,,}_results"
    for result in n20_table.tsv design.log primer_binding_sites.tsv primer_amplicons.tsv primer_openprimer_qc.tsv primer_pair_ranking.tsv; do
      [[ -s "$target_dir/$result" ]] || fail "missing or empty result for $gene: $result"
    done

    status="$(awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene { print $4 }' "$SUMMARY")"
    if [[ "$status" == "error" ]]; then
      rejected_count=$((rejected_count + 1))
      stage="$(awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene { print $5 }' "$SUMMARY")"
      [[ "$stage" == "primer_qc" || "$stage" == "homology_arms" ]] ||
        fail "$gene failed outside the documented primer/core-QC gate: $stage"
      if [[ "$stage" == "primer_qc" ]]; then
        grep -q $'\tprimer_qc\tTRY\t' "$target_dir/design.log" ||
          fail "design.log lacks primer_qc TRY for $gene"
        grep -q $'\tprimer_qc\tREJECTED\t' "$target_dir/design.log" ||
          fail "design.log lacks an explained candidate rejection for $gene"
      fi
      [[ -s "$target_dir/error.txt" ]] || fail "$gene lacks error.txt"
      [[ ! -d "$wet_target_dir" ]] ||
        fail "$gene has WetLab output after primer_qc rejection"
      try_count="$(grep -c $'\thomology_arms\tTRY\t' "$target_dir/design.log" || true)"
      [[ "$try_count" -gt 1 ]] ||
        fail "$gene did not continue to the next N20 set after rejection"
      printf 'DOCUMENTED DESIGN CHANGE: %s has no selected pair (stage=%s).\n' "$gene" "$stage"
      continue
    fi
    [[ "$status" == "ok" ]] || fail "$gene has unexpected summary status: $status"
    success_count=$((success_count + 1))
    awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene && $7 != "" && $7 != "NA" { found = 1 } END { exit !found }' "$SUMMARY" ||
      fail "$gene has no WetLab path in design_summary.tsv"
    for result in all_primers.fasta edited_genome.fasta edited_pTargets.fasta pcr_products.fasta pcr_products.tsv report.tsv; do
      [[ -s "$target_dir/$result" ]] || fail "missing successful result for $gene: $result"
    done
    for result in final_sequences.fasta final_sequences.txt wet_lab_report.txt edited_genome.fasta edited_pTargets.fasta pcr_products.fasta pcr_products.tsv; do
      [[ -s "$wet_target_dir/$result" ]] ||
        fail "missing or empty WetLab result for $gene: $result"
    done
    [[ ! -e "$wet_target_dir/design.log" ]] ||
      fail "WetLab unexpectedly contains a technical design.log for $gene"
    [[ ! -e "$wet_target_dir/n20_table.tsv" ]] ||
      fail "WetLab unexpectedly contains a CHOPCHOP table for $gene"
    fasta_records="$(grep -c '^>' "$wet_target_dir/final_sequences.fasta" || true)"
    txt_records="$(awk 'NR > 1 { count++ } END { print count + 0 }' "$wet_target_dir/final_sequences.txt")"
    [[ "$fasta_records" -eq "$txt_records" ]] ||
      fail "WetLab FASTA and TXT sequence sets differ for $gene"
    grep -q 'Без успешного нокаута' "$wet_target_dir/wet_lab_report.txt" ||
      fail "WetLab report lacks the unsuccessful-knockout PCR size for $gene"
    grep -q 'С успешным нокаутом' "$wet_target_dir/wet_lab_report.txt" ||
      fail "WetLab report lacks the successful-knockout PCR size for $gene"
    grep -q 'DECIPHER::AmplifyDNA' "$wet_target_dir/wet_lab_report.txt" ||
      fail "WetLab report lacks modelled PCR products for $gene"
    ptarget_records="$(grep -c '^>' "$wet_target_dir/edited_pTargets.fasta" || true)"
    [[ "$ptarget_records" -eq 1 ]] ||
      fail "WetLab does not contain one edited pTarget per selected N20 for $gene"
    grep -q '^ACTAGT' <(awk '!/^>/ { print; exit }' "$wet_target_dir/edited_pTargets.fasta") ||
      fail "edited pTarget does not start with one intact site1 for $gene"
    grep -q $'\tprimer_qc\tOK\t' "$target_dir/design.log" ||
      fail "design.log lacks primer_qc OK for $gene"
    if grep -q $'^primer_qc_fallback_used\tTRUE$' "$target_dir/report.tsv"; then
      fallback_count=$((fallback_count + 1))
      grep -q $'\tprimer_qc\tWARNING\t' "$target_dir/design.log" ||
        fail "design.log lacks a fallback warning for $gene"
      grep -q 'Предупреждения primer QC' "$wet_target_dir/wet_lab_report.txt" ||
        fail "WetLab report lacks fallback warnings for $gene"
    fi
    bash "$R_RUNNER" - "$target_dir" <<'EOF' || fail "selected primer QC trace is invalid for $gene"
args <- commandArgs(trailingOnly = TRUE)
target_dir <- args[[1]]
ranking <- read.delim(file.path(target_dir, "primer_pair_ranking.tsv"), check.names = FALSE)
amplicons <- read.delim(file.path(target_dir, "primer_amplicons.tsv"), check.names = FALSE)
selected <- ranking[ranking$selected %in% TRUE, , drop = FALSE]
stopifnot(nrow(selected) >= 2L)
physical_ids <- unique(na.omit(c(
  selected$left_pair_id,
  selected$right_pair_id,
  selected$pair_id[selected$reaction == "scrF_scrR"]
)))
selected_amplicons <- amplicons[amplicons$pair_id %in% physical_ids, , drop = FALSE]
stopifnot(all(c("LF_LR", "RF_RR", "scrF_scrR") %in% selected_amplicons$reaction))
stopifnot(all(vapply(
  c("LF_LR", "RF_RR", "scrF_scrR"),
  function(reaction) any(selected_amplicons$reaction == reaction & selected_amplicons$intended),
  logical(1)
)))
stopifnot(all(selected$structure_passed & selected$blocking_passed))
if (any(selected$fallback_selected)) {
  stopifnot(all(nzchar(selected$risk_warnings[selected$fallback_selected])))
} else {
  stopifnot(all(selected$strict_qc_passed))
}
EOF
  done

  [[ "$success_count" -gt 0 ]] ||
    fail "default filtering did not produce any design for the test genes"
  grep -Rq $'\tprimer_qc\tTRY\t' "$TECH_REPORT_DIR"/*_results/design.log ||
    fail "integration run did not exercise physical primer QC"
  grep -Rq $'\tprimer_qc\tREJECTED\t' "$TECH_REPORT_DIR"/*_results/design.log ||
    fail "integration run did not record a candidate rejection"

  printf 'TARGET SUMMARY: ok=%d expected_rejection=%d qc_fallback=%d total=3\n' \
    "$success_count" "$rejected_count" "$fallback_count"
  if [[ "$rejected_count" -gt 0 || "$fallback_count" -gt 0 ]]; then
    return 3
  fi
}

run_stage 1 'Prepare clean workspace' 0 prepare_workspace
run_stage 2 'R environment isolation' 0 bash "$TEST_DIR/test_r_environment.sh"
run_stage 3 'Unit tests' 0 bash "$R_RUNNER" "$TEST_DIR/test_unit.R"
run_stage 4 'Targeted regressions' 0 bash "$R_RUNNER" "$TEST_DIR/test_regressions.R"
run_stage 5 'Screening integration fixture' 0 bash "$R_RUNNER" "$TEST_DIR/test_screening_fixture.R"
run_stage 6 'End-to-end design and artifact validation' 3 run_integration

target_summary="$(grep '^TARGET SUMMARY:' "$TEST_LOG" | tail -n 1 || true)"
if [[ -n "$target_summary" ]]; then
  console "  $target_summary"
fi

if [[ "$verbose" -eq 0 ]]; then
  qc_violation_summary="$(
    awk '
      /^=== \[6\/6\].*: START ===$/ { in_integration = 1; next }
      /^=== \[6\/6\].*: END / { in_integration = 0 }
      in_integration && /Выбраны праймеры с QC-рисками/ && !seen[$0]++
    ' "$TEST_LOG" || true
  )"
  if [[ -n "$qc_violation_summary" ]]; then
    console 'QC threshold violations:'
    console "$qc_violation_summary"
  fi
fi

if [[ "$had_warnings" -eq 1 ]]; then
  status_line 'OVERALL RESULT: PASS WITH WARNINGS'
else
  status_line 'OVERALL RESULT: PASS'
fi
status_line "Log: $TEST_LOG"
