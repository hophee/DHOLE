#!/usr/bin/env bash
set -uo pipefail

readonly TEST_SCRIPT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)/$(basename -- "${BASH_SOURCE[0]}")"
readonly CONDA_ENV_NAME="oligo_design"

if [[ "${CONDA_DEFAULT_ENV:-}" != "$CONDA_ENV_NAME" ]]; then
  command -v conda >/dev/null 2>&1 || {
    printf 'TEST FAILED: conda is not available in PATH\n' >&2
    exit 1
  }
  exec conda run --no-capture-output --name "$CONDA_ENV_NAME" \
    bash "$TEST_SCRIPT" "$@"
fi

readonly TEST_DIR="$(dirname -- "$TEST_SCRIPT")"
readonly PROJECT_DIR="$(cd -- "$TEST_DIR/.." && pwd)"
readonly OUTPUT_DIR="$TEST_DIR/test_output"
readonly WET_LAB_DIR="$OUTPUT_DIR/WetLab"
readonly TECH_REPORT_DIR="$OUTPUT_DIR/TechReport"
readonly SUMMARY="$TECH_REPORT_DIR/design_summary.tsv"
readonly PLASMID="$TEST_DIR/test_target_plasmid.fasta"
readonly CAS_PLASMID="$TEST_DIR/test_cas_plasmid.fasta"
readonly TEST_GENES="recA,pta,hupB"

fail() {
  printf 'TEST FAILED: %s\n' "$1" >&2
  exit 1
}

[[ -s "$TEST_DIR/MG1655.fna" ]] || fail "MG1655.fna is missing or empty"
[[ -s "$TEST_DIR/MG1655.gff" ]] || fail "MG1655.gff is missing or empty"

cat > "$PLASMID" <<'EOF'
>synthetic_test_target_plasmid
GGGACTAGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTCCCCGGCACCGAGTCGGTGCTTTTTTTAGCGTCAACTCTGCAGAAAA
EOF

cat > "$CAS_PLASMID" <<'EOF'
>synthetic_test_cas_plasmid
TTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGT
EOF

rm -rf "$OUTPUT_DIR"
cd "$PROJECT_DIR" || fail "cannot enter project directory"

Rscript test/test_unit.R || fail "unit tests failed"
Rscript test/test_regressions.R || fail "regression tests failed"
Rscript test/test_screening_fixture.R || fail "screening fixture failed"

if ! Rscript oligo_designer.R \
  --genome "$TEST_DIR/MG1655.fna" \
  --genome-annotation "$TEST_DIR/MG1655.gff" \
  --annotation-format gff \
  --target-plasmid "$PLASMID" \
  --cas-plasmid "$CAS_PLASMID" \
  --output-dir "$OUTPUT_DIR" \
  --left-arm-max 350 \
  --right-arm-max 450 \
  --cds "$TEST_GENES"; then
  fail "oligo_designer.R returned a non-zero exit code"
fi

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
grep -q $'^primer3_buffer_divalent_salt_mm\t1.5$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks Primer3 buffer data"
grep -q $'^primer_qc_critical_3p_bases\t5$' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks primer QC defaults"
grep -q $'^openprimer_active_constraints\t' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks active openPrimeR constraints"
grep -q $'^cas_plasmid_file\t' "$TECH_REPORT_DIR/run_parameters.tsv" ||
  fail "run_parameters.tsv lacks the pCas input path"

for gene in recA pta hupB; do
  target_dir="$TECH_REPORT_DIR/${gene,,}_results"
  wet_target_dir="$WET_LAB_DIR/${gene,,}_results"
  for result in n20_table.tsv design.log primer_binding_sites.tsv primer_amplicons.tsv primer_openprimer_qc.tsv primer_pair_ranking.tsv; do
    [[ -s "$target_dir/$result" ]] || fail "missing or empty result for $gene: $result"
  done

  status="$(awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene { print $4 }' "$SUMMARY")"
  if [[ "$status" == "error" ]]; then
    stage="$(awk -F '\t' -v gene="$gene" 'NR > 1 && $1 == gene { print $5 }' "$SUMMARY")"
    [[ "$stage" == "primer_qc" || "$stage" == "homology_arms" ]] ||
      fail "$gene failed outside the documented strict primer_qc gate: $stage"
    if [[ "$stage" == "primer_qc" ]]; then
      grep -q $'\tprimer_qc\tTRY\t' "$target_dir/design.log" ||
        fail "design.log lacks primer_qc TRY for $gene"
      grep -q $'\tprimer_qc\tREJECTED\t' "$target_dir/design.log" ||
        fail "design.log lacks an explained candidate rejection for $gene"
    fi
    [[ -s "$target_dir/error.txt" ]] || fail "$gene lacks error.txt"
    [[ ! -d "$wet_target_dir" ]] ||
      fail "$gene has WetLab output after primer_qc rejection"
    try_count="$(grep -c $'\thomology_arms\tTRY\t' "$target_dir/design.log")"
    [[ "$try_count" -gt 1 ]] ||
      fail "$gene did not continue to the next N20 set after rejection"
    printf 'DOCUMENTED DESIGN CHANGE: %s has no selected pair (stage=%s).\n' "$gene" "$stage"
    continue
  fi
  [[ "$status" == "ok" ]] || fail "$gene has unexpected summary status: $status"
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
  fasta_records="$(grep -c '^>' "$wet_target_dir/final_sequences.fasta")"
  txt_records="$(awk 'NR > 1 { count++ } END { print count + 0 }' "$wet_target_dir/final_sequences.txt")"
  [[ "$fasta_records" -eq "$txt_records" ]] ||
    fail "WetLab FASTA and TXT sequence sets differ for $gene"
  grep -q 'Без успешного нокаута' "$wet_target_dir/wet_lab_report.txt" ||
    fail "WetLab report lacks the unsuccessful-knockout PCR size for $gene"
  grep -q 'С успешным нокаутом' "$wet_target_dir/wet_lab_report.txt" ||
    fail "WetLab report lacks the successful-knockout PCR size for $gene"
  grep -q 'DECIPHER::AmplifyDNA' "$wet_target_dir/wet_lab_report.txt" ||
    fail "WetLab report lacks modelled PCR products for $gene"
  ptarget_records="$(grep -c '^>' "$wet_target_dir/edited_pTargets.fasta")"
  [[ "$ptarget_records" -eq 1 ]] ||
    fail "WetLab does not contain one edited pTarget per selected N20 for $gene"
  grep -q '^ACTAGT' <(awk '!/^>/ { print; exit }' "$wet_target_dir/edited_pTargets.fasta") ||
    fail "edited pTarget does not start with one intact site1 for $gene"
  grep -q $'\tprimer_qc\tOK\t' "$target_dir/design.log" ||
    fail "design.log lacks primer_qc OK for $gene"
  Rscript - "$target_dir" <<'EOF' || fail "selected primer QC trace is invalid for $gene"
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
stopifnot(!any(selected_amplicons$off_target & !selected_amplicons$invalid_size))
EOF
done

grep -Rq $'\tprimer_qc\tTRY\t' "$TECH_REPORT_DIR"/*_results/design.log ||
  fail "integration run did not exercise physical primer QC"
grep -Rq $'\tprimer_qc\tREJECTED\t' "$TECH_REPORT_DIR"/*_results/design.log ||
  fail "integration run did not record a candidate rejection"

printf 'TEST PASSED: MG1655 recA, pta and hupB produced valid selected designs or explicit strict-QC rejection traces.\n'
