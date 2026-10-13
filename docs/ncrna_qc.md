# ncRNA/CDS quality-control audit

The strategy switch changes deletion geometry and search order. The following
checks use the same production functions and filtering configuration in both
paths.

| Check | Shared implementation | Policy |
|---|---|---|
| N20 off-target thresholds, self-complementarity, strand, genomic sequence | `filter_grnas`, `prepare_grna_pool` | Same filters; midpoint-distance filter is CDS geometry only |
| Primer length, Tm, GC, clamp, homopolymers, self/pair complementarity | `write_primer3_settings`, `callPrimer3` | Same generation settings at all filtering levels |
| Restriction sites in both arms, each N20/scaffold, assembly junctions, cassette and edited pTarget | `check_candidate_restriction_sites`, `model_edited_ptargets` | Mandatory at every level, including CTGCAG |
| Newly introduced genomic restriction sites | `check_genome_edit_restriction_sites` | Mandatory before candidate acceptance and again during output |
| Binding-site mismatches, critical 3′ ends, intended and off-target PCR products | `evaluate_candidate_reaction`, `cached_pair_specificity`, `evaluate_pair_specificity` | Same `evaluate_filtering_policy` gates |
| Annealing/full-primer thermodynamics, self/cross dimers, secondary structure, repeats, efficiency | `cached_openprimer_pair`, `evaluate_openprimer_pair` | Same openPrimeR configuration, caching and fallback policy |
| Screening primers in original and edited references | `write_design_outputs` | Same structural and QC gates, ranking and warnings |
| Final assembly, edited sequences and PCR products | `write_design_outputs`, `model_design_pcr_products`, `simulate_full_primer_pcr` | Same checks; a failed ncRNA output tries the next arm candidate before changing N20 |

CDS retains its inward search, internal arm boundaries, partial deletion and
phase-preserving bridge prefix. ncRNA requires a complete target deletion within
neighbor-buffer bounds and uses the full bridge. The existing independent
`--cds-fs` / `--ncrna-fs` flags still require a deletion divisible by three when
explicitly enabled. Arm length and actual N20-to-arm distances remain mandatory.

For ncRNA, each inner primer endpoint is generated once for each side, retaining
earlier candidates. Both windows expand independently by one nucleotide. Pairs
are tested by total added deletion length, so a 3-nt one-sided extension is tried
before a 2+2-nt extension. Within that length, the existing strict-QC preference
and risk ranking apply. The mandatory gates have the same meaning as CDS;
ncRNA does not postpone a passing smaller deletion to seek a larger strict-QC
candidate.

The Primer3 endpoint options can return primers with constraint violations, as
documented in the [Primer3 manual](https://primer3.org/manual.html#SEQUENCE_FORCE_LEFT_START).
The wrapper excludes forced candidates carrying `PRIMER_*_PROBLEMS`; ordinary
CDS generation is unchanged. This exclusion is separate from optional QC
fallbacks and cannot be relaxed by the filtering level.

`test/test_ncrna.R` exercises real geometry, shared restriction checks and shared
filtering policy, with synthetic Primer3 candidates and QC measurements. It
checks both strands, asymmetric and minimal total expansion, N20 offsets,
neighbor conflicts, contig edges, attempt limits, reuse of generated candidates,
output rejection, strategy resolution and failure reporting.
`test/test_restriction_sites.R` covers restriction junctions and complete
constructs; `test/test_qc_policy.R` covers detailed openPrimeR policy. Full
thermodynamic/tool integration additionally requires the project's Linux Conda
environments and external tools.
