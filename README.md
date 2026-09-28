# DHOLE

DHOLE designs CRISPR-Cas9 N20 oligos, homology-arm PCR primers, screening
primers, and an edited-genome model for bacterial CDS and ncRNA targets.

Primer selection is candidate-based: Primer3 proposes up to ten pairs per
homology arm and five screening pairs, then 2PAC applies structural rules,
exhaustive Biostrings specificity, openPrimeR QC, and a stable lexicographic
ranking. Three filtering levels control which primer-QC risks block selection.
A rejected first Primer3 row no longer rejects the target. If no ideal row is
available, the best candidate allowed by the selected level is returned with
explicit warnings.

## Installation

The environment is split to preserve the legacy Python 2 CHOPCHOP stack:

- `env.yml`: R, pak, compilers, Bowtie, OligoArrayAux, MAFFT, and the remaining
  command-line tools;
- `env_chopchop.yml`: CHOPCHOP and its Python 2 dependencies;
- `env_viennarna.yml`: ViennaRNA, exposed to the main environment by a small
  launcher created by `install.sh`.

Run `./install.sh` for a complete installation. The installer uses
`pak::pkg_install()` for CRAN/Bioconductor dependencies and explicitly targets
R's `.Library` inside `oligo_design`; it does not install them into a user
library. Missing openPrimeR constraints or executables are reported as QC
risks. Level 3 uses the Primer3 pair ΔTm as a fallback for its temperature
gate; it blocks selection only when that value is unavailable or exceeds the
hard limit.

Java remains required by the `rmelting` JAR behind `tools/melting-batch` and is
installed through Conda (`openjdk`). `pak` may still report the OS package
`java-11-openjdk-devel` as missing: `PKG_SYSREQS=false` disables installation
of OS packages, but not their reporting (see the
[pak configuration reference](https://pak.r-lib.org/reference/pak-config.html)).
The final installer check loads `rJava` and `rmelting`, initializes the JVM,
and checks for `java` and `javac`; a failed check makes installation fail.

Use `bash tools/run-r` in place of `Rscript` for this project. It selects
`oligo_design`, ignores R startup files, and excludes inherited `R_LIBS*`
paths. The installer, test runner, and MELTING package lookup use it too.
This prevents packages built for another R version in a shared user library
from causing errors such as `rlang.so: undefined symbol: R_MakeMissingBinding`.
`Rscript --vanilla` alone does not clear inherited library paths.

## Usage

```bash
bash tools/run-r oligo_designer.R \
  --genome genome.fasta \
  --genome-annotation genome.gff \
  --annotation-format gff \
  --target-plasmid pTarget.fasta \
  --filtering-level 2 \
  --site1 ACTAGT \
  --site2 CTGCAG \
  --cas-plasmid pCas.fasta \
  --output-dir results \
  --cds gene1,gene2 \
  --ncrna rna1,rna2
```

Targets are matched first by `locus_tag`, then by `gene`. `--cds` and
`--ncrna` accept comma- or space-separated values. Bakta TSV and GFF/GFF3
annotations are supported. The current implementation supports one complete
linear genome contig. Every FASTA record in pTarget and pCas is treated as a
separate circular specificity reference, while edited-pTarget modelling
requires exactly one pTarget record.

### Design parameters

| Argument | Default | Meaning |
|---|---:|---|
| `--filtering-level` | `2` | `1` lite/report, `2` default/warn, `3` hard/core |
| `--n20-mn` | `1` | Required N20 count |
| `--n20-strands` | `random` | `plus`, `minus`, `both`, or unconstrained `random` |
| `--n20-offtarget` | `0` | Maximum CHOPCHOP `MM0,MM1,...` values |
| `--site1` | `ACTAGT` (SpeI) | First restriction-site sequence in insert orientation |
| `--site2` | `CTGCAG` (PstI) | Second restriction-site sequence in insert orientation |
| `--ptarget-cassette-arc` | `shortest` | Arc between the sites used as the cassette: `shortest`, `forward`, or `reverse` |
| `--sgrna-annealing-temp-c` | `60` | Annealing temperature for the sgRNA-cassette PCR |
| `--cds-fs`, `--ncrna-fs` | off | Require deleted length divisible by three |
| `--left-arm-min/opt/max` | `300/350/400` | Left-arm structural limits |
| `--right-arm-min/opt/max` | `400/450/500` | Right-arm structural limits |
| `--n20-arm-min-distance` | `40` | Minimum N20-to-arm distance in nt |
| `--primer-max-mismatches` | `2` | Maximum mismatches per binding site |
| `--primer-critical-3p-bases` | `5` | Critical primer 3′ region |
| `--primer-max-3p-mismatches` | `0` | Allowed mismatches in that region |
| `--primer-min-product-size` | `50` | Minimum counted amplicon size |
| `--primer-max-product-size` | `2000` | Maximum counted amplicon size |
| `--primer-max-offtarget-products` | `0` | Preferred maximum non-intended amplicons |

The legacy Primer3 generation thresholds remain unchanged: primer length
`18/21/27` nt (min/opt/max), homopolymer maximum `5`, and pair Tm difference
maximum `8 °C`. Candidates outside the stricter openPrimeR profile remain in
the trace and may be used as a warned fallback according to the filtering
level. Primer3 buffer defaults
remain 50 mM monovalent salt, 1.5 mM Mg, 0.6 mM dNTP, and 50 nM DNA.
`run_parameters.tsv` records these values, the active constraints, their
effective limits, package versions, and tool paths.

## Selection policy

Arm geometry, feature bounds, N20 distance, frame, deletion rules, and Primer3
generation are identical at all filtering levels. The levels affect only PCR
specificity and openPrimeR QC:

1. `lite` — report only. Specificity and openPrimeR metrics rank candidates but
   do not block a structurally valid Primer3 pair.
2. `default` — warn and continue. Exactly one intended PCR product is required;
   off-target and openPrimeR failures are ranked risks and may be returned with
   warnings.
3. `hard` — core gates. Exactly one intended product, zero high-risk off-target
   products with perfect primer 3′ ends, and an evaluated pair ΔTm of at most
   `5 °C` are mandatory. Other openPrimeR and specificity criteria remain
   optimization criteria and may appear in a warned fallback.

At every level a candidate passing the former complete strict profile is
preferred. If none exists, the best candidate passing the level's core gates is
selected by intended-product availability, off-target risk, openPrimeR penalty,
dimer risk, Tm difference, Primer3 penalty, deleted nucleotides, and original
Primer3 order. Core gates cannot be compensated by a score.

## Output

```text
results/
├── WetLab/<target>_results/
│   ├── edited_genome.fasta
│   ├── edited_pTargets.fasta
│   ├── final_sequences.fasta
│   ├── final_sequences.txt
│   ├── pcr_products.fasta
│   ├── pcr_products.tsv
│   └── wet_lab_report.txt
└── TechReport/
    ├── design_summary.tsv
    ├── run_parameters.tsv
    └── <target>_results/
        ├── primer_binding_sites.tsv
        ├── primer_amplicons.tsv
        ├── primer_openprimer_qc.tsv
        ├── primer_pair_ranking.tsv
        ├── report.tsv
        └── design.log
```

The four primer QC tables preserve every evaluated binding site, amplicon,
openPrimeR metric/`EVAL_*` result, filtering level, strict/core gate, fallback
flag, warning, rank component, selection flag, and rejection reason.
`design.log` records `primer_qc TRY`, `REJECTED`, `OK`, and `WARNING`.
For every successful target, `wet_lab_report.txt` contains the complete final
oligo set, modelled-construction names and lengths, primer Tm values, expected
screening products for edited and unedited alleles, per-N20 distances to both
homology arms, screening off-target counts, and selected openPrimeR quality
metrics with readable labels. `edited_pTargets.fasta` contains one circular
pTarget model for every selected N20. Site matching checks both orientations
and the FASTA origin; each supplied site must occur exactly once physically.
For palindromic sites, the shorter circular arc is treated as the cassette by
default. Use `--ptarget-cassette-arc forward` or `reverse` when the cassette is
the longer arc or an explicit orientation is required. `forward` means the
input FASTA strand from `site1` to `site2`; `reverse` means the reverse-
complement strand from `site1` to `site2`.

`pcr_products.tsv` and the WetLab report contain one sgRNA-cassette PCR product
per N20, both homology-arm products, and screening products from the original
and edited genomes. Every row records its template location, product length,
full primer names, PCR simulation conditions, and sequence. Products are
generated by `DECIPHER::AmplifyDNA` with the complete service-tailed primers.
The restriction sites only delimit the cassette replaced in the edited-pTarget
model; there is no independent cassette-length threshold such as 59 nt. For
sgRNA-primer generation, the selected cassette must start with an existing N20
followed by the standard 83-nt SpCas9 scaffold used by the protocol.
Primer-binding sequences are taken
from the two ends of that scaffold, independently of any sequence between the
scaffold and `site2`. Non-template 5′ assembly tails are added separately and
are never searched for in the input plasmid. The complete binding pair must
define exactly one scaffold PCR product on the circular pTarget.
Because long 5-prime tails do not anneal during the initial cycle, DECIPHER is
given the post-first-cycle template in which those tails have been incorporated.
The reported location always refers to the original biological template.

WetLab output is created only after homology and screening pairs pass the
non-relaxable gates of the selected level. A fallback is marked in command
output, `design.log`, `report.tsv`, `primer_pair_ranking.tsv`, and
`wet_lab_report.txt`. Failed targets keep their technical trace and `error.txt`;
other targets continue.

Run tests with:

```bash
bash test/test_r_environment.sh
bash tools/run-r test/test_unit.R
bash test/test_run.sh
```

`test_run.sh` prints stage names, statuses, and durations. Full command output is
written to the timestamped log shown at startup; use `--verbose` to mirror it
to the console as well. A successful integration that uses QC fallbacks or an
expected target rejection is reported as `PASS WITH WARNINGS`.

### QC integration baseline

The MG1655 fixture deliberately limits arms to `350/450` nt. Its effective
high-stringency limits are: primer length `18..22`, GC ratio `0.4..0.6`, GC
clamp `1..3`, runs and repeats `0..4`, Tm `55..65 °C`, pair ΔTm `0..5 °C`,
self-dimer ΔG `>= -5`, cross-dimer ΔG `>= -7`, secondary-structure ΔG
`>= -1`, and primer efficiency `>= 0.001`.

`recA` and `hupB` may still exhaust structural homology-arm candidates; the
filtering level intentionally does not change that geometry. The default-mode
integration requires at least one test target to produce a complete design.
`test_screening_fixture.R` independently verifies strict selection, rejection
of Primer3 row 1, fallback warnings, and selection of row 2.
