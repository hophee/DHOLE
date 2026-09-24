# 2PAC

2PAC designs CRISPR-Cas9 N20 oligos, homology-arm PCR primers, screening
primers, and an edited-genome model for bacterial CDS and ncRNA targets.

Primer selection is candidate-based: Primer3 proposes up to ten pairs per
homology arm and five screening pairs, then 2PAC applies structural rules,
exhaustive Biostrings specificity, mandatory openPrimeR QC, and a stable
lexicographic ranking. A rejected first Primer3 row no longer rejects the
target. If all rows fail, 2PAC continues with the next admissible N20 set.

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
library. The pipeline fails closed when a required openPrimeR constraint or
executable is unavailable.

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
| `--n20-mn` | `1` | Required N20 count |
| `--n20-strands` | `random` | `plus`, `minus`, `both`, or unconstrained `random` |
| `--n20-offtarget` | `0` | Maximum CHOPCHOP `MM0,MM1,...` values |
| `--site1` | `ACTAGT` (SpeI) | First restriction-site sequence in insert orientation |
| `--site2` | `CTGCAG` (PstI) | Second restriction-site sequence in insert orientation |
| `--cds-fs`, `--ncrna-fs` | off | Require deleted length divisible by three |
| `--left-arm-min/opt/max` | `300/350/400` | Left-arm structural limits |
| `--right-arm-min/opt/max` | `400/450/500` | Right-arm structural limits |
| `--n20-arm-min-distance` | `40` | Minimum N20-to-arm distance in nt |
| `--primer-max-mismatches` | `2` | Maximum mismatches per binding site |
| `--primer-critical-3p-bases` | `5` | Critical primer 3′ region |
| `--primer-max-3p-mismatches` | `0` | Allowed mismatches in that region |
| `--primer-min-product-size` | `50` | Minimum counted amplicon size |
| `--primer-max-product-size` | `2000` | Maximum counted amplicon size |
| `--primer-max-offtarget-products` | `0` | Allowed non-intended amplicons |

The legacy Primer3 generation thresholds remain unchanged: primer length
`18/21/27` nt (min/opt/max), homopolymer maximum `5`, and pair Tm difference
maximum `8 °C`. Candidates outside the stricter openPrimeR profile are retained
in the trace and rejected at its explicit hard gate. Primer3 buffer defaults
remain 50 mM monovalent salt, 1.5 mM Mg, 0.6 mM dNTP, and 50 nM DNA.
`run_parameters.tsv` records these values, the active constraints, their
effective limits, package versions, and tool paths.

## Selection policy

Candidates pass four non-compensating levels:

1. Existing arm geometry, feature bounds, N20 distance, frame, and deletion
   rules.
2. An explicitly identified expected product and no forbidden off-target
   products across genome, pTarget, and pCas. Exhaustive binding-site pairing
   is authoritative; `matchProbePair()` is retained as a comparison generator.
3. Mandatory openPrimeR constraints. Annealing sequence is used for length,
   GC, Tm, efficiency/coverage; the ordered full oligo is used for self/cross
   dimers and secondary structure. Only primers in one physical PCR reaction
   are evaluated together.
4. Deterministic lexicographic ranking by off-target risk, soft failures,
   openPrimeR penalty, dimer risk, Tm difference, Primer3 penalty, deleted
   nucleotides, and original Primer3 order.

No soft score can compensate for a failed hard gate.

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
openPrimeR metric/`EVAL_*` result, gate, rank component, selection flag, and
rejection reason. `design.log` records `primer_qc TRY`, `REJECTED`, and `OK`.
For every successful target, `wet_lab_report.txt` contains the complete final
oligo set, modelled-construction names and lengths, primer Tm values, expected
screening products for edited and unedited alleles, per-N20 distances to both
homology arms, screening off-target counts, and selected openPrimeR quality
metrics with readable labels. `edited_pTargets.fasta` contains one circular
pTarget model for every selected N20. Site matching checks both orientations
and the FASTA origin; each supplied site must occur exactly once physically and
the two sites must share an orientation.

`pcr_products.tsv` and the WetLab report contain one sgRNA-cassette PCR product
per N20, both homology-arm products, and screening products from the original
and edited genomes. Every row records its template location, product length,
full primer names, PCR simulation conditions, and sequence. Products are
generated by `DECIPHER::AmplifyDNA` with the complete service-tailed primers.
Because long 5-prime tails do not anneal during the initial cycle, DECIPHER is
given the post-first-cycle template in which those tails have been incorporated.
The reported location always refers to the original biological template.

WetLab output is created only after homology and screening pairs pass all
gates. Failed targets keep their technical trace and `error.txt`; other targets
continue.

Run tests with:

```bash
bash test/test_r_environment.sh
bash tools/run-r test/test_unit.R
bash test/test_run.sh
```

### Strict-QC integration baseline

The MG1655 fixture deliberately limits arms to `350/450` nt. Its effective
high-stringency limits are: primer length `18..22`, GC ratio `0.4..0.6`, GC
clamp `1..3`, runs and repeats `0..4`, Tm `55..65 °C`, pair ΔTm `0..5 °C`,
self-dimer ΔG `>= -5`, cross-dimer ΔG `>= -7`, secondary-structure ΔG
`>= -1`, and primer efficiency `>= 0.001`.

With this fixture, `recA` and `hupB` currently exhaust structural homology-arm
candidates before openPrimeR; `pta` reaches `primer_qc` but has no pair passing
all hard gates. These are explicit strict-QC baseline outcomes, not successful
primer designs. `test_screening_fixture.R` independently verifies the complete
successful `scrF/scrR` path, including rejection of Primer3 row 1 and selection
of row 2.
