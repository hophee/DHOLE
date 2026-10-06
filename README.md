# DHOLE

DHOLE designs CRISPR-Cas9 N20 oligos, homology-arm PCR primers, screening
primers, and edited-genome and pTarget models for bacterial CDS and ncRNA targets.

## Installation

With Conda and Git available, run from the project directory:

```bash
bash install.sh
```

The installer creates the `oligo_design` and `oligo_design_chopchop`
environments and installs the required tools and R packages.
Use `bash tools/run-r` to run the project's R scripts.

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

Provide a single-contig genome FASTA, a Bakta TSV (`--annotation-format bakta`)
or GFF/GFF3 (`--annotation-format gff`) annotation, and pTarget and pCas FASTA
files. The genome is treated as linear and plasmids as circular. pTarget must
contain exactly one record for edited-plasmid modelling.

Select targets with `--cds`, `--ncrna`, or both. Targets are matched first by
`locus_tag`, then by `gene`; lists may be comma-separated or space-separated.

Each restriction site must occur once in pTarget. The selected cassette must
start with an existing N20 followed by the protocol's 83-nt SpCas9 scaffold.
`--ptarget-cassette-arc` selects the shortest arc by default; `forward` follows
the input strand from `site1` to `site2`, and `reverse` follows the
reverse-complement strand.

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

For all available options, run:

```bash
bash tools/run-r oligo_designer.R --help
```

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

`WetLab/<target>_results/` contains the final oligos, edited genome and pTarget
models, predicted PCR products, and a wet-lab report. One edited pTarget and
one sgRNA-cassette PCR product are provided per selected N20. PCR outputs also
include both homology arms and screening products for the original and edited
genomes. The report lists primer Tm values, product lengths, N20-to-arm
distances, and QC results.

`TechReport/` contains the design summary and run parameters. Per-target tables
record primer binding sites, predicted amplicons, QC metrics, candidate rankings,
and selection decisions; `design.log` records processing details.

WetLab files are produced for successful targets. Designs using fallback
candidates are marked in the reports. Failed targets retain technical reports
and `error.txt`; processing continues for the remaining targets.

Run tests with:

```bash
bash test/test_r_environment.sh
bash tools/run-r test/test_unit.R
bash test/test_run.sh
```

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



## Benchmark

Run `bash benchmark/run.sh [N=10] [M=10] [seed=1]` for repeated default-level-2
designs on random MG1655 CDS. Runtime and primer/filtering plots are saved to
a PDF; DHOLE outputs are archived automatically. See [benchmark/README.md](benchmark/README.md).

## License

DHOLE's original project code, including its scripts, tests, and configuration
files, is licensed under the GNU General Public License as published by the
Free Software Foundation, either version 2 of the License, or (at your option)
any later version (SPDX: `GPL-2.0-or-later`). See [LICENSE](LICENSE) for the
full GPL version 2 text.
