# DHOLE

DHOLE designs CRISPR-Cas9 N20 oligos, homology-arm PCR primers, screening
primers, and edited-genome and pTarget models for bacterial CDS and ncRNA targets.

## Installation

Use a Linux environment with Bash, Conda, Git, `curl`, `gzip`,
and GNU build tools available. The installer invokes `make` to compile Primer3;
the build tools and C/C++ compiler must be available in the shell running it.
Run from the project directory:

```bash
bash install.sh
```

The installer creates the `oligo_design` and `oligo_design_chopchop`
environments, or updates them with `conda env update --prune` if they already
exist. It installs R packages, builds Primer3, prepares CHOPCHOP and the tool
wrappers, and downloads the MG1655 genome and annotation used by the tests.
Use `bash tools/run-r` to run the project's R scripts.

## Usage

```bash
bash tools/run-r oligo_designer.R \
  --genome local/data/genome.fasta \
  --genome-annotation local/data/genome.gff \
  --annotation-format gff \
  --target-plasmid local/data/pTarget.fasta \
  --filtering-level 2 \
  --site1 ACTAGT \
  --site2 CTGCAG \
  --cas-plasmid local/data/pCas.fasta \
  --output-dir local/runs/results \
  --cds gene1,gene2 \
  --ncrna rna1,rna2
```

Provide a single-contig genome FASTA, a Bakta TSV (`--annotation-format bakta`)
or GFF/GFF3 (`--annotation-format gff`) annotation, and a pTarget FASTA.
`--cas-plasmid` is optional. When supplied, its sequences are included in primer
specificity checks. The genome is treated as linear and plasmids as circular.
pTarget must contain exactly one record for edited-plasmid modelling.

### Annotation and targets

Bakta TSV input must provide `start`, `stop`, `strand`, `gene`, `locus_tag`,
and `sequence_id` or `seqid` after column-name normalization. The `type` column
is optional. Both name columns must exist, even if individual values are empty.
For example, a tab-separated header and record are:

```text
sequence_id	type	start	stop	strand	locus_tag	gene
chr1	CDS	101	1000	+	LOC_0001	gene1
```

GFF/GFF3 input uses the standard nine columns. DHOLE reads `seqid`, `type`,
`start`, `end`, `strand`, and the `gene` and `locus_tag` attributes; it does not
use `ID` or `Name` to look up targets. At least one of `gene` or `locus_tag`
must identify each requested target. For example:

```text
chr1	Bakta	CDS	101	1000	.	+	0	locus_tag=LOC_0001;gene=gene1
```

Coordinates are 1-based and inclusive, with `+` or `-` for the target strand.
Annotation sequence IDs must match the genome FASTA record ID. GFF files may
contain an embedded FASTA section; DHOLE reads genome sequences from `--genome`.

Select targets with `--gene` (automatic), `--cds`, or `--ncrna`. Targets are matched first by
`locus_tag`, then by `gene`; lists may be comma-separated or space-separated.
`--gene` uses the structured annotation `type`, ignoring case and surrounding
whitespace: `CDS` selects CDS geometry, any other type selects ncRNA geometry.
If the entire annotation lacks type information, all automatic targets use CDS.
An empty or unknown type in an otherwise typed annotation uses ncRNA and is
reported. Explicit options override `--gene`; aliases resolving to the same
record are designed once. Conflicting `--cds`/`--ncrna` assignments are errors.

`--n20-strands plus` and `minus` refer to the genomic reference strands,
regardless of the target's strand. With more than one requested N20, `both`
requires at least one guide on each strand; with one N20, either strand is
accepted. `random` leaves the strand unconstrained. Guide sets are enumerated
deterministically from the ranked candidates; this mode does not draw a random
sample.

Each restriction site must occur once in pTarget. The selected cassette must
start with its original N20 followed by the selected `--sgrna-scaffold`.
The default scaffold is 83 nt; alternatives may differ in sequence
and length. Choose a scaffold compatible with your Cas9: DHOLE checks sequence
and assembly consistency, not experimental activity.
`--ptarget-cassette-arc` selects the shortest arc by default; `forward` follows
the input strand from `site1` to `site2`, and `reverse` follows the
reverse-complement strand.

### Design parameters

| Argument | Default | Meaning |
|---|---:|---|
| `--n20-mn` | `1` | Required N20 count |
| `--n20-strands` | `random` | `plus`, `minus`, `both`, or unconstrained `random` |
| `--n20-offtarget` | `0` | Maximum CHOPCHOP `MM0,MM1,...` values |
| `--n20-mid-closeness-max` | `0.18` | CDS candidates require `mid_closeness <= threshold`; `0` disables this filter. Does not filter ncRNA |
| `--site1` | `ACTAGT` (SpeI) | First restriction-site sequence in insert orientation |
| `--site2` | `CTGCAG` (PstI) | Second restriction-site sequence in insert orientation |
| `--ptarget-cassette-arc` | `shortest` | Arc between the sites used as the cassette: `shortest`, `forward`, or `reverse` |
| `--bridge-sequence` | `ATGACTGCCCGCAAG` | Bridge in target orientation, at least 3 nt; CDS uses a phase-preserving prefix, ncRNA uses the full sequence |
| `--sgrna-scaffold` | 83-nt sequence below | Constant sgRNA DNA sequence without N20; must match the selected pTarget cassette after its original N20 |
| `--sgrna-forward-annealing` | First 36 scaffold nt | Forward primer annealing region, written 5′→3′ |
| `--sgrna-reverse-annealing` | Reverse complement of last 23 scaffold nt | Reverse primer annealing region, written 5′→3′ |
| `--sgrna-reverse-overhang` | `AGTTGACGCT` | Reverse primer 5′ tail; its reverse complement supplies the left-arm overlap |
| `--sgrna-annealing-temp-c` | `60` | Annealing temperature for the sgRNA-cassette PCR |
| `--cds-fs`, `--ncrna-fs` | off | Require deleted length divisible by three |
| `--left-arm-min/opt/max` | `300/350/400` | Left-arm structural limits |
| `--right-arm-min/opt/max` | `400/450/500` | Right-arm structural limits |
| `--n20-arm-min-distance` | `40` | Minimum N20-to-arm distance in nt |
| `--ncrna-neighbor-buffer` | `20` | Minimum retained gap between the deletion and neighboring annotations |
| `--ncrna-max-n20-attempts` | `3` | Maximum N20 sets tried for ncRNA, including the initial set |

Service sequences accept only A/C/G/T and are normalized to uppercase. The
default scaffold is:

```text
GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTT
```

For a scaffold shorter than 36 or 23 nt, explicitly supply the corresponding
annealing sequence. The chosen primers must uniquely amplify the entire
selected scaffold from circular pTarget. DHOLE selects the new N20 from the
genomic target candidates.

### Assembly and oligonucleotide names

The insert, including its restriction sites, has this order in target orientation:

```text
site1 | N20 | sgRNA scaffold | left-arm overlap | left arm | bridge | right arm | site2
```

Each N20 produces a separate pTarget model with the same homology arms.
The overlap between the sgRNA PCR product and the left arm is the reverse
complement of `--sgrna-reverse-overhang`. The bridge joins the two arms.
Shared overlaps occur once in the assembled sequence. Genome editing replaces
the deleted interval with the bridge; the sgRNA cassette stays in pTarget.

Primer names use the target's annotated gene name when available, followed by
an underscore and one of these suffixes. All oligonucleotide sequences are
written 5′ to 3′.

| Suffix | Purpose | Full oligonucleotide |
|---|---|---|
| `sgF<N>` | Forward sgRNA-cassette primer for N20 number N | `ACG` + site1 + N20 + forward scaffold annealing region |
| `sgR` | Shared reverse sgRNA-cassette primer | Reverse overhang + reverse scaffold annealing region |
| `LF` | Left-arm forward primer | Left-arm overlap + annealing region |
| `LR` | Left-arm reverse primer | Reverse complement of bridge + annealing region |
| `RF` | Right-arm forward primer | Bridge + annealing region |
| `RR` | Right-arm reverse primer | `ACG` + reverse complement of site2 + annealing region |
| `scrF`, `scrR` | Genome-screening primers | Annealing region only |

The annealing region binds the PCR template. The 5′ tails supply the guide,
restriction sites, or overlaps. Left and right arms follow target orientation;
for a minus-strand target their order is reversed relative to genomic coordinates.

For CDS, a bridge of length `L` and deletion of length `d` uses the first
`L - ((L - d) mod 3)` bases, removing at most two bases from its 3′ end. This
preserves `(inserted length - deleted length) mod 3 = 0`. The default bridge
therefore gives lengths 15, 13, and 14 for deletion remainders 0, 1, and 2.
The bridge is reverse-complemented when inserted into the genomic FASTA for
a minus-strand target. `--cds-fs` and `--ncrna-fs` independently impose their
requirement that the deletion length be divisible by three. Junction codons
and possible stop codons are not evaluated.

### CDS search

CDS designs delete an internal part of the annotated feature. The inner ends
of the selected homology arms must lie within the CDS, and the deleted length
must be smaller than the feature length. The outer ends may extend beyond it.

The initial search interval depends on CDS length:

| CDS length | Initial interval |
|---|---|
| Up to 500 nt | Entire annotated CDS |
| 501 to 1499 nt | 501-nt interval centred at `start + length %/% 2` |
| At least 1500 nt | Central interval after removing `floor((end - start) / 3)` nt from each end |

CHOPCHOP searches this interval once. For each N20 set, the arm search starts
with templates of the preferred arm lengths on either side. Each further
attempt moves both interval boundaries inward by one nucleotide and increases
the template lengths up to the configured maxima. Primer3 may trim these
templates; DHOLE checks arm lengths and N20 distances on the resulting products.
The search stops when a pair is accepted or the next interval is shorter than
30% of the CDS length. This interval cutoff is not a minimum final deletion size.
Candidate selection and fallback rules are described under [Quality control](#quality-control).

`mid_closeness` is the absolute distance between the integer midpoint of
the feature and the integer midpoint of N20, divided by feature length.
Candidates are sorted by this metric even when its filter is disabled.
`--n20-arm-min-distance` independently controls the N20-to-arm offset.

### ncRNA search

ncRNA search starts with deletion of the entire annotated target, then expands
the available boundary positions outward. Previously generated arm candidates
remain available; either side can continue when the other reaches its limit.
Candidates are tried in increasing total extension (left plus right), with the
same QC gates as CDS. Within an extension, the QC ranking applies.
The inclusive deletion bounds are `left_neighbor_end + buffer + 1` and
`right_neighbor_start - buffer - 1`. Neighbors on either strand of the same
contig count; missing neighbors use the contig edge without a buffer. Overlaps
with other features or insufficient gaps fail explicitly. Region/contig
container records are excluded. Arms may include retained neighboring genes.

CHOPCHOP runs once over the largest permitted deletion interval. Each N20 set
is searched through all allowed boundary extensions before moving to the next
set, up to the attempt limit. Primer3 candidates and specificity and
openPrimeR caches are reused. Output or screening QC failures also try further
arm candidates for that set.

For all available options, run:

```bash
bash tools/run-r oligo_designer.R --help
```

## Quality control

### Primer generation and strict QC

Primer3 generates candidates with defaults of `18/21/27` nt for minimum,
preferred, and maximum primer length, Tm `55/60/65 °C`, GC `30..70%`, and pair
ΔTm up to `8 °C`. These generation settings apply at every filtering level.

The subsequent openPrimeR strict profile uses primer length `18..22`, GC ratio
`0.4..0.6`, GC clamp `0..3`, runs and repeats `0..4`, Tm `55..65 °C`, pair ΔTm
`0..5 °C`, self-dimer ΔG `>= -5`, cross-dimer ΔG `>= -7`, secondary-structure
ΔG `>= -1`, and primer efficiency `>= 0.001`. A GC clamp of zero passes this
constraint but produces a `low_gc_clamp` warning and is a late ranking
criterion. Such a warning can occur in a strict-passing design without fallback.

### Specificity and filtering parameters

| Argument | Default | Meaning |
|---|---:|---|
| `--filtering-level` | `2` | `1` lite/report, `2` default/warn, `3` hard/core |
| `--primer-max-mismatches` | `2` | Maximum mismatches per binding site |
| `--primer-critical-3p-bases` | `5` | Critical primer 3′ region |
| `--primer-max-3p-mismatches` | `0` | Allowed mismatches in that region |
| `--primer-min-product-size` | `50` | Minimum counted amplicon size |
| `--primer-max-product-size` | `2000` | Maximum counted amplicon size |
| `--primer-max-offtarget-products` | `0` | Preferred maximum non-intended amplicons |
| `--primer-hard-max-tm-diff` | `5` | Maximum pair ΔTm in °C required by level 3 |

### Selection policy

Arm geometry, feature bounds, N20 distance, frame, deletion rules, and Primer3
generation are identical at all filtering levels. The levels affect only PCR
specificity and openPrimeR QC:

1. `lite` (report only). Specificity and openPrimeR metrics rank candidates but
   do not block a structurally valid Primer3 pair.
2. `default` (warn and continue). Exactly one intended PCR product is required;
   off-target and openPrimeR failures are ranked risks and may be returned with
   warnings.
3. `hard` (core gates). Exactly one intended product, zero high-risk off-target
   products with perfect primer 3′ ends, and a finite pair ΔTm no greater than
   `--primer-hard-max-tm-diff` (default `5 °C`) are mandatory. Other openPrimeR
   and specificity criteria are optimization criteria and may appear in a
   warned fallback.

If openPrimeR evaluation is unavailable, DHOLE records a warning and uses the
Primer3 pair ΔTm for the level-3 check. A hard-mode design can therefore pass
the core gates without a completed openPrimeR evaluation. Changing the hard
ΔTm threshold does not change the openPrimeR strict-profile limit.

Within the candidates being compared, a pair passing both specificity and the
complete openPrimeR strict profile is preferred. Other eligible pairs are
ranked by intended-product availability, off-target risk, openPrimeR penalty,
dimer risk, Tm difference, Primer3 penalty, deleted nucleotides, low-GC-clamp
count, and original Primer3 order. Core gates cannot be compensated by a score.

For a CDS N20 set, `lite` accepts a pair from the first search window with an
eligible candidate, even if that pair requires fallback. Levels 2 and 3 defer
fallback candidates while checking later windows for a strict-passing pair.
If no strict pair is found, they select the best deferred candidate for that
N20 set. Search ends at the first N20 set that completes a design; strict and
fallback candidates are not compared across all N20 sets.

ncRNA prioritizes the smallest total extension passing these same mandatory
gates; strict-profile preference and risk ranking apply within that extension.

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

### File formats and reports

| File | Format and use |
|---|---|
| `final_sequences.fasta` | Named full oligonucleotide sequences, including 5′ tails |
| `final_sequences.txt` | Tab-separated `name`, `sequence`, and `purpose` columns with a header; the same oligos as the FASTA, suitable for preparing an order |
| `edited_genome.fasta`, `edited_pTargets.fasta` | Sequence models for inspecting the genomic edit and plasmid constructs |
| `pcr_products.fasta`, `pcr_products.tsv` | Product sequences; the TSV also includes names, descriptions, locations, lengths, primers, and PCR conditions |
| `wet_lab_report.txt` | Human-readable design and PCR details, including primer metrics and QC warnings |
| `design_summary.tsv` | Tabular per-target status, failure stage/reason, and result paths |
| `run_parameters.tsv` | Tab-separated `parameter` and `value` columns with a header |
| `<target>_results/report.tsv` | Two tab-separated fields per line, parameter name and value, without a header |
| Primer binding, amplicon, QC, and ranking TSVs | Tables with headers for examining candidate checks and selection decisions |

`TechReport/` contains the design summary and run parameters. Per-target tables
record primer binding sites, predicted amplicons, QC metrics, candidate rankings,
and selection decisions; `design.log` records processing details.
`target_strategies.tsv` records resolved targets, annotation types, and selection
reasons.

`run_parameters.tsv` records the original bridge and all resolved sgRNA
sequences, including annealing regions and assembly overlap. A completed
target's `report.tsv` records `effective_bridge` in target orientation,
`effective_bridge_length`, `bridge_trimmed_nt`, and `genomic_bridge` in genomic
FASTA orientation.

After feature lookup, `report.tsv` records the target bounds, annotation type,
and chosen strategy. ncRNA reports also record allowed deletion bounds,
neighbors, and buffer. Actual deletion bounds are added after successful
selection, along with per-side extension for ncRNA. N20 search records its
attempt count and termination reason; before it starts, these are zero and
`not_started`. Failed targets retain the diagnostics available at the failure
stage and the failure reason. An early lookup error has no feature or deletion
bounds. ncRNA expansion windows and reached limits are recorded in `design.log`.

WetLab files are produced for successful targets. Designs using fallback
candidates are marked in the reports. Failed targets retain technical reports
and `error.txt`; processing continues for the remaining targets.

### Unique primer collection

Each completed run writes `WetLab/unique_primers.tsv`, combining oligonucleotides
from successful targets into one row per uppercase sequence. This file is
created even without `--primer-table`; it contains only the header when there
are no successful targets and no imported entries.

To extend a collection, add `--primer-table previous_results/WetLab/unique_primers.tsv`
to the next design command. DHOLE merges the imported entries with the new
oligos and writes the combined table into the new output directory. Reimporting
the same entries does not duplicate sequences or source pairs.

| Required column | Contents |
|---|---|
| `sequence` | Full oligonucleotide sequence, used for deduplication |
| `name` | Distinct oligonucleotide names, separated by semicolons |
| `source_gene` | Original requested target identifiers |
| `genome` | Source genome FASTA filename without its final extension |
| `gene` | Target labels; a target found in multiple source genomes is labelled `source_gene_genome` |

The semicolon-separated entries in `source_gene`, `genome`, and `gene` correspond
position by position. Input tables must include all five columns, nonempty
values, valid DNA sequences, and consistent source lists. `final_sequences.txt`
has a different schema and cannot be passed directly as `--primer-table`.

### Reusing an output directory

Genome indexes are reused when the FASTA checksum and Bowtie version marker
match and all index files are present and nonempty. Otherwise DHOLE rebuilds
them. Target directory names are lowercased, with characters outside
`A-Z`, `a-z`, `0-9`, `_`, `.`, and `-` replaced by underscores.

A rerun removes the seven listed WetLab files for each requested target before
designing it again. It resets that target's `design.log` and `report.tsv` and
removes its previous `error.txt`. The whole output directory is not cleared:
other technical files and directories for targets omitted from the new run can
remain. Use a new output directory for a separate run, or consult the current
`design_summary.tsv` to distinguish its results from older files. Use
`--primer-table` explicitly to carry a primer collection into a new run.

## Testing

Run tests with:

```bash
bash test/test_r_environment.sh
bash tools/run-r test/test_unit.R
bash tools/run-r test/test_restriction_sites.R
bash tools/run-r test/test_ncrna.R
bash test/test_run.sh
```

### QC integration baseline

The MG1655 fixture limits the left arm to `300..350` nt and the right arm to
`400..450` nt. It uses the generation settings and strict QC profile described
under [Quality control](#quality-control), with filtering level 2.

`recA` and `hupB` may exhaust structural homology-arm candidates; the
filtering level intentionally does not change that geometry. The default-mode
integration requires at least one test target to produce a complete design.
`test_screening_fixture.R` independently verifies strict selection, rejection
of Primer3 row 1, fallback warnings, and selection of row 2.

## License

DHOLE's original project code, including its scripts, tests, and configuration
files, is licensed under the GNU General Public License as published by the
Free Software Foundation, either version 2 of the License, or (at your option)
any later version (SPDX: `GPL-2.0-or-later`). See [LICENSE](LICENSE) for the
full GPL version 2 text.

`callPrimer3.R` is derived from third-party code. Its file header states that
redistribution permission is unresolved and excludes it from DHOLE's
`GPL-2.0-or-later` grant. That exception applies to the wrapper even though it
is included in this repository. Downloaded dependencies, including Primer3
and CHOPCHOP, retain their own licenses.
