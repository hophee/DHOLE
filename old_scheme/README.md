# 2PAC legacy scheme

## Installation

Use `../install.sh` to install most dependencies. The script is safe to rerun:
it updates the `chopchop` Conda environment and reuses existing Git checkouts.

The CHOPCHOP upstream repository is cloned from
[GitHub](https://github.com/JokingHero/chopchop).

You will also need to install the following R packages: `dplyr`, `readr`,
`Biostrings`, `ape`.

## Usage

```bash
./old_scheme/chopa.sh \
  --genome genome.fasta \
  --genome_annotation genome_annotation.tsv \
  --gene_name gene_name \
  --target_plasmid pTarget.fasta \
  --output_dir output_directory
```

- Bakta TSV and GFF/GFF3 annotation are supported through
  `--annotation_format bakta|gff`.
- The target is looked up first in `locus_tag`, then in `gene`.

## Output

- `all_primers.fasta`: set of designed primers
- `edited_genome.fasta` and `pTargets.fasta`: genome and pTarget sequences after editing
- `n20_table.tsv`: chopchop.py output
- `report.txt`: melting temperatures and PCR product lengths
- `PCR_product_sequnces.fasta`: predicted PCR products
