#!/usr/bin/env bash
set -uo pipefail

usage() {
  cat <<'EOF'
Usage:
  old_scheme/po_genam.sh --output_dir DIR --genome FASTA --genome_annotation TSV_OR_GFF \
    --target_plasmid FASTA --gene_list FILE [--cas_plasmid FASTA] [--annotation_format bakta|gff]
EOF
}

output_dir=""; genome=""; genome_annotation=""; target_plasmid=""; gene_list=""
cas_plasmid=""; annotation_format="bakta"
while [ "$#" -gt 0 ]; do
  case "$1" in
    --output_dir|--output-dir) output_dir="$2" ;;
    --genome) genome="$2" ;;
    --genome_annotation|--genome-annotation) genome_annotation="$2" ;;
    --target_plasmid|--target-plasmid) target_plasmid="$2" ;;
    --gene_list|--gene-list) gene_list="$2" ;;
    --cas_plasmid|--cas-plasmid) cas_plasmid="$2" ;;
    --annotation_format|--annotation-format) annotation_format="$2" ;;
    --help|-h) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
  [ "$#" -ge 2 ] || { echo "Missing value for $1" >&2; exit 2; }
  shift 2
done
for name in output_dir genome genome_annotation target_plasmid gene_list; do
  [ -n "${!name}" ] || { echo "Missing --${name}" >&2; usage >&2; exit 2; }
done

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$output_dir"
exec > >(tee -a "$output_dir/batch.log") 2>&1
while IFS= read -r gene || [ -n "$gene" ]; do
  [ -n "$gene" ] || continue
  args=(--genome "$genome" --genome_annotation "$genome_annotation" --gene_name "$gene" --target_plasmid "$target_plasmid" --output_dir "$output_dir/${gene,,}_results" --annotation_format "$annotation_format")
  [ -n "$cas_plasmid" ] && args+=(--cas_plasmid "$cas_plasmid")
  "$script_dir/chopa.sh" "${args[@]}"
done < "$gene_list"
