#!/usr/bin/env python3
"""Export successful benchmark deletions in AutoESDCas' upstream-sequence format."""
import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

HEADER = ["Name", "Sequence upstream of the manipulation site (>100bp)",
          "Reference sequence", "Inserted sequence", "Manipulation type"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def read_table(path, delimiter=","):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def read_report(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    require(all(len(row) == 2 for row in rows), f"Malformed report: {path}")
    require(len({row[0] for row in rows}) == len(rows), f"Duplicate report keys: {path}")
    return dict(rows)


def read_fasta(path):
    records = {}
    key = None
    with path.open(encoding="utf-8-sig") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                key = line[1:].split()[0]
                require(key not in records, f"Duplicate FASTA ID: {key}")
                records[key] = []
            else:
                require(key is not None, f"Sequence before FASTA header: {path}")
                records[key].append(line.upper())
    result = {key: "".join(parts) for key, parts in records.items()}
    require(result and all(result.values()), f"Empty FASTA: {path}")
    return result


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def extract_deletion(genome, report, edited, strand, upstream_length):
    start, end = int(report["deletion_start"]), int(report["deletion_end"])
    require(1 <= start <= end <= len(genome), "Deletion outside genome")
    require(end - start + 1 == int(report["deleted_nt"]), "deleted_nt mismatch")
    require(strand in ("+", "-"), "Invalid target strand")
    bridge = report["genomic_bridge"].upper()
    expected = report["effective_bridge"].upper()
    if strand == "-":
        expected = reverse_complement(expected)
    require(bridge == expected, "Bridge orientation mismatch")
    require(edited == genome[:start - 1] + bridge + genome[end:], "Edited genome mismatch")
    require(upstream_length > 100, "Upstream sequence must be >100 bp")
    require(start - 1 >= upstream_length, "Insufficient upstream sequence; no circular wrapping")
    upstream = genome[start - 1 - upstream_length:start - 1]
    reference = genome[start - 1:end]
    require(set(upstream + reference + bridge) <= set("ACGT"), "Non-ACGT bases in exported sequences")
    hits = 0
    for query in (upstream, reverse_complement(upstream)):
        position = genome.find(query)
        while position >= 0:
            hits += 1
            position = genome.find(query, position + 1)
    require(hits == 1, "Upstream sequence is not unique; increase --upstream")
    return start, end, upstream, reference, bridge


def export(results, genome_path, output, upstream_length=200):
    require(not output.exists(), f"Output already exists: {output}")
    genomes = read_fasta(genome_path)
    genome_md5 = hashlib.md5(genome_path.read_bytes()).hexdigest()
    require(len(genomes) == 1, "Expected the single-contig MG1655 FASTA")
    attempts = read_table(results / "timings.csv")
    require(attempts, "Empty timings.csv")
    require(len({row["attempt"] for row in attempts}) == len(attempts), "Duplicate attempt IDs")
    exports = {"deletion": [], "bridge": []}
    manifest = []
    for row in attempts:
        if row["status"] != "ok":
            print(f"Skipping {row['attempt']}: {row['status']}", file=sys.stderr)
            continue
        require(re.fullmatch(r"run_\d+_gene_\d+", row["attempt"]), "Invalid attempt ID")
        run = results / "raw" / row["attempt"]
        summary = read_table(run / "TechReport/design_summary.tsv", "\t")
        require(len(summary) == 1 and summary[0]["status"] == "ok" and
                summary[0]["gene"] == row["target"], f"Summary mismatch: {run}")
        folder = summary[0]["output_dir"].replace("\\", "/").rstrip("/").split("/")[-1]
        require(folder not in ("", ".", ".."), "Invalid target directory")
        tech = run / "TechReport" / folder
        report = read_report(tech / "report.tsv")
        require(report["target"] == row["target"] and report["class"] == "cds", "Target mismatch")
        require(int(report["target_end"]) - int(report["target_start"]) + 1 == int(row["gene_nt"]),
                "Gene length mismatch")
        require(row["contig"] in genomes, f"Unknown contig: {row['contig']}")
        models = read_fasta(tech / "edited_genome.fasta")
        require(len(models) == 1, "Expected one edited genome")
        start, end, upstream, reference, bridge = extract_deletion(
            genomes[row["contig"]], report, next(iter(models.values())), row["strand"], upstream_length)
        require(start == int(row["deletion_start"]) and end == int(row["deletion_end"]) and
                len(reference) == int(row["deleted_nt"]), "Timing/report deletion mismatch")
        for variant in exports:
            name = f"{row['attempt']}_{variant}"
            inserted, operation = ("-", "deletion") if variant == "deletion" else (bridge, "substitution")
            exports[variant].append([name, upstream, reference, inserted, operation])
            manifest.append(dict(Name=name, variant=variant, attempt=row["attempt"], run=row["run"],
                target=row["target"], gene=row["gene"], contig=row["contig"], strand=row["strand"],
                gene_nt=row["gene_nt"], deletion_start=start, deletion_end=end, deleted_nt=len(reference),
                genomic_bridge=bridge, genome_md5=genome_md5,
                result_dir=f"raw/{row['attempt']}", target_dir=folder))
    require(manifest, "No successful designs to export")
    # Validate all results before writing; never overwrite a previous export.
    output.mkdir(parents=True)
    for variant, rows in exports.items():
        with (output / f"{variant}.csv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(rows)
    with (output / "manifest.csv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest)
    print(f"Exported {len(manifest) // 2} designs in each variant to {output}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path, help="directory produced by run_dhole.R")
    parser.add_argument("--genome", type=Path,
                        default=Path(__file__).resolve().parents[1] / "local/data/MG1655.fna")
    parser.add_argument("--output", type=Path, help="new directory (default: RESULTS/autoesdcas)")
    parser.add_argument("--upstream", type=int, default=200)
    args = parser.parse_args()
    export(args.results, args.genome, args.output or args.results / "autoesdcas", args.upstream)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"ERROR: {error}")
