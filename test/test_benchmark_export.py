"""Synthetic export checks; no external packages or DHOLE installation required."""
import csv
import importlib.util
import random
import tempfile
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location(
    "export_autoesdcas", Path(__file__).resolve().parents[1] / "benchmark/export_autoesdcas.py")
exporter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(exporter)


def table(path, rows, delimiter=","):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


class ExportTest(unittest.TestCase):
    def setUp(self):
        rng = random.Random(4)
        self.genome = "".join(rng.choices("ACGT", k=1200))
        self.report = dict(target="b0001", **{"class": "cds"}, target_start="301", target_end="700",
            deletion_start="401", deletion_end="600", deleted_nt="200",
            effective_bridge="ATGACTGCCCGCAA", genomic_bridge="ATGACTGCCCGCAA")

    def model(self, report):
        return self.genome[:400] + report["genomic_bridge"] + self.genome[600:]

    def test_coordinates_bridge_and_invalid_models(self):
        for strand in ("+", "-"):
            report = dict(self.report)
            if strand == "-":
                report["genomic_bridge"] = exporter.reverse_complement(report["effective_bridge"])
            result = exporter.extract_deletion(self.genome, report, self.model(report), strand, 200)
            self.assertEqual(result, (401, 600, self.genome[200:400], self.genome[400:600], report["genomic_bridge"]))
        for updates in ({"deleted_nt": "199"}, {"deletion_start": "0"}, {"genomic_bridge": "AAA"}):
            with self.assertRaises(ValueError):
                exporter.extract_deletion(self.genome, dict(self.report, **updates), self.model(self.report), "+", 200)
        with self.assertRaisesRegex(ValueError, "Edited genome mismatch"):
            exporter.extract_deletion(self.genome, self.report, "N" + self.model(self.report)[1:], "+", 200)
        with self.assertRaisesRegex(ValueError, "not unique"):
            genome = "A" * 1200
            exporter.extract_deletion(genome, self.report, genome[:400] + self.report["genomic_bridge"] + genome[600:], "+", 200)

    def test_complete_export_keeps_repetitions_and_skips_failure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            genome = root / "genome.fna"
            genome.write_text(">chr\n" + self.genome + "\n", encoding="utf-8")
            rows = []
            for run, strand in ((1, "+"), (2, "-")):
                attempt = f"run_{run:03d}_gene_0001"
                tech = root / "raw" / attempt / "TechReport"
                target = tech / "b0001_results"
                target.mkdir(parents=True)
                report = dict(self.report)
                if strand == "-":
                    report["genomic_bridge"] = exporter.reverse_complement(report["effective_bridge"])
                (target / "report.tsv").write_text("".join(f"{k}\t{v}\n" for k, v in report.items()), encoding="utf-8")
                (target / "edited_genome.fasta").write_text(">edited\n" + self.model(report) + "\n", encoding="utf-8")
                table(tech / "design_summary.tsv", [dict(gene="b0001", status="ok", output_dir="/old/path/b0001_results")], "\t")
                rows.append(dict(attempt=attempt, run=run, target="b0001", gene="test", contig="chr", strand=strand,
                    gene_nt=400, status="ok", deletion_start=401, deletion_end=600, deleted_nt=200))
            rows.append(dict(rows[0], attempt="run_003_gene_0001", status="error"))
            table(root / "timings.csv", rows)
            output = root / "autoesdcas"
            exporter.export(root, genome, output)
            deletion = exporter.read_table(output / "deletion.csv")
            bridge = exporter.read_table(output / "bridge.csv")
            manifest = exporter.read_table(output / "manifest.csv")
            self.assertEqual(len(deletion), 2)
            self.assertEqual(len(manifest), 4)
            self.assertEqual(list(deletion[0]), exporter.HEADER)
            self.assertEqual(deletion[0]["Inserted sequence"], "-")
            self.assertEqual(deletion[0]["Reference sequence"], self.genome[400:600])
            self.assertEqual(bridge[1]["Inserted sequence"], exporter.reverse_complement(self.report["effective_bridge"]))
            self.assertEqual(bridge[0]["Manipulation type"], "substitution")
            with self.assertRaisesRegex(ValueError, "already exists"):
                exporter.export(root, genome, output)
            (root / "raw/run_002_gene_0001/TechReport/b0001_results/edited_genome.fasta").write_text(">wrong\nACGT\n")
            with self.assertRaisesRegex(ValueError, "Edited genome mismatch"):
                exporter.export(root, genome, root / "invalid")
            self.assertFalse((root / "invalid").exists())


if __name__ == "__main__":
    unittest.main()
