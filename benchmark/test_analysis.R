# Run from the repository root: bash tools/run-r benchmark/test_analysis.R
source("benchmark/benchmark.R")

local({
  out <- tempfile("dhole-benchmark-test-")
  dir.create(out)
  stopifnot(startsWith(normalizePath(out), normalizePath(tempdir())))
  on.exit(unlink(out, recursive = TRUE))
  gff <- file.path(out, "test.gff")
  writeLines(c("##gff-version 3",
    "chr\ttest\tCDS\t1\t30\t.\t+\t0\tgene=a;locus_tag=b0001",
    "chr\ttest\tCDS\t40\t60\t.\t-\t0\tgene=a;locus_tag=b0001",
    "chr\ttest\tCDS\t70\t90\t.\t+\t0\tgene=b",
    "chr\ttest\trRNA\t100\t120\t.\t+\t.\tlocus_tag=rna",
    "##FASTA", ">chr", "ACGT"), gff)
  stopifnot(identical(cds_ids(gff), c("b0001", "b")))

  tech <- file.path(out, "raw", "run_001", "TechReport")
  target <- file.path(tech, "a_results")
  wet <- file.path(out, "raw", "run_001", "WetLab", "a_results")
  dir.create(target, recursive = TRUE)
  dir.create(wet, recursive = TRUE)
  write_tsv(data.frame(gene = c("a", "b"), status = c("ok", "error"),
    stage = c(NA, "homology"), reason = c(NA, "no candidates"),
    output_dir = c(target, NA), wet_lab_dir = c(wet, NA)), file.path(tech, "design_summary.tsv"))
  writeLines(c("Report", "name\tpurpose\tannealing_sequence\ttm_c",
    "LF\tleft\tGGAA\t60", "LR\tleft\tCCAA\t61",
    "RF\tright\tAAAA\t62", "RR\tright\tCCCC\t63",
    "SF\tscreen\tATGC\t64", "SR\tscreen\tGCGC\t65", "", "Other sections"),
    file.path(wet, "wet_lab_report.txt"))
  write_tsv(data.frame(pair_id = c("p1", "p2"), structure_passed = c(TRUE, FALSE),
    selected = c(TRUE, FALSE), rejection_reason = c("", "structural_constraints")),
    file.path(target, "primer_pair_ranking.tsv"))
  sampled <- data.frame(run = c(1L, 1L, 2L, 2L), gene = c("a", "b", "a", "c"))
  result <- collect_results(file.path(out, "raw"), sampled)
  stopifnot(identical(result$targets$status, c("ok", "error", "missing_result", "missing_result")),
    identical(result$primers$gc_percent, c(50, 50, 0, 100, 50, 100)),
    identical(result$primers$tm_c, 60:65), nrow(result$rankings) == 2L)
  timings <- data.frame(run = 1:2, real_s = c(10, 30), cpu_s = c(20, 60), exit_code = c(0L, 1L))
  analyse(out, sampled, timings)
  means <- read_tsv(file.path(out, "mean_times.tsv"))
  rates <- read_tsv(file.path(out, "success_rates.tsv"))
  stopifnot(means$seconds[means$metric == "Real"] == 20,
    means$seconds[means$metric == "CPU (user + system)"] == 40,
    rates$percent[rates$category == "Success"] == 25,
    file.info(file.path(out, "analysis.pdf"))$size > 1000)
  # A batch with no summary or QC tables must still produce a report.
  analyse(out, sampled[sampled$run == 2L, ], timings[2L, ])
  stopifnot(file.info(file.path(out, "analysis.pdf"))$size > 1000)
})
cat("Benchmark analysis checks passed\n")
