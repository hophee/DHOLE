#!/usr/bin/env Rscript
source("benchmark/run_dhole.R")
source("benchmark/compare_primers.R")

expect_error <- function(expr, pattern) {
  error <- tryCatch({ force(expr); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(pattern, conditionMessage(error)))
}

stopifnot(identical(benchmark_counts(character()), c(10L, 10L, 1L)))
for (bad in c("0", "-1", "2.5", "NA", "2147483648")) expect_error(benchmark_counts(bad), "positive integers")
annotation <- data.frame(type = c("CDS", "CDS", "ncRNA"), locus_tag = c("b2", "b1", "rna"),
  gene = c("two", "one", "rna"), start = c(100, 300, 600), stop = c(200, 499, 700))
stopifnot(identical(sample_cds(annotation, 2, 42), sample_cds(annotation[3:1, ], 2, 42)),
          identical(sort(sample_cds(annotation, 2, 1)$gene_nt), c(101, 200)))
expect_error(sample_cds(annotation, 3, 1), "available CDS")
stopifnot(nrow(sample_cds(rbind(annotation, annotation[1, ]), 2, 1)) == 2L)
split_cds <- rbind(transform(annotation[1, ], type = "gene", start = 50, stop = 250),
                   annotation[1, ], transform(annotation[1, ], start = 220, stop = 250))
stopifnot(sample_cds(split_cds, 1, 1)$gene_nt == 201)
stopifnot(read_benchmark_result("nonexistent", "b1", 7)$status == "process_error",
          read_benchmark_result("nonexistent", "b1", 0)$status == "missing_result")

local({
tmp <- tempfile("benchmark-test-")
dir.create(tmp)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
set.seed(4)
genome <- paste(sample(c("A", "C", "G", "T"), 1200, replace = TRUE), collapse = "")
genome_path <- file.path(tmp, "genome.fna")
writeLines(c(">chr", genome), genome_path)
seqs <- c(substr(genome, 251, 270), rc(substr(genome, 381, 400)),
          substr(genome, 601, 620), rc(substr(genome, 751, 770)))
tail_base <- setdiff(c("A", "C", "G", "T"), substr(genome, 250, 250))[[1]]
tailed <- paste0("ACG", tail_base, seqs[[1]])
stopifnot(annealing_suffix(tailed, genome) == seqs[[1]])
expect_error(annealing_suffix("NNACGT", genome), "A/C/G/T")
expect_error(annealing_suffix(strrep("A", 20), genome), "No exact")
expect_error(annealing_suffix(seqs[[1]], paste0(genome, genome)), "Ambiguous")

purposes <- c("left_arm_forward_primer", "left_arm_reverse_primer", "right_arm_forward_primer",
              "right_arm_reverse_primer", "screening_forward_primer", "screening_reverse_primer")
manifest <- list()
auto <- file.path(tmp, "auto")
for (run in 1:2) {
  attempt <- sprintf("run_%03d_gene_0001", run)
  wet <- file.path(tmp, "raw", attempt, "WetLab/b1_results")
  dir.create(wet, recursive = TRUE)
  metrics <- data.frame(name = paste0("p", 1:6), purpose = purposes,
    annealing_sequence = c(if (run == 1) seqs else rev(seqs), seqs[1:2]), tm_c = 55:60)
  write.table(metrics, file.path(wet, "wet_lab_report.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  d <- read_dhole_primers(file.path(wet, "wet_lab_report.txt"), if (run == 1) "+" else "-")
  stopifnot(identical(d$sequence, seqs), identical(d$reported_tm_c, if (run == 1) as.numeric(55:58) else as.numeric(58:55)))
  for (variant in c("deletion", "bridge")) {
    manifest[[length(manifest) + 1L]] <- data.frame(Name = paste0(attempt, "_", variant),
      variant = variant, attempt = attempt, run = run, target = "b1", contig = "chr",
      strand = if (run == 1) "+" else "-", deletion_start = 401, deletion_end = 600,
      result_dir = paste0("raw/", attempt), target_dir = "b1_results", genome_md5 = unname(tools::md5sum(genome_path)))
  }
}
manifest <- do.call(rbind, manifest)
manifest_path <- file.path(tmp, "manifest.csv")
write.csv(manifest, manifest_path, row.names = FALSE)
for (variant in c("deletion", "bridge")) {
  dir.create(file.path(auto, variant), recursive = TRUE)
  for (arm in c("UHA", "DHA")) {
    m <- manifest[manifest$variant == variant, ]
    if (variant == "bridge" && arm == "DHA") m <- m[1, , drop = FALSE]
    pair <- if (arm == "UHA") c(tailed, seqs[[2]]) else seqs[3:4]
    x <- data.frame(ID = paste0(m$Name, ";chr:400-600"),
      PRIMER_LEFT_WHOLE_SEQUENCE = pair[[1]], PRIMER_RIGHT_WHOLE_SEQUENCE = pair[[2]])
    if (arm == "DHA") names(x) <- c("Region", "primer_f_seq_(5'-3')", "primer_r_seq_(5'-3')")
    write.csv(x, file.path(auto, variant, paste0("Primer_", arm, ".csv")), row.names = FALSE)
  }
}
parsed <- read_auto_arm(file.path(auto, "bridge/Primer_DHA.csv"), "DHA")
stopifnot(nrow(parsed) == 1L, parsed$region == "chr:400-600")
empty <- parsed[FALSE, c("ID", "forward", "reverse")]
names(empty) <- c("ID", "PRIMER_LEFT_WHOLE_SEQUENCE", "PRIMER_RIGHT_WHOLE_SEQUENCE")
write.csv(empty, file.path(tmp, "empty.csv"), row.names = FALSE)
stopifnot(nrow(read_auto_arm(file.path(tmp, "empty.csv"), "UHA")) == 0L)

buffer <- c(monovalent_salt_mm = 50, divalent_salt_mm = 1.5, dntp_mm = 0.6, dna_nm = 50)
oligotm <- Sys.getenv("DHOLE_TEST_OLIGOTM", "primer3/src/oligotm")
if (file.exists(oligotm)) {
  stopifnot(abs(calculate_tm("GCTCTTTTGGTTATGTTGGTATGC", oligotm, buffer) - 58.897687) < 1e-5)
  out <- compare(tmp, manifest_path, auto, file.path(tmp, "comparison"), genome_path, oligotm, buffer)
  stopifnot(nrow(out$primers) == 16L, sum(out$primers$status == "matched") == 14L,
    nrow(out$pairs) == 8L, sum(out$summary$matched_pairs) == 7L,
    all(out$primers$tm_difference_c[out$primers$status == "matched"] == 0),
    all(out$pairs$delta_tm_difference_c[!is.na(out$pairs$delta_tm_difference_c)] == 0))
  toy <- out$primers[1:2, ]
  toy$dhole_tm_c <- c(60, 57)
  toy$auto_tm_c <- c(63, 55)
  paired <- pair_comparison(toy)
  stopifnot(paired$dhole_delta_tm_c == 3, paired$auto_delta_tm_c == 8, paired$delta_tm_difference_c == 5)
  expect_error(compare(tmp, manifest_path, auto, file.path(tmp, "comparison"), genome_path, oligotm, buffer), "already exists")
  bad <- read.csv(file.path(auto, "deletion/Primer_UHA.csv"), check.names = FALSE)
  bad$ID[[1]] <- sub("400-600", "401-600", bad$ID[[1]])
  write.csv(bad, file.path(auto, "deletion/Primer_UHA.csv"), row.names = FALSE)
  expect_error(compare(tmp, manifest_path, auto, file.path(tmp, "bad"), genome_path, oligotm, buffer), "coordinates differ")
  writeLines(c(">chr", paste0("A", genome)), genome_path)
  expect_error(compare(tmp, manifest_path, auto, file.path(tmp, "wrong"), genome_path, oligotm, buffer), "Genome differs")
} else message("SKIP: oligotm integration; set DHOLE_TEST_OLIGOTM or run install.sh")

# Exercise the real benchmark runner and process failures with a synthetic designer.
project <- file.path(tmp, "project with spaces")
dir.create(file.path(project, "benchmark"), recursive = TRUE)
dir.create(file.path(project, "local/data"), recursive = TRUE)
file.copy("benchmark/run_dhole.R", file.path(project, "benchmark/run_dhole.R"))
fixture <- data.frame(type = "CDS", locus_tag = c("b1", "b2", "b3"), gene = c("one", "two", "three"),
  start = 301, stop = 700, seqid = "chr", strand = "+")
write.table(fixture, file.path(project, "local/data/MG1655.gff"), sep = "\t", row.names = FALSE, quote = FALSE)
writeLines(c(">chr", genome), file.path(project, "local/data/MG1655.fna"))
writeLines(c(
  'read_genome_annotation <- function(path, format) read.delim(path)',
  'main <- function(args) {',
  ' target <- args[[match("--cds", args) + 1L]]',
  ' out <- args[[match("--output-dir", args) + 1L]]',
  ' if (target == "b2") quit(status = 7)',
  ' tech <- file.path(out, "TechReport")',
  ' dest <- file.path(tech, paste0(target, "_results"))',
  ' dir.create(dest, recursive = TRUE)',
  ' summary <- data.frame(gene = target, status = if (target == "b3") "error" else "ok", reason = "fixture", output_dir = dest)',
  ' write.table(summary, file.path(tech, "design_summary.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)',
  ' writeLines(c("deletion_start\t401", "deletion_end\t600", "deleted_nt\t200"), file.path(dest, "report.tsv"))',
  '}'), file.path(project, "oligo_designer.R"))
output <- file.path(project, "results")
argv <- c("--vanilla", file.path(project, "benchmark/run_dhole.R"), "2", "3", "42", output)
log <- file.path(tmp, "runner.log")
status <- system2(file.path(R.home("bin"), "Rscript"), vapply(argv, shQuote, character(1)), stdout = log, stderr = log)
if (status != 0) stop(paste(readLines(log), collapse = "\n"))
timing <- read.csv(file.path(output, "timings.csv"))
stopifnot(nrow(timing) == 6L, all(table(timing$status) == 2L), all(timing$elapsed_s >= 0),
  all(timing$gene_nt == 400), all(timing$deleted_nt[timing$status == "ok"] == 200),
  all(is.na(timing$deleted_nt[timing$status != "ok"])), all(timing$exit_code[timing$target == "b2"] == 7L),
  identical(timing$target[1:3], timing$target[4:6]), nrow(read.csv(file.path(output, "run_timings.csv"))) == 2L)
})
cat("Benchmark R tests passed.\n")
