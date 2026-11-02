#!/usr/bin/env Rscript

benchmark_counts <- function(args) {
  values <- c(10L, 10L, 1L)
  for (i in seq_len(min(length(args), 3L))) {
    x <- suppressWarnings(as.numeric(args[[i]]))
    if (!is.finite(x) || x < 1 || x != floor(x) || x > .Machine$integer.max)
      stop("N1, N2 and seed must be positive integers")
    values[[i]] <- as.integer(x)
  }
  values
}

sample_cds <- function(annotation, n, seed) {
  x <- annotation[toupper(trimws(annotation$type)) %in% "CDS", , drop = FALSE]
  x$target <- ifelse(!is.na(x$locus_tag) & nzchar(x$locus_tag), x$locus_tag, x$gene)
  x <- x[!is.na(x$target) & nzchar(x$target), , drop = FALSE]
  targets <- unique(x$target)
  if (any(grepl("[,[:space:]]", targets))) stop("CDS identifiers contain CLI separators")
  # Match DHOLE's locus_tag-first lookup, including the gene record for split/alternative CDS.
  index <- match(targets, annotation$locus_tag)
  index[is.na(index)] <- match(targets[is.na(index)], annotation$gene)
  x <- annotation[index, , drop = FALSE]
  x$target <- targets
  x <- x[order(x$target, method = "radix"), , drop = FALSE]
  if (n > nrow(x)) stop("N2 exceeds the available CDS count: ", nrow(x))
  set.seed(seed, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
  x <- x[sample.int(nrow(x), n), , drop = FALSE]
  x$gene_nt <- x$stop - x$start + 1L
  x
}

read_benchmark_result <- function(directory, target, exit_code) {
  result <- list(status = "process_error", reason = paste("DHOLE exit code", exit_code),
                 deletion_start = NA_integer_, deletion_end = NA_integer_, deleted_nt = NA_integer_)
  if (exit_code != 0L) return(result)
  path <- file.path(directory, "TechReport", "design_summary.tsv")
  if (!file.exists(path)) {
    result$status <- "missing_result"
    result$reason <- "design_summary.tsv was not written"
    return(result)
  }
  rows <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  rows <- rows[rows$gene == target, , drop = FALSE]
  if (nrow(rows) != 1L) stop("Expected one summary row for ", target)
  result$status <- rows$status[[1]]
  result$reason <- rows$reason[[1]]
  if (result$status == "ok") {
    report <- read.delim(file.path(directory, "TechReport", basename(rows$output_dir[[1]]), "report.tsv"),
                         header = FALSE, quote = "", comment.char = "", col.names = c("key", "value"))
    for (key in c("deletion_start", "deletion_end", "deleted_nt")) {
      value <- suppressWarnings(as.integer(report$value[report$key == key]))
      if (length(value) != 1L || is.na(value)) stop("Missing/invalid ", key)
      result[[key]] <- value
    }
    if (result$deletion_start < 1L || result$deleted_nt < 1L ||
        result$deletion_end - result$deletion_start + 1L != result$deleted_nt)
      stop("Inconsistent deletion coordinates")
  }
  result
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("-h", "--help"))) {
    cat("Usage: bash tools/run-r benchmark/run_dhole.R [N1=10] [N2=10] [seed=1] [output_dir]\n",
        "Same sampled MG1655 CDS in each repetition; one DHOLE process per CDS.\n",
        "DHOLE_TARGET_PLASMID / DHOLE_CAS_PLASMID optionally select real FASTA inputs.\n")
    return(invisible(NULL))
  }
  if (length(args) > 4L) stop("Expected N1, N2, seed, output_dir; see --help")
  counts <- benchmark_counts(args)
  n1 <- counts[[1]]; n2 <- counts[[2]]; seed <- counts[[3]]
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  genome <- file.path(root, "local/data/MG1655.fna")
  gff <- file.path(root, "local/data/MG1655.gff")
  for (path in c(genome, gff))
    if (!file.exists(path) || file.info(path)$size == 0) stop("Missing input: ", path, "; run install.sh")
  dhole <- new.env()
  source(file.path(root, "oligo_designer.R"), local = dhole)
  sampled <- sample_cds(dhole$read_genome_annotation(gff, "gff"), n2, seed)
  out <- if (length(args) == 4L) args[[4]] else file.path(root, "local/runs/benchmark",
    paste0("results-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid()))
  if (file.exists(out)) stop("Output already exists: ", out)
  plasmid <- Sys.getenv("DHOLE_TARGET_PLASMID")
  cas <- Sys.getenv("DHOLE_CAS_PLASMID")
  for (path in c(plasmid, cas)[nzchar(c(plasmid, cas))])
    if (!file.exists(path)) stop("Missing plasmid FASTA: ", path)
  if (!dir.create(out, recursive = TRUE)) stop("Cannot create output: ", out)
  out <- normalizePath(out)
  if (!nzchar(plasmid)) {
    plasmid <- file.path(out, "synthetic_pTarget.fasta")
    # Same synthetic cassette as test/test_run.sh; use a real pTarget for experimental comparisons.
    writeLines(c(">synthetic_test_target_plasmid", paste0(
      "GCAGGGGACTAGTACGTACGTACGTACGTACGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGAATTCTCTAGAGTCGACCTGCAG",
      strrep("A", 200))), plasmid)
  }
  plasmid <- normalizePath(plasmid)
  if (nzchar(cas)) cas <- normalizePath(cas)
  write.csv(sampled, file.path(out, "sampled_cds.csv"), row.names = FALSE)
  inputs <- c(genome, gff, plasmid, cas[nzchar(cas)])
  writeLines(c(sprintf("N1=%d N2=%d seed=%d filtering_level=2", n1, n2, seed),
    "Timing: one fresh DHOLE process per gene, including startup and index generation.",
    "Sampling: one fixed CDS panel, repeated N1 times; per-process R seed is fixed.",
    paste(names(tools::md5sum(inputs)), tools::md5sum(inputs)),
    capture.output(sessionInfo())), file.path(out, "metadata.txt"))
  worker <- file.path(out, "worker.R")
  writeLines(c("args <- commandArgs(trailingOnly = TRUE)",
    "source(args[[1]])", "set.seed(as.integer(args[[2]]))", "main(args[-c(1, 2)])"), worker)
  timing_file <- file.path(out, "timings.csv")
  for (run in seq_len(n1)) {
    for (i in seq_len(n2)) {
      target <- sampled$target[[i]]
      attempt <- sprintf("run_%03d_gene_%04d", run, i)
      relative_dir <- file.path("raw", attempt)
      directory <- file.path(out, relative_dir)
      dir.create(directory, recursive = TRUE)
      argv <- c("--vanilla", worker, file.path(root, "oligo_designer.R"), as.character(seed),
        "--genome", genome, "--genome-annotation", gff, "--annotation-format", "gff",
        "--target-plasmid", plasmid, if (nzchar(cas)) c("--cas-plasmid", cas),
        "--filtering-level", "2", "--output-dir", directory, "--cds", target)
      writeLines(argv, file.path(directory, "command_args.txt"))
      message(sprintf("[%d/%d, %d/%d] %s", run, n1, i, n2, target))
      started <- proc.time()[["elapsed"]]
      exit_code <- system2(file.path(R.home("bin"), "Rscript"), vapply(argv, shQuote, character(1)),
        stdout = file.path(directory, "console.log"), stderr = file.path(directory, "console.log"))
      elapsed <- proc.time()[["elapsed"]] - started
      result <- tryCatch(read_benchmark_result(directory, target, exit_code), error = function(e)
        list(status = "invalid_result", reason = conditionMessage(e), deletion_start = NA_integer_,
             deletion_end = NA_integer_, deleted_nt = NA_integer_))
      row <- data.frame(attempt = attempt, run = run, target = target, gene = sampled$gene[[i]],
        contig = sampled$seqid[[i]], strand = sampled$strand[[i]], gene_nt = sampled$gene_nt[[i]],
        elapsed_s = elapsed, exit_code = exit_code, as.data.frame(result), result_dir = relative_dir)
      write.table(row, timing_file, sep = ",", row.names = FALSE,
        col.names = !file.exists(timing_file), append = file.exists(timing_file), na = "")
    }
  }
  timing <- read.csv(timing_file)
  totals <- aggregate(elapsed_s ~ run, timing, sum)
  totals$attempts <- n2
  totals$successful <- vapply(totals$run, function(r) sum(timing$status[timing$run == r] == "ok"), integer(1))
  write.csv(totals, file.path(out, "run_timings.csv"), row.names = FALSE)
  message("Results: ", out, "\nSuccessful designs: ", sum(timing$status == "ok"), "/", nrow(timing))
  invisible(out)
}

if (sys.nframe() == 0L) main()
