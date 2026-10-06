#!/usr/bin/env Rscript

read_tsv <- function(path) {
  if (!file.exists(path)) return(data.frame())
  read.delim(path, quote = "", comment.char = "", check.names = FALSE,
             stringsAsFactors = FALSE)
}

write_tsv <- function(x, path) {
  write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

cds_ids <- function(path) {
  lines <- readLines(path, warn = FALSE)
  fasta <- match("##FASTA", lines)
  if (!is.na(fasta)) lines <- head(lines, fasta - 1L)
  lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
  fields <- strsplit(lines, "\t", fixed = TRUE)
  ids <- vapply(fields, function(row) {
    if (length(row) != 9L || row[[3L]] != "CDS") return(NA_character_)
    attrs <- strsplit(row[[9L]], ";", fixed = TRUE)[[1L]]
    for (key in c("locus_tag=", "gene=")) {
      value <- attrs[startsWith(attrs, key)]
      if (length(value)) return(utils::URLdecode(substring(value[[1L]], nchar(key) + 1L)))
    }
    NA_character_
  }, character(1))
  unique(ids[!is.na(ids) & nzchar(ids)])
}

collect_results <- function(raw, sampled) {
  targets <- sampled
  targets$status <- "missing_result"
  targets$stage <- "run_failure_or_incomplete"
  targets$reason <- NA_character_
  primers <- rankings <- list()
  for (run in unique(sampled$run)) {
    tech <- file.path(raw, sprintf("run_%03d", run), "TechReport")
    summary <- read_tsv(file.path(tech, "design_summary.tsv"))
    for (i in which(targets$run == run)) {
      gene <- targets$gene[[i]]
      hit <- if (nrow(summary)) match(gene, summary$gene) else NA_integer_
      if (!is.na(hit)) {
        targets[i, c("status", "stage", "reason")] <-
          summary[hit, c("status", "stage", "reason")]
      }
      # DHOLE records the actual directory (which can use a resolved gene name).
      dir <- if (!is.na(hit)) summary$output_dir[[hit]] else NA_character_
      if (is.na(dir) || !dir.exists(dir)) next
      ranking <- read_tsv(file.path(dir, "primer_pair_ranking.tsv"))
      if (nrow(ranking)) {
        ranking$run <- run
        ranking$gene <- gene
        rankings[[length(rankings) + 1L]] <- ranking
      }
      if (targets$status[[i]] != "ok") next
      wet <- summary$wet_lab_dir[[hit]]
      report <- readLines(file.path(wet, "wet_lab_report.txt"), warn = FALSE)
      header <- which(startsWith(report, "name\tpurpose\tannealing_sequence\ttm_c"))
      if (length(header) != 1L) stop("Primer metrics table missing: ", wet)
      # Six final PCR primers: four homology primers and two screening primers.
      metrics <- read.delim(text = paste(report[header + 0:6], collapse = "\n"),
                            quote = "", check.names = FALSE)
      seqs <- toupper(metrics$annealing_sequence)
      metrics$gc_percent <- 100 * nchar(gsub("[^GC]", "", seqs)) / nchar(seqs)
      metrics$run <- run
      metrics$gene <- gene
      primers[[length(primers) + 1L]] <- metrics
    }
  }
  # Ranking schemas can differ between early failures and completed targets.
  bind <- function(xs) {
    if (!length(xs)) return(data.frame())
    columns <- unique(unlist(lapply(xs, names)))
    do.call(rbind, lapply(xs, function(x) {
      for (name in setdiff(columns, names(x))) x[[name]] <- NA
      x[columns]
    }))
  }
  list(targets = targets, primers = bind(primers), rankings = bind(rankings))
}

analyse <- function(out, sampled, timings) {
  library(ggplot2)
  result <- collect_results(file.path(out, "raw"), sampled)
  for (name in names(result)) write_tsv(result[[name]], file.path(out, paste0(name, ".tsv")))
  times <- rbind(data.frame(run = timings$run, metric = "Real", seconds = timings$real_s),
                 data.frame(run = timings$run, metric = "CPU (user + system)", seconds = timings$cpu_s))
  means <- aggregate(seconds ~ metric, times, mean)
  write_tsv(means, file.path(out, "mean_times.tsv"))
  pdf(file.path(out, "analysis.pdf"), width = 10, height = 7)
  on.exit(dev.off(), add = TRUE)
  theme_set(theme_bw(base_size = 12))
  note <- sprintf("%d runs x %d CDS; level 2; all attempts included; nonzero exits: %d",
                  nrow(timings), nrow(sampled) / nrow(timings), sum(timings$exit_code != 0))
  print(ggplot(times, aes(metric, seconds)) + geom_boxplot(outlier.shape = NA) +
          geom_point(position = position_jitter(width = 0.08, seed = 1), alpha = 0.7) +
          geom_point(data = means, colour = "red", size = 3) +
          geom_text(data = means, aes(label = sprintf("Mean: %.2f s", seconds)),
                    vjust = -1, colour = "red") +
          scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
          labs(title = "DHOLE runtime per invocation", subtitle = note, x = NULL, y = "Seconds"))
  counts <- as.data.frame(table(factor(ifelse(result$targets$status == "ok", "Success", "Failure"),
                                       levels = c("Success", "Failure"))))
  names(counts) <- c("category", "count")
  counts$percent <- 100 * counts$count / nrow(sampled)
  write_tsv(counts, file.path(out, "success_rates.tsv"))
  print(ggplot(counts, aes(category, percent, fill = category)) + geom_col() +
          geom_text(aes(label = sprintf("%.1f%% (n=%d)", percent, count)), vjust = -0.4) +
          ylim(0, 110) + guides(fill = "none") +
          labs(title = "Complete primer designs", subtitle = note, x = NULL, y = "% of sampled CDS attempts"))
  empty_plot <- function(title) {
    ggplot() + annotate("text", x = 0, y = 0, label = "No observations available") +
      theme_void() + labs(title = title)
  }
  for (metric in c("tm_c", "gc_percent")) {
    title <- if (metric == "tm_c") "Primer Tm (C)" else "Primer GC (%)"
    p <- result$primers
    if (nrow(p)) {
      p$value <- as.numeric(p[[metric]])
      p <- p[is.finite(p$value), ]
    }
    print(if (!nrow(p)) empty_plot(title) else
      ggplot(p, aes(value)) + geom_histogram(bins = 30, fill = "steelblue", colour = "white") +
        labs(title = title, subtitle = "Final PCR primers across all successful attempts; annealing regions only",
             x = title, y = "Primer observations"))
  }
  ranks <- result$rankings
  if (nrow(ranks)) {
    gates <- intersect(c("structure_passed", "specificity_passed", "openprimer_passed",
                         "strict_qc_passed", "blocking_passed", "selected", "fallback_selected"), names(ranks))
    gate_rows <- do.call(rbind, lapply(gates, function(gate) {
      data.frame(gate = gate, category = ifelse(is.na(ranks[[gate]]), "Not evaluated",
                                               ifelse(ranks[[gate]], "Pass / yes", "Fail / no")))
    }))
    gate_counts <- as.data.frame(table(gate_rows))
    write_tsv(gate_counts, file.path(out, "filter_counts.tsv"))
    print(ggplot(gate_counts, aes(gate, Freq, fill = category)) + geom_col(position = "dodge") +
            coord_flip() + labs(title = "Candidate filtering", subtitle = "All evaluated candidate rows, including failed targets",
                                x = NULL, y = "Candidate rows", fill = NULL))
    reasons <- trimws(unlist(strsplit(ranks$rejection_reason[!is.na(ranks$rejection_reason)], ";", fixed = TRUE)))
    reasons <- reasons[nzchar(reasons)]
    if (length(reasons)) {
      reason_counts <- as.data.frame(table(reasons))
      names(reason_counts) <- c("reason", "count")
      write_tsv(reason_counts, file.path(out, "rejection_counts.tsv"))
      print(ggplot(reason_counts, aes(reorder(reason, count), count)) + geom_col(fill = "tomato") +
              coord_flip() + labs(title = "Candidate rejection categories", x = NULL, y = "Occurrences (categories may overlap)"))
    }
  } else print(empty_plot("Candidate filtering"))
  failed <- result$targets[result$targets$status != "ok", ]
  if (nrow(failed)) {
    failed$stage[is.na(failed$stage)] <- "unknown"
    print(ggplot(failed, aes(stage)) + geom_bar(fill = "tomato") + coord_flip() +
            labs(title = "Failed targets by pipeline stage", x = NULL, y = "CDS attempts"))
  }
  invisible(result)
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) && args[[1L]] %in% c("-h", "--help")) {
    cat("Usage: bash benchmark/run.sh [N=10] [M=10] [seed=1]\n",
        "Requires test/MG1655.fna and test/MG1655.gff; writes benchmark/results-*/\n")
    return(invisible(NULL))
  }
  if (length(args) > 3L) stop("Expected at most N, M, seed; use --help")
  values <- c(10L, 10L, 1L)
  for (i in seq_along(args)) {
    x <- suppressWarnings(as.numeric(args[[i]]))
    if (!is.finite(x) || x < 1 || x > .Machine$integer.max || x != floor(x))
      stop("N, M and seed must be positive integers")
    values[[i]] <- as.integer(x)
  }
  n <- values[[1L]]; m <- values[[2L]]; seed <- values[[3L]]
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1L]])
  root <- dirname(dirname(normalizePath(script)))
  setwd(root)
  genome <- file.path(root, "test", "MG1655.fna")
  annotation <- file.path(root, "test", "MG1655.gff")
  for (path in c(genome, annotation)) {
    if (!file.exists(path) || file.info(path)$size == 0) stop("Missing or empty input: ", path)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 in oligo_design")
  if (!nzchar(Sys.which("bash")) || !nzchar(Sys.which("tar"))) stop("bash and tar are required")
  ids <- cds_ids(annotation)
  if (m > length(ids)) stop("M exceeds available unique CDS: ", length(ids))
  if (any(grepl("[,[:space:]]", ids))) stop("CDS identifiers contain CLI separators")
  set.seed(seed)
  sampled <- do.call(rbind, lapply(seq_len(n), function(run) {
    data.frame(run = run, gene = sample(ids, m, replace = FALSE))
  }))
  out <- file.path(root, "benchmark", paste0("results-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", Sys.getpid()))
  if (!dir.create(out)) stop("Cannot create output directory: ", out)
  raw <- file.path(out, "raw")
  dir.create(raw)
  write_tsv(sampled, file.path(out, "sampled_cds.tsv"))
  writeLines(c(sprintf("N=%d M=%d seed=%d filtering_level=2", n, m, seed),
               paste(names(tools::md5sum(c(genome, annotation))), tools::md5sum(c(genome, annotation))),
               capture.output(sessionInfo())), file.path(out, "metadata.txt"))
  # Same synthetic specificity references as test/test_run.sh.
  plasmid <- file.path(raw, "test_target_plasmid.fasta")
  cas <- file.path(raw, "test_cas_plasmid.fasta")
  writeLines(c(">synthetic_test_target_plasmid", paste0(
    "GCAGGGGACTAGTACGTACGTACGTACGTACGTGTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGCTTTTTTTGAATTCTCTAGAGTCGACCTGCAG",
    strrep("A", 200))), plasmid)
  writeLines(c(">synthetic_test_cas_plasmid",
    "TTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGTTGCAAGCTTAGGCTAACGT"), cas)
  timings <- data.frame()
  for (run in seq_len(n)) {
    message(sprintf("[%d/%d] Designing %d CDS; output: %s", run, n, m, out))
    run_dir <- file.path(raw, sprintf("run_%03d", run))
    dir.create(run_dir)
    command <- c(file.path(R.home("bin"), "Rscript"), "--vanilla", file.path(root, "oligo_designer.R"),
                 "--genome", genome, "--genome-annotation", annotation, "--annotation-format", "gff",
                 "--target-plasmid", plasmid, "--cas-plasmid", cas, "--filtering-level", "2",
                 "--output-dir", run_dir, "--cds", paste(sampled$gene[sampled$run == run], collapse = ","))
    command <- paste(vapply(command, shQuote, character(1), type = "sh"), collapse = " ")
    timing <- file.path(run_dir, "time.tsv")
    log <- file.path(run_dir, "console.log")
    # Bash time includes waited-for child processes (Primer3, CHOPCHOP, etc.).
    shell <- paste0("export LC_ALL=C; TIMEFORMAT=$'%R\\t%U\\t%S'; { time ", command,
                    " >", shQuote(log), " 2>&1; } 2>", shQuote(timing))
    writeLines(shell, file.path(run_dir, "command.sh"))
    status <- system2("bash", c(shQuote(file.path(run_dir, "command.sh"))))
    measured <- scan(timing, quiet = TRUE)
    if (length(measured) != 3L || any(!is.finite(measured))) stop("Invalid timing: ", timing)
    timings <- rbind(timings, data.frame(run = run, exit_code = status, real_s = measured[[1L]],
                        user_s = measured[[2L]], system_s = measured[[3L]], cpu_s = sum(measured[2:3])))
    write_tsv(timings, file.path(out, "timings.tsv"))
  }
  analyse(out, sampled, timings)
  archive <- file.path(out, "dhole-results.tar.gz")
  status <- system2("tar", c("-czf", shQuote(archive), "-C", shQuote(out), "raw"))
  if (status != 0L) stop("Archiving failed; raw results retained: ", raw)
  message("Report: ", file.path(out, "analysis.pdf"), "\nArchive: ", archive)
}

if (sys.nframe() == 0L) main()
