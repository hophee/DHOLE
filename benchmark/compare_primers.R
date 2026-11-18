#!/usr/bin/env Rscript

read_csv_text <- function(path) {
  read.csv(path, check.names = FALSE, colClasses = "character", na.strings = "",
           fileEncoding = "UTF-8-BOM")
}

check_dna <- function(sequence) {
  sequence <- toupper(trimws(sequence))
  if (anyNA(sequence) || any(!grepl("^[ACGT]+$", sequence))) stop("Expected nonempty A/C/G/T sequences")
  sequence
}

rc <- function(sequence) as.character(Biostrings::reverseComplement(Biostrings::DNAString(sequence)))

annealing_suffix <- function(primer, template) {
  primer <- check_dna(primer)
  if (nchar(primer) < 15L) stop("Primer shorter than 15 nt")
  template <- Biostrings::DNAString(template)
  # Exact longest 3' match measures contiguous template binding, including coincidentally matching tail bases.
  for (length in seq.int(nchar(primer), 15L)) {
    suffix <- substring(primer, nchar(primer) - length + 1L)
    hits <- Biostrings::countPattern(suffix, template, fixed = TRUE)
    if (hits > 1L) stop("Ambiguous annealing suffix in genomic flank")
    if (hits == 1L) return(suffix)
  }
  stop("No exact annealing suffix of at least 15 nt in the expected genomic flank")
}

read_auto_arm <- function(path, arm) {
  sheet <- paste0("Primer_", arm)
  location <- paste0(path, " [", sheet, "]")
  x <- tryCatch(
    as.data.frame(readxl::read_xlsx(path, sheet = sheet, col_types = "text",
                                  trim_ws = FALSE, .name_repair = "minimal")),
    error = function(e) stop("Cannot read ", location, ": ", conditionMessage(e), call. = FALSE)
  )
  modern <- c("ID", "PRIMER_LEFT_WHOLE_SEQUENCE", "PRIMER_RIGHT_WHOLE_SEQUENCE")
  original <- c("Region", "primer_f_seq_(5'-3')", "primer_r_seq_(5'-3')")
  columns <- if (all(modern %in% names(x))) modern else original
  if (!all(columns %in% names(x))) stop("Unsupported AutoESDCas columns: ", location)
  x <- x[, columns, drop = FALSE]
  names(x) <- c("ID", "forward", "reverse")
  if (anyNA(x$ID) || any(!grepl("^[^;]+;[^;]+:[0-9]+-[0-9]+$", x$ID)))
    stop("Expected Name;contig:start-end IDs in ", location)
  x$Name <- sub(";.*$", "", x$ID)
  x$region <- sub("^[^;]+;", "", x$ID)
  x$arm <- rep(arm, nrow(x))
  if (anyDuplicated(x$Name)) stop("Multiple designs for one Name in ", location, "; select one explicitly")
  x$forward <- check_dna(x$forward)
  x$reverse <- check_dna(x$reverse)
  x
}

read_dhole_primers <- function(path, strand) {
  lines <- readLines(path, warn = FALSE)
  header <- which(lines == "name\tpurpose\tannealing_sequence\ttm_c")
  if (length(header) != 1L || header + 6L > length(lines)) stop("Missing primer metrics table: ", path)
  x <- read.delim(text = paste(lines[header + 0:6], collapse = "\n"), quote = "", check.names = FALSE)
  purposes <- c("left_arm_forward_primer", "left_arm_reverse_primer",
                "right_arm_forward_primer", "right_arm_reverse_primer")
  if (anyDuplicated(x$purpose) || !all(purposes %in% x$purpose)) stop("Invalid primer purposes: ", path)
  x <- x[match(purposes, x$purpose), , drop = FALSE]
  if (!strand %in% c("+", "-")) stop("Invalid strand")
  # Minus-strand target orientation reverses both arm order and primer direction.
  if (strand == "-") x <- x[4:1, , drop = FALSE]
  data.frame(arm = rep(c("UHA", "DHA"), each = 2), direction = rep(c("F", "R"), 2),
             sequence = check_dna(x$annealing_sequence), reported_tm_c = as.numeric(x$tm_c))
}

calculate_tm <- function(sequences, executable, buffer) {
  sequences <- unique(check_dna(sequences))
  args <- c("-mv", buffer[["monovalent_salt_mm"]], "-dv", buffer[["divalent_salt_mm"]],
            "-n", buffer[["dntp_mm"]], "-d", buffer[["dna_nm"]], "-tp", "1", "-sc", "1")
  values <- vapply(sequences, function(sequence) {
    output <- system2(executable, c(args, sequence), stdout = TRUE, stderr = TRUE)
    value <- suppressWarnings(as.numeric(output))
    if (!is.null(attr(output, "status")) || length(value) != 1L || !is.finite(value))
      stop("oligotm failed for ", sequence, ": ", paste(output, collapse = " "))
    value
  }, numeric(1))
  setNames(values, sequences)
}

pair_comparison <- function(primers) {
  keys <- c("Name", "variant", "attempt", "run", "target", "arm")
  values <- c("dhole_tm_c", "auto_tm_c")
  pairs <- merge(primers[primers$direction == "F", c(keys, values)],
                 primers[primers$direction == "R", c(keys, values)], by = keys, suffixes = c("_F", "_R"))
  pairs$dhole_delta_tm_c <- abs(pairs$dhole_tm_c_F - pairs$dhole_tm_c_R)
  pairs$auto_delta_tm_c <- abs(pairs$auto_tm_c_F - pairs$auto_tm_c_R)
  pairs$delta_tm_difference_c <- pairs$auto_delta_tm_c - pairs$dhole_delta_tm_c
  pairs
}

plot_primer_comparison <- function(primers, pairs, variants, path) {
  outline <- "#252A34"
  fills <- c(UHA = "#0072B2", DHA = "#D55E00")
  shapes <- c(deletion = 21, bridge = 24)
  labels <- c(deletion = "AutoESDCas-deletion", bridge = "AutoESDCas-insertion (bridge)")
  auto_labels <- if (length(variants) == 1L) "AutoESDCas" else labels[variants]

  grDevices::png(path, width = 1400, height = 740, res = 140)
  on.exit(grDevices::dev.off(), add = TRUE)
  layout(matrix(c(1, 2, 3, 3), nrow = 2, byrow = TRUE), heights = c(1, 0.12))
  par(mar = c(4, 4, 2.5, 1), las = 1, pty = "s", mgp = c(2.4, 0.7, 0))

  panel <- function(rows, auto, dhole, title, from_zero = FALSE) {
    matched <- is.finite(rows[[auto]]) & is.finite(rows[[dhole]])
    x <- rows[[auto]][matched]
    y <- rows[[dhole]][matched]
    limits <- if (length(x)) range(c(x, y)) else c(0, 1)
    if (from_zero) limits <- c(0, max(limits))
    if (diff(limits) == 0) limits <- if (from_zero) c(0, 1) else limits + c(-0.5, 0.5)
    plot(x, y, type = "n", xlim = limits, ylim = limits, asp = 1,
         xlab = "AutoESDCas (\u00b0C)", ylab = "DHOLE (\u00b0C)", main = title)
    abline(0, 1, col = "#A7ADB5", lty = 2)
    points(x, y, pch = shapes[rows$variant[matched]], col = outline,
           bg = fills[rows$arm[matched]], cex = 1.05)
    if (!length(x)) text(mean(limits), limits[1] + 0.65 * diff(limits), "No matched results")
  }

  panel(primers, "auto_tm_c", "dhole_tm_c", "Tm per primer")
  panel(pairs, "auto_delta_tm_c", "dhole_delta_tm_c", "Absolute \u0394Tm per pair", from_zero = TRUE)
  par(mar = rep(0, 4), pty = "m")
  plot.new()
  legend("center", c(auto_labels, names(fills)), horiz = TRUE, bty = "n", text.width = NA,
         pch = c(shapes[variants], 22, 22), col = outline,
         pt.bg = c(rep("white", length(variants)), fills), pt.cex = 1.15)
}

publish_comparison <- function(staged, out) {
  lock <- paste0(out, ".lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("Cannot lock output directory: ", out)
  on.exit(unlink(lock, recursive = TRUE, expand = FALSE), add = TRUE)
  if (file.exists(out) && !dir.exists(out)) stop("Output path is a file: ", out)

  archive <- NULL
  if (dir.exists(out)) {
    parent <- dirname(out)
    prefix <- paste0(basename(out), "-")
    entries <- list.files(parent, all.files = TRUE, no.. = TRUE)
    directories <- entries[dir.exists(file.path(parent, entries))]
    suffixes <- substring(directories[startsWith(directories, prefix)], nchar(prefix) + 1L)
    number <- max(0, as.numeric(suffixes[grepl("^[0-9]+$", suffixes)]))
    repeat {
      next_number <- number + 1
      if (!is.finite(next_number) || next_number <= number)
        stop("Cannot increment the output archive suffix for ", out)
      archive <- paste0(out, "-", sprintf("%.0f", next_number))
      if (!file.exists(archive) && !basename(archive) %in% entries) break
      number <- next_number
    }
    if (!file.rename(out, archive)) stop("Cannot move previous results to ", archive)
  }
  if (!file.rename(staged, out)) {
    if (!is.null(archive) && !file.exists(out) && file.rename(archive, out)) archive <- NULL
    stop("Cannot publish results to ", out,
         if (!is.null(archive)) paste0(". Previous results are in ", archive))
  }
  if (!is.null(archive)) message("Previous results: ", archive)
}

compare <- function(results, manifest_path, auto_paths, out, genome_path, oligotm, buffer) {
  if (basename(out) %in% c("", ".", "..")) stop("Specify a named output directory: ", out)
  manifest <- read_csv_text(manifest_path)
  required <- c("Name", "variant", "attempt", "run", "target", "contig", "strand", "deletion_start",
                "deletion_end", "result_dir", "target_dir", "genome_md5")
  if (!all(required %in% names(manifest)) || !nrow(manifest) || anyDuplicated(manifest$Name) ||
      anyNA(manifest[, required]) || !all(manifest$variant %in% c("deletion", "bridge")))
    stop("Invalid export manifest")
  missing_variants <- setdiff(names(auto_paths), manifest$variant)
  if (length(missing_variants))
    stop("Requested variants absent from the manifest: ", paste(missing_variants, collapse = ", "))
  manifest <- manifest[manifest$variant %in% names(auto_paths), , drop = FALSE]
  if (!all(manifest$genome_md5 == unname(tools::md5sum(genome_path))))
    stop("Genome differs from the FASTA used for AutoESDCas export")
  genomes <- Biostrings::readDNAStringSet(genome_path)
  names(genomes) <- sub(" .*", "", names(genomes))
  if (length(genomes) != 1L || !all(manifest$contig == names(genomes))) stop("Genome/manifest contig mismatch")
  genome <- as.character(genomes[[1]])
  auto <- list()
  for (variant in names(auto_paths)) {
    path <- auto_paths[[variant]]
    for (arm in c("UHA", "DHA")) {
      rows <- read_auto_arm(path, arm)
      if (any(!rows$Name %in% manifest$Name[manifest$variant == variant]))
        stop("AutoESDCas names absent from the ", variant, " manifest rows: ", path, " [Primer_", arm, "]")
      rows$variant <- rep(variant, nrow(rows))
      auto[[length(auto) + 1L]] <- rows
    }
  }
  auto <- do.call(rbind, auto)
  all_rows <- list()
  for (i in seq_len(nrow(manifest))) {
    m <- manifest[i, ]
    start <- as.integer(m$deletion_start); end <- as.integer(m$deletion_end)
    if (is.na(start) || is.na(end) || start < 2L || end < start || end >= nchar(genome))
      stop("Invalid deletion bounds for ", m$Name)
    expected_region <- paste0(m$contig, ":", start - 1L, "-", end)
    dhole <- read_dhole_primers(file.path(results, m$result_dir, "WetLab", m$target_dir, "wet_lab_report.txt"), m$strand)
    for (arm in c("UHA", "DHA")) {
      a <- auto[auto$Name == m$Name & auto$arm == arm, , drop = FALSE]
      if (nrow(a) && a$region != expected_region) stop("AutoESDCas deletion coordinates differ: ", m$Name)
      flank <- if (arm == "UHA") substr(genome, 1L, start - 1L) else substring(genome, end + 1L)
      for (direction in c("F", "R")) {
        d <- dhole[dhole$arm == arm & dhole$direction == direction, ]
        full <- sequence <- NA_character_
        status <- "missing_auto_design"
        reason <- NA_character_
        if (nrow(a)) {
          full <- a[[if (direction == "F") "forward" else "reverse"]]
          sequence <- tryCatch(annealing_suffix(full, if (direction == "F") flank else rc(flank)),
            error = function(e) { reason <<- conditionMessage(e); NA_character_ })
          status <- if (is.na(sequence)) "invalid_auto_annealing" else "matched"
        }
        all_rows[[length(all_rows) + 1L]] <- data.frame(
          m[, c("Name", "variant", "attempt", "run", "target")], arm = arm, direction = direction,
          dhole_sequence = d$sequence, dhole_reported_tm_c = d$reported_tm_c,
          auto_full_sequence = full, auto_annealing_sequence = sequence, status = status, reason = reason)
      }
    }
  }
  primers <- do.call(rbind, all_rows)
  sequences <- c(primers$dhole_sequence, primers$auto_annealing_sequence)
  tm <- calculate_tm(sequences[!is.na(sequences)], oligotm, buffer)
  primers$dhole_tm_c <- unname(tm[primers$dhole_sequence])
  primers$auto_tm_c <- unname(tm[primers$auto_annealing_sequence])
  primers$tm_difference_c <- primers$auto_tm_c - primers$dhole_tm_c
  pairs <- pair_comparison(primers)
  summary <- do.call(rbind, lapply(split(pairs, interaction(pairs$variant, pairs$arm, drop = TRUE)), function(x) {
    paired <- is.finite(x$auto_delta_tm_c)
    data.frame(variant = x$variant[[1]], arm = x$arm[[1]], expected_pairs = nrow(x), matched_pairs = sum(paired),
      mean_dhole_delta_tm_c = if (any(paired)) mean(x$dhole_delta_tm_c[paired]) else NA_real_,
      mean_auto_delta_tm_c = if (any(paired)) mean(x$auto_delta_tm_c[paired]) else NA_real_,
      mean_delta_tm_difference_c = if (any(paired)) mean(x$delta_tm_difference_c[paired]) else NA_real_)
  }))
  parent <- dirname(out)
  if (!dir.exists(parent)) dir.create(parent, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(parent)) stop("Cannot create output parent directory: ", parent)
  parent <- normalizePath(parent, winslash = "/", mustWork = TRUE)
  staged <- tempfile(".comparison-", tmpdir = parent)
  if (!dir.create(staged)) stop("Cannot create staging directory: ", staged)
  on.exit(unlink(staged, recursive = TRUE, expand = FALSE), add = TRUE)
  write.csv(primers, file.path(staged, "primer_comparison.csv"), row.names = FALSE, na = "")
  write.csv(pairs, file.path(staged, "pair_comparison.csv"), row.names = FALSE, na = "")
  write.csv(summary, file.path(staged, "summary.csv"), row.names = FALSE, na = "")
  plot_primer_comparison(primers, pairs, names(auto_paths), file.path(staged, "primer_comparison.png"))
  writeLines(c("Both tools: annealing sequences, Primer3 oligotm -tp 1 -sc 1.",
    paste(names(buffer), buffer, sep = "="), paste("oligotm=", normalizePath(oligotm)),
    paste("genome_md5=", tools::md5sum(genome_path)),
    "AutoESDCas: longest exact 3-prime genomic match >=15 nt; no mismatches allowed.",
    "Delta Tm = abs(Tm_F - Tm_R) within each UHA or DHA pair of the same gene and attempt.",
    "Differences = AutoESDCas minus DHOLE. Summary means use matched pairs only.",
    "Repeated attempts are retained; they are not independent biological replicates."), file.path(staged, "methods.txt"))
  publish_comparison(staged, file.path(parent, basename(out)))
  message("Results: ", out, "\nMatched primers: ", sum(primers$status == "matched"), "/", nrow(primers))
  invisible(list(primers = primers, pairs = pairs, summary = summary))
}

parse_comparison_args <- function(args) {
  parser <- argparser::arg_parser(
    "Compare DHOLE with AutoESDCas. Pass --dhole-results and at least one of --deletion or --bridge.",
    name = "bash tools/run-r benchmark/compare_primers.R",
    hide.opts = TRUE
  )
  parser <- argparser::add_argument(parser,
    c("--dhole-results", "--deletion", "--bridge", "--output", "--genome", "--oligotm"),
    help = c("DHOLE benchmark directory containing autoesdcas/manifest.csv",
             "AutoESDCas Design_result.xlsx for deletion",
             "AutoESDCas Design_result.xlsx for bridge",
             "Output directory (default: DHOLE_RESULTS/comparison). Existing results are saved with a numbered suffix.",
             "Reference FASTA (default: local/data/MG1655.fna)",
             "Primer3 executable (default: primer3/src/oligotm)"))
  options <- argparser::parse_args(parser, args)
  if (is.na(options$dhole_results) || !nzchar(options$dhole_results))
    stop("--dhole-results is required; see --help", call. = FALSE)
  if (is.na(options$deletion) && is.na(options$bridge))
    stop("Pass --deletion, --bridge, or both; see --help", call. = FALSE)
  options
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  options <- parse_comparison_args(args)
  if (!requireNamespace("readxl", quietly = TRUE))
    stop("Missing R package readxl in oligo_design; see benchmark/README.md", call. = FALSE)
  results <- options$dhole_results
  auto_paths <- c(deletion = options$deletion, bridge = options$bridge)
  auto_paths <- auto_paths[!is.na(auto_paths)]
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  dhole <- new.env()
  source(file.path(root, "oligo_designer.R"), local = dhole)
  out <- if (is.na(options$output)) file.path(results, "comparison") else options$output
  genome <- if (is.na(options$genome)) file.path(root, "local/data/MG1655.fna") else options$genome
  oligotm <- if (is.na(options$oligotm)) file.path(root, "primer3/src/oligotm") else options$oligotm
  if (!file.exists(oligotm)) stop("Missing Primer3 oligotm: ", oligotm, "; run install.sh")
  compare(results, file.path(results, "autoesdcas/manifest.csv"), auto_paths, out, genome,
          oligotm, dhole$primer3_buffer_parameters())
}

if (sys.nframe() == 0L) main()
