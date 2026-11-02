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
  x <- read_csv_text(path)
  modern <- c("ID", "PRIMER_LEFT_WHOLE_SEQUENCE", "PRIMER_RIGHT_WHOLE_SEQUENCE")
  original <- c("Region", "primer_f_seq_(5'-3')", "primer_r_seq_(5'-3')")
  columns <- if (all(modern %in% names(x))) modern else original
  if (!all(columns %in% names(x))) stop("Unsupported AutoESDCas columns: ", path)
  x <- x[, columns, drop = FALSE]
  names(x) <- c("ID", "forward", "reverse")
  if (anyNA(x$ID) || any(!grepl("^[^;]+;[^;]+:[0-9]+-[0-9]+$", x$ID)))
    stop("Expected Name;contig:start-end IDs in ", path)
  x$Name <- sub(";.*$", "", x$ID)
  x$region <- sub("^[^;]+;", "", x$ID)
  x$arm <- rep(arm, nrow(x))
  if (anyDuplicated(x$Name)) stop("Multiple designs for one Name in ", path, "; select one explicitly")
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

compare <- function(results, manifest_path, auto_dir, out, genome_path, oligotm, buffer) {
  if (file.exists(out)) stop("Output already exists: ", out)
  manifest <- read_csv_text(manifest_path)
  required <- c("Name", "variant", "attempt", "run", "target", "contig", "strand", "deletion_start",
                "deletion_end", "result_dir", "target_dir", "genome_md5")
  if (!all(required %in% names(manifest)) || !nrow(manifest) || anyDuplicated(manifest$Name) ||
      anyNA(manifest[, required]) || !all(manifest$variant %in% c("deletion", "bridge")))
    stop("Invalid export manifest")
  if (!all(manifest$genome_md5 == unname(tools::md5sum(genome_path))))
    stop("Genome differs from the FASTA used for AutoESDCas export")
  genomes <- Biostrings::readDNAStringSet(genome_path)
  names(genomes) <- sub(" .*", "", names(genomes))
  if (length(genomes) != 1L || !all(manifest$contig == names(genomes))) stop("Genome/manifest contig mismatch")
  genome <- as.character(genomes[[1]])
  auto <- list()
  for (variant in c("deletion", "bridge")) {
    for (arm in c("UHA", "DHA")) {
      path <- file.path(auto_dir, variant, paste0("Primer_", arm, ".csv"))
      rows <- read_auto_arm(path, arm)
      if (any(!rows$Name %in% manifest$Name[manifest$variant == variant]))
        stop("AutoESDCas names absent from the matching manifest variant: ", path)
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
  if (!dir.create(out, recursive = TRUE)) stop("Cannot create output: ", out)
  write.csv(primers, file.path(out, "primer_comparison.csv"), row.names = FALSE, na = "")
  write.csv(pairs, file.path(out, "pair_comparison.csv"), row.names = FALSE, na = "")
  write.csv(summary, file.path(out, "summary.csv"), row.names = FALSE, na = "")
  writeLines(c("Both tools: annealing sequences, Primer3 oligotm -tp 1 -sc 1.",
    paste(names(buffer), buffer, sep = "="), paste("oligotm=", normalizePath(oligotm)),
    paste("genome_md5=", tools::md5sum(genome_path)),
    "AutoESDCas: longest exact 3-prime genomic match >=15 nt; no mismatches allowed.",
    "Delta Tm = abs(Tm_F - Tm_R) within each UHA or DHA pair of the same gene and attempt.",
    "Differences = AutoESDCas minus DHOLE. Summary means use matched pairs only.",
    "Repeated attempts are retained; they are not independent biological replicates."), file.path(out, "methods.txt"))
  message("Results: ", out, "\nMatched primers: ", sum(primers$status == "matched"), "/", nrow(primers))
  invisible(list(primers = primers, pairs = pairs, summary = summary))
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  if (any(args %in% c("-h", "--help"))) {
    cat("Usage: bash tools/run-r benchmark/compare_primers.R RESULTS AUTO_RESULTS [OUTPUT] [GENOME] [OLIGOTM]\n",
        "AUTO_RESULTS/{deletion,bridge}/{Primer_UHA,Primer_DHA}.csv: exported AutoESDCas sheets.\n",
        "Manifest: RESULTS/autoesdcas/manifest.csv. See benchmark/README.md.\n")
    return(invisible(NULL))
  }
  if (length(args) < 2L || length(args) > 5L) stop("Expected RESULTS AUTO_RESULTS; see --help")
  script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[[1]])
  root <- normalizePath(file.path(dirname(script), ".."), mustWork = TRUE)
  dhole <- new.env()
  source(file.path(root, "oligo_designer.R"), local = dhole)
  out <- if (length(args) >= 3L) args[[3]] else file.path(args[[1]], "comparison")
  genome <- if (length(args) >= 4L) args[[4]] else file.path(root, "local/data/MG1655.fna")
  oligotm <- if (length(args) >= 5L) args[[5]] else file.path(root, "primer3/src/oligotm")
  if (!file.exists(oligotm)) stop("Missing Primer3 oligotm: ", oligotm, "; run install.sh")
  compare(args[[1]], file.path(args[[1]], "autoesdcas/manifest.csv"), args[[2]], out, genome,
          oligotm, dhole$primer3_buffer_parameters())
}

if (sys.nframe() == 0L) main()
