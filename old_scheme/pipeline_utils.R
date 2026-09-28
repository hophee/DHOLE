# Shared command-line and genome-annotation helpers for the ChopChop pipelines.

parse_named_args <- function(args, required, optional = character()) {
  allowed <- c(required, optional)
  values <- setNames(vector("list", length(allowed)), allowed)
  names(values) <- allowed

  i <- 1
  while (i <= length(args)) {
    option <- args[[i]]
    if (!startsWith(option, "--")) {
      stop(sprintf("Ожидался именованный аргумент, получено: %s", option), call. = FALSE)
    }

    name <- gsub("-", "_", sub("^--", "", option))
    if (!name %in% allowed) {
      stop(sprintf("Неизвестный аргумент: %s", option), call. = FALSE)
    }
    if (!is.null(values[[name]])) {
      stop(sprintf("Аргумент %s указан более одного раза", option), call. = FALSE)
    }
    if (i == length(args) || startsWith(args[[i + 1]], "--")) {
      stop(sprintf("Для аргумента %s не задано значение", option), call. = FALSE)
    }

    values[[name]] <- args[[i + 1]]
    i <- i + 2
  }

  missing <- required[vapply(values[required], is.null, logical(1))]
  if (length(missing) > 0) {
    stop(
      sprintf("Не заданы обязательные аргументы: %s", paste0("--", gsub("_", "-", missing), collapse = ", ")),
      call. = FALSE
    )
  }

  values
}

extract_gff_attribute <- function(attributes, attribute_name) {
  pattern <- paste0("(?:^|;)\\s*", attribute_name, "=([^;]+)")
  matches <- regexec(pattern, attributes, perl = TRUE)
  values <- regmatches(attributes, matches)
  result <- vapply(values, function(x) {
    if (length(x) < 2) NA_character_ else utils::URLdecode(trimws(x[[2]]))
  }, character(1))
  result
}

read_genome_annotation <- function(path, format = "bakta") {
  format <- tolower(format)

  if (format == "bakta") {
    lines <- readLines(path)
    hash_lines <- grep("^#", lines)
    skip_lines <- if (length(hash_lines) > 0) max(hash_lines) - 1 else 0
    annotation <- readr::read_tsv(path, skip = skip_lines, show_col_types = FALSE) |>
      janitor::clean_names()
  } else if (format == "gff") {
    gff <- ape::read.gff(path)
    required_columns <- c("type", "start", "end", "strand", "attributes")
    missing_columns <- setdiff(required_columns, names(gff))
    if (length(missing_columns) > 0) {
      stop(
        sprintf("GFF не содержит обязательные колонки: %s", paste(missing_columns, collapse = ", ")),
        call. = FALSE
      )
    }

    annotation <- data.frame(
      seqid = as.character(gff$seqid),
      type = as.character(gff$type),
      start = as.numeric(gff$start),
      stop = as.numeric(gff$end),
      strand = as.character(gff$strand),
      gene = extract_gff_attribute(gff$attributes, "gene"),
      locus_tag = extract_gff_attribute(gff$attributes, "locus_tag"),
      stringsAsFactors = FALSE
    )
  } else {
    stop("Неподдерживаемый формат аннотации. Допустимы: bakta, gff", call. = FALSE)
  }

  required_columns <- c("type", "start", "stop", "strand", "gene", "locus_tag")
  missing_columns <- setdiff(required_columns, names(annotation))
  if (length(missing_columns) > 0) {
    stop(
      sprintf("Аннотация не содержит обязательные колонки: %s", paste(missing_columns, collapse = ", ")),
      call. = FALSE
    )
  }

  if (!"seqid" %in% names(annotation) && "sequence_id" %in% names(annotation)) {
    annotation$seqid <- as.character(annotation$sequence_id)
  }
  annotation
}

find_target_feature <- function(annotation, target_name) {
  locus_matches <- which(!is.na(annotation$locus_tag) & annotation$locus_tag == target_name)
  if (length(locus_matches) > 0) {
    return(locus_matches[[1]])
  }

  gene_matches <- which(!is.na(annotation$gene) & annotation$gene == target_name)
  if (length(gene_matches) > 0) {
    return(gene_matches[[1]])
  }

  stop(
    sprintf("Не найден ген/feature с locus_tag или gene: %s", target_name),
    call. = FALSE
  )
}
