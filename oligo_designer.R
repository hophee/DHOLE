#!/usr/bin/env Rscript
# Batch oligo designer. Functions in this file can also be imported by Shiny.

suppressPackageStartupMessages(library(Biostrings))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(argparser))

.dhole_project_dir <- local({
  source_files <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
  script <- if (length(source_files)) tail(source_files, 1L)[[1]] else {
    sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[[1]])
  }
  dirname(normalizePath(script, mustWork = TRUE))
})
source(file.path(.dhole_project_dir, "callPrimer3.R"), local = TRUE)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

normalize_restriction_site <- function(site, argument_name) {
  site <- toupper(trimws(as.character(site)))
  if (length(site) != 1L || is.na(site) || !nzchar(site) ||
    !grepl("^[ACGT]+$", site)) {
    stop(
      sprintf("%s должен быть непустой последовательностью A/C/G/T", argument_name),
      call. = FALSE
    )
  }
  site
}

reverse_complement_string <- function(sequence) {
  as.character(reverseComplement(DNAString(toupper(sequence))))
}

design_bridge <- function(design_class, deleted_nt) {
  substr("ATGACTGCCCGCAAG", 1L,
         if (design_class == "cds") 15L - deleted_nt %% 3L else 15L)
}

homology_primer_pair <- function(primer_row, side, bridge, site2) {
  c(
    paste0(if (side == "left") "AGCGTCAACT" else bridge,
           primer_row$PRIMER_LEFT_SEQUENCE[[1]]),
    paste0(if (side == "left") reverse_complement_string(bridge) else {
      paste0("ACG", reverse_complement_string(site2))
    }, primer_row$PRIMER_RIGHT_SEQUENCE[[1]])
  )
}

circular_match_positions <- function(sequence, motif) {
  sequence <- toupper(as.character(sequence))
  motif <- toupper(as.character(motif))
  sequence_length <- nchar(sequence)
  motif_length <- nchar(motif)
  if (!sequence_length || !motif_length || motif_length > sequence_length) {
    return(integer())
  }
  extended <- paste0(
    sequence,
    if (motif_length > 1L) substr(sequence, 1L, motif_length - 1L) else ""
  )
  hits <- start(matchPattern(DNAString(motif), DNAString(extended), fixed = TRUE))
  as.integer(unique(hits[hits > 0L & hits <= sequence_length]))
}

find_oriented_restriction_pair <- function(plasmid, site1, site2) {
  sequence <- toupper(as.character(plasmid))
  site1 <- normalize_restriction_site(site1, "site1")
  site2 <- normalize_restriction_site(site2, "site2")
  if (site1 == site2) {
    stop("site1 и site2 должны быть разными", call. = FALSE)
  }

  locate_site <- function(site, label) {
    forward <- circular_match_positions(sequence, site)
    reverse <- circular_match_positions(
      sequence,
      reverse_complement_string(site)
    )
    positions <- sort(unique(c(forward, reverse)))
    if (length(positions) != 1L) {
      stop(
        sprintf(
          "%s должен встречаться в кольцевой pTarget ровно один раз; найдено: %d",
          label,
          length(positions)
        ),
        call. = FALSE
      )
    }
    list(
      position = positions[[1]],
      strands = c(
        if (positions[[1]] %in% forward) "+",
        if (positions[[1]] %in% reverse) "-"
      )
    )
  }

  first <- locate_site(site1, "site1")
  second <- locate_site(site2, "site2")
  orientations <- intersect(first$strands, second$strands)
  if (!length(orientations)) {
    stop("site1 и site2 имеют несовместимую ориентацию", call. = FALSE)
  }
  cassette <- NULL
  if (length(orientations) > 1L) {
    cassette <- locate_circular_pcr_template(
      plasmid,
      "GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCT",
      "AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC",
      nchar(sequence)
    )
    orientation <- cassette$strand
  } else {
    orientation <- orientations[[1]]
  }
  oriented <- if (orientation == "+") {
    sequence
  } else {
    reverse_complement_string(sequence)
  }
  site1_start <- circular_match_positions(oriented, site1)
  site2_start <- circular_match_positions(oriented, site2)
  if (length(site1_start) != 1L || length(site2_start) != 1L) {
    stop("Не удалось однозначно ориентировать пару сайтов", call. = FALSE)
  }

  rotated <- paste0(
    substr(oriented, site1_start, nchar(oriented)),
    if (site1_start > 1L) substr(oriented, 1L, site1_start - 1L) else ""
  )
  site2_rotated <- ((site2_start - site1_start) %% nchar(oriented)) + 1L
  replaced_length <- site2_rotated + nchar(site2) - 1L
  if (!is.null(cassette)) {
    cassette_start <- circular_match_positions(rotated, as.character(cassette$sequence))
    if (length(cassette_start) != 1L || cassette_start <= nchar(site1) ||
        cassette_start + length(cassette$sequence) - 1L >= site2_rotated) {
      stop("sgRNA-кассета не лежит между site1 и site2", call. = FALSE)
    }
  }
  if (
    site2_rotated <= nchar(site1) ||
      replaced_length >= nchar(oriented)
  ) {
    stop(
      "Сайты рестрикции перекрываются или не оставляют pTarget backbone",
      call. = FALSE
    )
  }
  list(
    backbone = substr(rotated, replaced_length + 1L, nchar(rotated)),
    orientation = orientation,
    site1_start = first$position,
    site2_start = second$position
  )
}

model_edited_ptargets <- function(
  plasmid,
  sgrna_products,
  left_arm_product,
  bridge,
  right_arm_product,
  site1 = "ACTAGT",
  site2 = "CTGCAG",
  name_prefix = "pTarget"
) {
  site1 <- normalize_restriction_site(site1, "site1")
  site2 <- normalize_restriction_site(site2, "site2")
  sgrna_products <- toupper(as.character(sgrna_products))
  left_arm_product <- toupper(as.character(left_arm_product))
  right_arm_product <- toupper(as.character(right_arm_product))
  left_overlap <- "AGCGTCAACT"
  if (
    !length(sgrna_products) ||
      length(left_arm_product) != 1L ||
      length(right_arm_product) != 1L ||
      any(!endsWith(sgrna_products, left_overlap)) ||
      !startsWith(left_arm_product, left_overlap) ||
      !endsWith(left_arm_product, bridge) ||
      !startsWith(right_arm_product, bridge)
  ) {
    stop("PCR-продукты не содержат ожидаемые перекрытия сборки", call. = FALSE)
  }
  pair <- find_oriented_restriction_pair(plasmid, site1, site2)
  assembled <- paste0(
    sgrna_products,
    substr(left_arm_product, nchar(left_overlap) + 1L, nchar(left_arm_product)),
    substr(right_arm_product, nchar(bridge) + 1L, nchar(right_arm_product))
  )
  inserts <- substr(assembled, 4L, nchar(assembled) - 3L)
  if (
    any(!startsWith(inserts, site1)) ||
      any(!endsWith(inserts, site2))
  ) {
    stop("PCR-сборка не ограничена ожидаемыми site1/site2", call. = FALSE)
  }
  edited <- DNAStringSet(paste0(inserts, pair$backbone))
  names(edited) <- paste0(name_prefix, "_pTarget_N20_", seq_along(edited))
  list(sequences = edited, restriction_pair = pair)
}

locate_circular_pcr_template <- function(
  plasmid,
  forward_primer,
  reverse_primer,
  max_product_size
) {
  sequence <- toupper(as.character(plasmid))
  forward_primer <- toupper(as.character(forward_primer))
  reverse_binding <- reverse_complement_string(reverse_primer)
  sequence_length <- nchar(sequence)
  candidates <- list()
  for (strand in c("+", "-")) {
    oriented <- if (strand == "+") sequence else {
      reverse_complement_string(sequence)
    }
    forward_starts <- circular_match_positions(oriented, forward_primer)
    reverse_starts <- circular_match_positions(oriented, reverse_binding)
    for (forward_start in forward_starts) {
      for (reverse_start in reverse_starts) {
        offset <- (reverse_start - forward_start) %% sequence_length
        product_length <- offset + nchar(reverse_binding)
        if (
          offset < nchar(forward_primer) ||
            product_length > max_product_size ||
            product_length > sequence_length
        ) {
          next
        }
        oriented_end <- ((forward_start + product_length - 2L) %%
          sequence_length) + 1L
        original_start <- if (strand == "+") {
          forward_start
        } else {
          sequence_length - forward_start + 1L
        }
        original_end <- if (strand == "+") {
          oriented_end
        } else {
          sequence_length - oriented_end + 1L
        }
        candidates[[length(candidates) + 1L]] <- list(
          sequence = DNAString(substr(
            paste0(oriented, oriented),
            forward_start,
            forward_start + product_length - 1L
          )),
          start = original_start,
          end = original_end,
          strand = strand,
          wraps_origin = oriented_end < forward_start
        )
      }
    }
  }
  if (length(candidates) != 1L) {
    stop(
      sprintf(
        "Ожидался один sgRNA PCR-продукт на кольцевой pTarget; найдено: %d",
        length(candidates)
      ),
      call. = FALSE
    )
  }
  candidates[[1]]
}

simulate_full_primer_pcr <- function(
  reaction,
  description,
  location,
  template_amplicon,
  forward_name,
  reverse_name,
  full_forward,
  full_reverse,
  annealing_forward,
  annealing_reverse,
  annealing_temp_c,
  buffer = primer3_buffer_parameters(),
  max_product_size = 2000L
) {
  if (!requireNamespace("DECIPHER", quietly = TRUE)) {
    stop("R-пакет DECIPHER не установлен", call. = FALSE)
  }
  full_forward <- toupper(as.character(full_forward))
  full_reverse <- toupper(as.character(full_reverse))
  annealing_forward <- toupper(as.character(annealing_forward))
  annealing_reverse <- toupper(as.character(annealing_reverse))
  template_amplicon <- toupper(as.character(template_amplicon))
  if (
    !endsWith(full_forward, annealing_forward) ||
      !endsWith(full_reverse, annealing_reverse) ||
      nchar(template_amplicon) <
        nchar(annealing_forward) + nchar(annealing_reverse) ||
      !startsWith(template_amplicon, annealing_forward) ||
      !endsWith(
        template_amplicon,
        reverse_complement_string(annealing_reverse)
      )
  ) {
    stop(
      sprintf("Некорректные праймеры или шаблон для PCR %s", reaction),
      call. = FALSE
    )
  }

  interior_start <- nchar(annealing_forward) + 1L
  interior_end <- nchar(template_amplicon) - nchar(annealing_reverse)
  interior <- if (interior_start <= interior_end) {
    substr(template_amplicon, interior_start, interior_end)
  } else {
    ""
  }
  # AmplifyDNA cannot match long non-annealing 5' tails on the initial template.
  # After cycle one those tails are incorporated, which is the template modelled.
  primed_template <- paste0(
    full_forward,
    interior,
    reverse_complement_string(full_reverse)
  )
  primers <- DNAStringSet(c(forward = full_forward, reverse = full_reverse))
  templates <- DNAStringSet(c(post_first_cycle = primed_template))
  effective_max <- max(as.integer(max_product_size), nchar(primed_template))
  configure_openprimer_environment()
  products <- DECIPHER::AmplifyDNA(
    primers,
    templates,
    maxProductSize = effective_max,
    annealingTemp = annealing_temp_c,
    P = buffer[["dna_nm"]] * 1e-9,
    ions = 0.2,
    includePrimers = TRUE,
    minEfficiency = 0.001,
    processors = 1L
  )
  expected <- which(as.character(products) == primed_template)
  if (length(products) != 1L || length(expected) != 1L) {
    stop(
      sprintf("DECIPHER не вернул однозначный PCR-продукт для %s", reaction),
      call. = FALSE
    )
  }
  product <- products[expected]
  conditions <- paste0(
    "annealing=", format(annealing_temp_c, trim = TRUE), " C; ",
    "primer=", buffer[["dna_nm"]], " nM each; ",
    "monovalent=", buffer[["monovalent_salt_mm"]], " mM; ",
    "Mg=", buffer[["divalent_salt_mm"]], " mM; ",
    "dNTP=", buffer[["dntp_mm"]], " mM; ",
    "DECIPHER ions=0.2 M"
  )
  data.frame(
    name = reaction,
    description = description,
    location = location,
    length_bp = length(product[[1]]),
    primers = paste(forward_name, reverse_name, sep = " + "),
    pcr_conditions = conditions,
    sequence = as.character(product[[1]]),
    stringsAsFactors = FALSE
  )
}

model_design_pcr_products <- function(
  input,
  feature,
  pair,
  sgrnas,
  sgrna_reverse,
  arm_primers,
  screening,
  screening_primers,
  edited_genome,
  screening_product_sizes
) {
  max_product_size <- input$parameters$primer_qc$max_product_size
  buffer <- input$parameters$primer3_buffer
  sgrna_annealing <- substr(
    as.character(sgrnas[[1]]),
    3L + nchar(input$parameters$site1) + 20L + 1L,
    length(sgrnas[[1]])
  )
  ptarget_template <- locate_circular_pcr_template(
    input$target_plasmid_sequence,
    sgrna_annealing,
    as.character(sgrna_reverse[[1]]),
    max_product_size
  )
  ptarget_location <- paste0(
    input$target_plasmid_name %||% "pTarget",
    ":",
    ptarget_template$start,
    "->",
    ptarget_template$end,
    " (",
    ptarget_template$strand,
    ", circular",
    if (ptarget_template$wraps_origin) ", crosses FASTA origin" else "",
    ")"
  )
  rows <- lapply(seq_along(sgrnas), function(i) {
    simulate_full_primer_pcr(
      paste0("sgRNA_N20_", i),
      paste0("ПЦР sgRNA-кассеты для N20_", i),
      ptarget_location,
      ptarget_template$sequence,
      names(sgrnas)[[i]],
      names(sgrna_reverse)[[1]],
      sgrnas[[i]],
      sgrna_reverse[[1]],
      sgrna_annealing,
      sgrna_reverse[[1]],
      60,
      buffer,
      max_product_size
    )
  })

  genome_name <- input$genome_contig %||% "genome"
  for (i in seq_len(2L)) {
    template <- input$genome[
      pair$genome_start[[i]]:pair$genome_end[[i]]
    ]
    if (feature$strand == "-") {
      template <- reverseComplement(template)
    }
    reaction <- c("left_homology_arm", "right_homology_arm")[[i]]
    rows[[length(rows) + 1L]] <- simulate_full_primer_pcr(
      reaction,
      paste(
        "ПЦР",
        c("левого", "правого")[[i]],
        "плеча гомологии"
      ),
      sprintf(
        "%s:%d-%d (%s)",
        genome_name,
        pair$genome_start[[i]],
        pair$genome_end[[i]],
        feature$strand
      ),
      template,
      names(arm_primers)[[2L * i - 1L]],
      names(arm_primers)[[2L * i]],
      arm_primers[[2L * i - 1L]],
      arm_primers[[2L * i]],
      pair$PRIMER_LEFT_SEQUENCE[[i]],
      pair$PRIMER_RIGHT_SEQUENCE[[i]],
      round(min(
        pair$PRIMER_LEFT_TM[[i]],
        pair$PRIMER_RIGHT_TM[[i]]
      ) - 3, 1),
      buffer,
      max_product_size
    )
  }

  screening_start <- screening$genome_start[[1]]
  original_end <- screening$genome_end[[1]]
  edited_end <- screening_start +
    screening_product_sizes[["successful_insertion_bp"]] - 1L
  if (edited_end > length(edited_genome[[1]])) {
    stop("Скрининговый PCR-продукт выходит за границу edited genome", call. = FALSE)
  }
  screening_temp <- round(min(
    screening$PRIMER_LEFT_TM[[1]],
    screening$PRIMER_RIGHT_TM[[1]]
  ) - 3, 1)
  screening_templates <- list(
    original_genome = input$genome[screening_start:original_end],
    edited_genome = edited_genome[[1]][screening_start:edited_end]
  )
  screening_ends <- c(original_end, edited_end)
  screening_descriptions <- c(
    "Скрининговая ПЦР исходного генома",
    "Скрининговая ПЦР редактированного генома"
  )
  for (i in seq_along(screening_templates)) {
    template_name <- names(screening_templates)[[i]]
    rows[[length(rows) + 1L]] <- simulate_full_primer_pcr(
      paste0("screening_", template_name),
      screening_descriptions[[i]],
      sprintf(
        "%s:%d-%d (+)",
        if (i == 1L) genome_name else "edited_genome",
        screening_start,
        screening_ends[[i]]
      ),
      screening_templates[[i]],
      names(screening_primers)[[1]],
      names(screening_primers)[[2]],
      screening_primers[[1]],
      screening_primers[[2]],
      screening$PRIMER_LEFT_SEQUENCE[[1]],
      screening$PRIMER_RIGHT_SEQUENCE[[1]],
      screening_temp,
      buffer,
      max_product_size
    )
  }
  bind_rows(rows)
}

extract_gff_attribute <- function(attributes, attribute_name) {
  pattern <- paste0("(?:^|;)\\s*", attribute_name, "=([^;]+)")
  matches <- regexec(pattern, attributes, perl = TRUE)
  values <- regmatches(attributes, matches)
  vapply(
    values,
    function(x) {
      if (length(x) < 2) NA_character_ else utils::URLdecode(trimws(x[[2]]))
    },
    character(1)
  )
}

read_genome_annotation <- function(path, format = "bakta") {
  format <- tolower(format)
  if (format == "bakta") {
    lines <- readLines(path)
    hash_lines <- grep("^#", lines)
    skip_lines <- if (length(hash_lines)) max(hash_lines) - 1 else 0
    annotation <- read_tsv(path, skip = skip_lines, show_col_types = FALSE) |>
      janitor::clean_names()
  } else if (format == "gff") {
    gff <- ape::read.gff(path)
    required <- c("seqid", "type", "start", "end", "strand", "attributes")
    missing <- setdiff(required, names(gff))
    if (length(missing)) {
      stop(
        sprintf(
          "GFF не содержит обязательные колонки: %s",
          paste(missing, collapse = ", ")
        ),
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
    stop(
      "Неподдерживаемый формат аннотации. Допустимы: tsv (bakta), gff",
      call. = FALSE
    )
  }
  required <- c("type", "start", "stop", "strand", "gene", "locus_tag")
  missing <- setdiff(required, names(annotation))
  if (length(missing)) {
    stop(
      sprintf(
        "Аннотация не содержит обязательные колонки: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if (!"seqid" %in% names(annotation) && "sequence_id" %in% names(annotation)) {
    annotation$seqid <- as.character(annotation$sequence_id)
  }
  if (!"seqid" %in% names(annotation)) {
    stop(
      "Аннотация не содержит колонку seqid или sequence_id",
      call. = FALSE
    )
  }
  annotation
}

find_target_feature <- function(annotation, target_name) {
  by_locus_tag <- which(
    !is.na(annotation$locus_tag) & annotation$locus_tag == target_name
  )
  if (length(by_locus_tag)) {
    return(by_locus_tag[[1]])
  }
  by_gene <- which(!is.na(annotation$gene) & annotation$gene == target_name)
  if (length(by_gene)) {
    return(by_gene[[1]])
  }
  stop(
    sprintf("Не найден ген/feature с locus_tag или gene: %s", target_name),
    call. = FALSE
  )
}

parse_designer_args <- function(args) {
  parser <- arg_parser("Batch oligo designer")
  parser <- add_argument(
    parser,
    "--genome",
    help = "Path to the genome FASTA file"
  )
  parser <- add_argument(
    parser,
    "--genome-annotation",
    help = "Path to the genome annotation file"
  )
  parser <- add_argument(
    parser,
    "--target-plasmid",
    help = "Path to the target plasmid FASTA file"
  )
  # TODO: map enzyme names to recognition sites via a maintained site dictionary.
  parser <- add_argument(
    parser,
    "--site1",
    help = "First restriction-site sequence in insert orientation",
    default = "ACTAGT"
  )
  parser <- add_argument(
    parser,
    "--site2",
    help = "Second restriction-site sequence in insert orientation",
    default = "CTGCAG"
  )
  parser <- add_argument(
    parser,
    "--output-dir",
    help = "Directory for design results"
  )
  parser <- add_argument(
    parser,
    "--cds",
    help = "CDS targets (comma-separated or space-separated)",
    nargs = Inf
  )
  parser <- add_argument(
    parser,
    "--ncrna",
    help = "ncRNA targets (comma-separated or space-separated)",
    nargs = Inf
  )
  parser <- add_argument(
    parser,
    "--cas-plasmid",
    help = "Path to the Cas plasmid FASTA file"
  )
  parser <- add_argument(
    parser,
    "--annotation-format",
    help = "Annotation format: bakta or gff",
    default = "gff"
  )
  parser <- add_argument(
    parser,
    "--chopchop-script",
    help = "Path to chopchop.py",
    default = "chopchop/chopchop.py"
  )
  parser <- add_argument(
    parser,
    "--chopchop-python",
    help = "Python 2 executable for CHOPCHOP",
    default = "chopchop-python"
  )
  parser <- add_argument(
    parser,
    "--primer3",
    help = "Path to primer3_core",
    default = "primer3/src/primer3_core"
  )
  parser <- add_argument(
    parser,
    "--filtering-level",
    help = "Primer QC: 1=lite/report, 2=default/warn, 3=hard/core",
    default = 2L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--n20-mn",
    help = "Minimum number of N20 sequences",
    default = 1L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--n20-strands",
    help = "N20 strand mode: plus, minus, both, or random",
    default = "random"
  )
  parser <- add_argument(
    parser,
    "--n20-offtarget",
    help = "Comma-separated maximum values for MM0, MM1, ...",
    default = "0"
  )
  parser <- add_argument(
    parser,
    "--cds-fs",
    help = "Require the deleted CDS fragment length to be divisible by three",
    flag = TRUE
  )
  parser <- add_argument(
    parser,
    "--ncrna-fs",
    help = "Require the deleted ncRNA fragment length to be divisible by three",
    flag = TRUE
  )
  parser <- add_argument(
    parser,
    "--left-arm-min",
    help = "Minimum left homology arm length",
    default = 300L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--left-arm-opt",
    help = "Preferred initial left homology arm length",
    default = 350L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--left-arm-max",
    help = "Maximum left homology arm length",
    default = 400L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--right-arm-min",
    help = "Minimum right homology arm length",
    default = 400L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--right-arm-opt",
    help = "Preferred initial right homology arm length",
    default = 450L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--right-arm-max",
    help = "Maximum right homology arm length",
    default = 500L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--n20-arm-min-distance",
    help = "Minimum distance from every N20 to both homology arms",
    default = 40L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--primer-max-mismatches",
    help = "Maximum mismatches per primer binding site",
    default = 2L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--primer-critical-3p-bases",
    help = "Length of the critical primer 3-prime region",
    default = 5L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--primer-max-3p-mismatches",
    help = "Maximum mismatches in the critical 3-prime region",
    default = 0L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--primer-min-product-size",
    help = "Minimum counted PCR product size",
    default = 50L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--primer-max-product-size",
    help = "Maximum counted PCR product size",
    default = 2000L,
    type = "integer"
  )
  parser <- add_argument(
    parser,
    "--primer-max-offtarget-products",
    help = "Maximum allowed non-intended PCR products",
    default = 0L,
    type = "integer"
  )

  normalize_scalar <- function(value) {
    if (is.null(value) || length(value) != 1L || is.na(value[[1]])) {
      return(character())
    }
    value <- trimws(as.character(value[[1]]))
    if (nzchar(value)) value else character()
  }
  normalize_targets <- function(value) {
    if (is.null(value) || !length(value)) {
      return(character())
    }
    value <- as.character(value[!is.na(value)])
    targets <- trimws(unlist(
      strsplit(value, ",", fixed = TRUE),
      use.names = FALSE
    ))
    unique(targets[nzchar(targets)])
  }
  parse_offtarget_thresholds <- function(value) {
    value <- normalize_scalar(value)
    if (!length(value)) {
      stop("--n20-offtarget не может быть пустым", call. = FALSE)
    }
    fields <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
    if (!length(fields) || any(!nzchar(fields))) {
      stop(
        "--n20-offtarget должен быть списком целых чисел через запятую",
        call. = FALSE
      )
    }
    thresholds <- suppressWarnings(as.numeric(fields))
    if (
      anyNA(thresholds) ||
        any(thresholds < 0) ||
        any(thresholds != as.integer(thresholds))
    ) {
      stop(
        "--n20-offtarget допускает только неотрицательные целые числа",
        call. = FALSE
      )
    }
    as.integer(thresholds)
  }
  validate_arm_lengths <- function(side, minimum, preferred, maximum) {
    values <- c(minimum, preferred, maximum)
    if (anyNA(values) || any(values < 1L) || !identical(values, sort(values))) {
      stop(
        sprintf(
          "Для %s плеча требуется 0 < min <= opt <= max",
          side
        ),
        call. = FALSE
      )
    }
    as.integer(values)
  }

  parsed <- parse_args(parser, args)
  annotation_format <- tolower(normalize_scalar(parsed$annotation_format))
  if (!annotation_format %in% c("bakta", "gff")) {
    stop(
      "Некорректный --annotation-format. Допустимы: bakta, gff",
      call. = FALSE
    )
  }
  n20_mn <- as.integer(parsed$n20_mn)
  if (length(n20_mn) != 1L || is.na(n20_mn) || n20_mn < 1L) {
    stop("--n20-mn должен быть положительным целым числом", call. = FALSE)
  }
  n20_strands <- tolower(normalize_scalar(parsed$n20_strands))
  strand_aliases <- c(
    "+" = "plus",
    "plus" = "plus",
    "-" = "minus",
    "minus" = "minus",
    "both" = "both",
    "any" = "random",
    "random" = "random"
  )
  if (
    length(n20_strands) != 1L ||
      !n20_strands %in% names(strand_aliases)
  ) {
    stop(
      "Некорректный --n20-strands. Допустимы: plus, minus, both, random",
      call. = FALSE
    )
  }
  n20_strands <- unname(strand_aliases[[n20_strands]])
  filtering_level <- as.integer(parsed$filtering_level)
  if (
    length(filtering_level) != 1L ||
      is.na(filtering_level) ||
      !filtering_level %in% 1:3
  ) {
    stop("--filtering-level должен быть равен 1, 2 или 3", call. = FALSE)
  }
  left_arm <- validate_arm_lengths(
    "левого",
    parsed$left_arm_min,
    parsed$left_arm_opt,
    parsed$left_arm_max
  )
  right_arm <- validate_arm_lengths(
    "правого",
    parsed$right_arm_min,
    parsed$right_arm_opt,
    parsed$right_arm_max
  )
  n20_arm_min_distance <- as.integer(parsed$n20_arm_min_distance)
  if (
    length(n20_arm_min_distance) != 1L ||
      is.na(n20_arm_min_distance) ||
      n20_arm_min_distance < 0L
  ) {
    stop(
      "--n20-arm-min-distance должен быть неотрицательным целым числом",
      call. = FALSE
    )
  }
  qc_values <- c(
    max_mismatches = parsed$primer_max_mismatches,
    critical_3p_bases = parsed$primer_critical_3p_bases,
    max_3p_mismatches = parsed$primer_max_3p_mismatches,
    min_product_size = parsed$primer_min_product_size,
    max_product_size = parsed$primer_max_product_size,
    max_allowed_offtarget_products = parsed$primer_max_offtarget_products
  )
  qc_values <- setNames(as.integer(qc_values), names(qc_values))
  if (
    anyNA(qc_values) ||
      any(qc_values < 0L) ||
      qc_values[["critical_3p_bases"]] < 1L ||
      qc_values[["min_product_size"]] < 1L ||
      qc_values[["min_product_size"]] > qc_values[["max_product_size"]] ||
      qc_values[["max_3p_mismatches"]] > qc_values[["critical_3p_bases"]]
  ) {
    stop("Некорректные параметры primer specificity QC", call. = FALSE)
  }

  values <- list(
    genome = normalize_scalar(parsed$genome),
    genome_annotation = normalize_scalar(parsed$genome_annotation),
    target_plasmid = normalize_scalar(parsed$target_plasmid),
    site1 = normalize_restriction_site(parsed$site1, "--site1"),
    site2 = normalize_restriction_site(parsed$site2, "--site2"),
    output_dir = normalize_scalar(parsed$output_dir),
    cds = normalize_targets(parsed$cds),
    ncrna = normalize_targets(parsed$ncrna),
    cas_plasmid = normalize_scalar(parsed$cas_plasmid),
    annotation_format = annotation_format,
    chopchop_script = normalize_scalar(parsed$chopchop_script),
    chopchop_python = normalize_scalar(parsed$chopchop_python),
    primer3 = normalize_scalar(parsed$primer3),
    filtering_level = filtering_level,
    n20_mn = n20_mn,
    n20_strands = n20_strands,
    n20_offtarget = parse_offtarget_thresholds(parsed$n20_offtarget),
    cds_fs = isTRUE(parsed$cds_fs),
    ncrna_fs = isTRUE(parsed$ncrna_fs),
    left_arm = setNames(left_arm, c("min", "opt", "max")),
    right_arm = setNames(right_arm, c("min", "opt", "max")),
    n20_arm_min_distance = n20_arm_min_distance,
    primer_qc = utils::modifyList(primer_qc_defaults(), as.list(qc_values))
  )
  if (values$site1 == values$site2) {
    stop("--site1 и --site2 должны быть разными", call. = FALSE)
  }

  required <- c("genome", "genome_annotation", "target_plasmid", "output_dir")
  required_options <- c(
    "--genome",
    "--genome-annotation",
    "--target-plasmid",
    "--output-dir"
  )
  missing <- required_options[
    vapply(values[required], length, integer(1)) != 1L
  ]
  if (length(missing)) {
    stop(
      sprintf("Не заданы: %s", paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  if (!length(values$cds) && !length(values$ncrna)) {
    stop("Нужен хотя бы один из аргументов --cds или --ncrna", call. = FALSE)
  }
  values
}

run_tool <- function(command, args, stdout = "", stderr = "") {
  status <- system2(command, args = args, stdout = stdout, stderr = stderr)
  if (status != 0) {
    stop(
      sprintf("External command failed (%s): %s", status, command),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

primer3_buffer_parameters <- function() {
  c(
    monovalent_salt_mm = 50,
    divalent_salt_mm = 1.5,
    dntp_mm = 0.6,
    dna_nm = 50
  )
}

primer_qc_defaults <- function() {
  list(
    max_mismatches = 2L,
    critical_3p_bases = 5L,
    max_3p_mismatches = 0L,
    min_product_size = 50L,
    max_product_size = 2000L,
    max_allowed_offtarget_products = 0L,
    expected_coordinate_tolerance = 0L,
    expected_size_tolerance = 0L,
    hard_max_tm_diff = 5,
    primer_efficiency_min = 0.001,
    openprimer_profile = "C_Taq_PCR_high_stringency.xml",
    openprimer_critical_annealing = c(
      "primer_length",
      "gc_ratio",
      "gc_clamp",
      "no_runs",
      "no_repeats",
      "melting_temp_range",
      "melting_temp_diff"
    ),
    openprimer_critical_full = c(
      "self_dimerization",
      "cross_dimerization",
      "secondary_structure"
    ),
    openprimer_soft_constraints = character()
  )
}

filtering_level_name <- function(level) {
  names <- c("lite", "default", "hard")
  level <- as.integer(level)
  if (length(level) != 1L || is.na(level) || !level %in% seq_along(names)) {
    stop("filtering_level должен быть равен 1, 2 или 3", call. = FALSE)
  }
  names[[level]]
}

make_specificity_references <- function(
  genome,
  target_plasmid,
  cas_plasmid = NULL
) {
  inputs <- list(
    genome = list(path = genome, topology = "linear"),
    target_plasmid = list(path = target_plasmid, topology = "circular"),
    cas_plasmid = list(path = cas_plasmid, topology = "circular")
  )
  records <- list()
  for (reference_type in names(inputs)) {
    input <- inputs[[reference_type]]
    if (is.null(input$path) || !length(input$path) || is.na(input$path) ||
      !nzchar(input$path)) {
      next
    }
    sequences <- readDNAStringSet(input$path)
    if (!length(sequences)) {
      stop(
        sprintf("FASTA не содержит записей: %s", input$path),
        call. = FALSE
      )
    }
    source_names <- names(sequences)
    missing_names <- is.na(source_names) | !nzchar(source_names)
    source_names[missing_names] <- paste0("record_", which(missing_names))
    for (i in seq_along(sequences)) {
      records[[length(records) + 1L]] <- data.frame(
        reference_id = paste(reference_type, source_names[[i]], sep = "::"),
        reference_type = reference_type,
        topology = input$topology,
        source_name = source_names[[i]],
        source_length = length(sequences[[i]]),
        source_path = input$path,
        stringsAsFactors = FALSE
      )
      records[[length(records)]]$sequence <- I(list(sequences[[i]]))
    }
  }
  references <- bind_rows(records)
  if (!nrow(references) || anyDuplicated(references$reference_id)) {
    stop("Не удалось создать уникальный набор specificity references", call. = FALSE)
  }
  references
}

reference_digest <- function(references) {
  digest::digest(list(
    references$reference_id,
    references$reference_type,
    references$topology,
    vapply(references$sequence, as.character, character(1))
  ))
}

extend_reference_for_pcr <- function(reference, max_product_size, primer_length) {
  sequence <- reference$sequence[[1]]
  original_length <- reference$source_length[[1]]
  if (reference$topology[[1]] != "circular") {
    return(sequence)
  }
  extension_length <- min(
    original_length,
    max_product_size + primer_length - 1L
  )
  DNAString(paste0(
    as.character(sequence),
    substr(as.character(sequence), 1L, extension_length)
  ))
}

iupac_mismatch_positions <- function(pattern, observed) {
  pattern <- strsplit(toupper(pattern), "", fixed = TRUE)[[1]]
  observed <- strsplit(toupper(observed), "", fixed = TRUE)[[1]]
  codes <- IUPAC_CODE_MAP
  compatible <- mapply(
    function(x, y) {
      x_bases <- strsplit(unname(codes[[x]]), "", fixed = TRUE)[[1]]
      y_bases <- strsplit(unname(codes[[y]]), "", fixed = TRUE)[[1]]
      length(intersect(x_bases, y_bases)) > 0L
    },
    pattern,
    observed,
    USE.NAMES = FALSE
  )
  which(!compatible)
}

enumerate_primer_binding_sites <- function(
  primer_sequence,
  primer_id,
  references,
  config = primer_qc_defaults()
) {
  primer_sequence <- toupper(as.character(primer_sequence))
  primer_length <- nchar(primer_sequence)
  if (primer_length < 1L) {
    stop("Пустая последовательность праймера", call. = FALSE)
  }
  strand_patterns <- list(
    `+` = primer_sequence,
    `-` = as.character(reverseComplement(DNAString(primer_sequence)))
  )
  sites <- list()
  for (i in seq_len(nrow(references))) {
    reference <- references[i, , drop = FALSE]
    subject <- extend_reference_for_pcr(
      reference,
      config$max_product_size,
      primer_length
    )
    original_length <- reference$source_length[[1]]
    for (strand in names(strand_patterns)) {
      pattern <- strand_patterns[[strand]]
      hits <- matchPattern(
        DNAString(pattern),
        subject,
        algorithm = "auto",
        max.mismatch = config$max_mismatches,
        with.indels = FALSE,
        fixed = FALSE
      )
      if (!length(hits)) {
        next
      }
      for (j in seq_along(hits)) {
        raw_start <- start(hits)[[j]]
        raw_end <- end(hits)[[j]]
        if (reference$topology[[1]] == "circular" &&
          raw_start > original_length + config$max_product_size) {
          next
        }
        observed <- as.character(subseq(subject, raw_start, raw_end))
        pattern_mismatches <- iupac_mismatch_positions(pattern, observed)
        primer_mismatches <- if (strand == "+") {
          pattern_mismatches
        } else {
          primer_length - pattern_mismatches + 1L
        }
        primer_mismatches <- sort(primer_mismatches)
        critical_start <- max(1L, primer_length - config$critical_3p_bases + 1L)
        normalized_start <- if (reference$topology[[1]] == "circular") {
          (raw_start - 1L) %% original_length + 1L
        } else {
          raw_start
        }
        normalized_end <- if (reference$topology[[1]] == "circular") {
          (raw_end - 1L) %% original_length + 1L
        } else {
          raw_end
        }
        sites[[length(sites) + 1L]] <- data.frame(
          site_id = paste(
            primer_id,
            reference$reference_id[[1]],
            strand,
            raw_start,
            raw_end,
            sep = ":"
          ),
          primer_id = primer_id,
          reference_id = reference$reference_id[[1]],
          reference_type = reference$reference_type[[1]],
          start = normalized_start,
          end = normalized_end,
          strand = strand,
          mismatches = length(primer_mismatches),
          mismatch_positions = paste(primer_mismatches, collapse = ","),
          mismatches_3p = sum(primer_mismatches >= critical_start),
          passes_3p = sum(primer_mismatches >= critical_start) <=
            config$max_3p_mismatches,
          circular_wrap = raw_end > original_length,
          duplicated_circular_hit = raw_start > original_length,
          raw_start = raw_start,
          raw_end = raw_end,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(sites)) {
    return(data.frame(
      site_id = character(), primer_id = character(),
      reference_id = character(), reference_type = character(),
      start = integer(), end = integer(), strand = character(),
      mismatches = integer(), mismatch_positions = character(),
      mismatches_3p = integer(), passes_3p = logical(),
      circular_wrap = logical(), duplicated_circular_hit = logical(),
      raw_start = integer(), raw_end = integer(),
      stringsAsFactors = FALSE
    ))
  }
  bind_rows(sites)
}

empty_amplicon_table <- function() {
  data.frame(
    reference_id = character(), reference_type = character(),
    start = integer(), end = integer(), product_size = integer(),
    sequence = character(), circular_wrap = logical(),
    invalid_size = logical(), forward_site_id = character(),
    reverse_site_id = character(), forward_mismatches_3p = integer(),
    reverse_mismatches_3p = integer(),
    duplicated_circular_hit = logical(), stringsAsFactors = FALSE
  )
}

enumerate_primer_amplicons <- function(sites, references, config) {
  products <- list()
  for (reference_id in unique(sites$reference_id)) {
    reference <- references[
      references$reference_id == reference_id,
      ,
      drop = FALSE
    ]
    reference_sites <- sites[
      sites$reference_id == reference_id & sites$passes_3p,
      ,
      drop = FALSE
    ]
    forward_sites <- reference_sites[reference_sites$primer_id == "forward", ]
    reverse_sites <- reference_sites[reference_sites$primer_id == "reverse", ]
    if (!nrow(forward_sites) || !nrow(reverse_sites)) {
      next
    }
    subject <- extend_reference_for_pcr(
      reference,
      config$max_product_size,
      max(forward_sites$raw_end - forward_sites$raw_start + 1L)
    )
    original_length <- reference$source_length[[1]]
    for (fi in seq_len(nrow(forward_sites))) {
      for (ri in seq_len(nrow(reverse_sites))) {
        forward_site <- forward_sites[fi, , drop = FALSE]
        reverse_site <- reverse_sites[ri, , drop = FALSE]
        if (forward_site$strand[[1]] == "+" &&
          reverse_site$strand[[1]] == "-") {
          raw_start <- forward_site$raw_start[[1]]
          raw_end <- reverse_site$raw_end[[1]]
          product_orientation <- "+"
        } else if (forward_site$strand[[1]] == "-" &&
          reverse_site$strand[[1]] == "+") {
          raw_start <- reverse_site$raw_start[[1]]
          raw_end <- forward_site$raw_end[[1]]
          product_orientation <- "-"
        } else {
          next
        }
        if (raw_start > raw_end) {
          next
        }
        product_size <- raw_end - raw_start + 1L
        if (reference$topology[[1]] == "circular" &&
          (raw_start > original_length || product_size > original_length)) {
          next
        }
        if (raw_end > length(subject)) {
          next
        }
        products[[length(products) + 1L]] <- data.frame(
          reference_id = reference_id,
          reference_type = reference$reference_type[[1]],
          start = if (reference$topology[[1]] == "circular") {
            (raw_start - 1L) %% original_length + 1L
          } else raw_start,
          end = if (reference$topology[[1]] == "circular") {
            (raw_end - 1L) %% original_length + 1L
          } else raw_end,
          product_size = product_size,
          sequence = if (product_size < config$min_product_size ||
                         product_size > config$max_product_size) {
            NA_character_
          } else if (product_orientation == "+") {
            as.character(subseq(subject, raw_start, raw_end))
          } else {
            as.character(reverseComplement(subseq(subject, raw_start, raw_end)))
          },
          circular_wrap = raw_end > original_length,
          invalid_size = product_size < config$min_product_size ||
            product_size > config$max_product_size,
          forward_site_id = forward_site$site_id[[1]],
          reverse_site_id = reverse_site$site_id[[1]],
          forward_mismatches_3p = forward_site$mismatches_3p[[1]],
          reverse_mismatches_3p = reverse_site$mismatches_3p[[1]],
          duplicated_circular_hit = FALSE,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(products)) {
    return(empty_amplicon_table())
  }
  products <- bind_rows(products)
  key <- paste(
    products$reference_id,
    products$start,
    products$end,
    products$product_size,
    sep = ":"
  )
  products$duplicated_circular_hit <- duplicated(key)
  products[!products$duplicated_circular_hit, , drop = FALSE]
}

match_probe_pair_all_references <- function(
  forward_probe,
  reverse_probe,
  references,
  config = primer_qc_defaults()
) {
  products <- list()
  max_primer_length <- max(nchar(c(forward_probe, reverse_probe)))
  for (i in seq_len(nrow(references))) {
    reference <- references[i, , drop = FALSE]
    subject <- extend_reference_for_pcr(
      reference,
      config$max_product_size,
      max_primer_length
    )
    views <- matchProbePair(
      Fprobe = DNAString(forward_probe),
      Rprobe = DNAString(reverse_probe),
      subject = subject,
      algorithm = "auto",
      max.mismatch = config$max_mismatches,
      with.indels = FALSE,
      fixed = FALSE
    )
    if (!length(views)) {
      next
    }
    original_length <- reference$source_length[[1]]
    for (j in seq_along(views)) {
      raw_start <- start(views)[[j]]
      raw_end <- end(views)[[j]]
      product_size <- width(views)[[j]]
      if (reference$topology[[1]] == "circular" &&
        (raw_start > original_length || product_size > original_length)) {
        next
      }
      products[[length(products) + 1L]] <- data.frame(
        reference_id = reference$reference_id[[1]],
        reference_type = reference$reference_type[[1]],
        start = if (reference$topology[[1]] == "circular") {
          (raw_start - 1L) %% original_length + 1L
        } else raw_start,
        end = if (reference$topology[[1]] == "circular") {
          (raw_end - 1L) %% original_length + 1L
        } else raw_end,
        product_size = product_size,
        sequence = as.character(views[[j]]),
        circular_wrap = raw_end > original_length,
        invalid_size = product_size < config$min_product_size ||
          product_size > config$max_product_size,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!length(products)) {
    return(empty_amplicon_table()[, c(
      "reference_id", "reference_type", "start", "end", "product_size",
      "sequence", "circular_wrap", "invalid_size"
    )])
  }
  products <- bind_rows(products)
  key <- paste(
    products$reference_id,
    products$start,
    products$end,
    products$product_size,
    sep = ":"
  )
  products[!duplicated(key), , drop = FALSE]
}

expected_product_match <- function(amplicons, expected_product, config) {
  if (!nrow(amplicons)) {
    return(logical())
  }
  coordinate_tolerance <- expected_product$coordinate_tolerance %||%
    config$expected_coordinate_tolerance
  size_tolerance <- expected_product$size_tolerance %||%
    config$expected_size_tolerance
  start_range <- expected_product$start_range %||% c(
    expected_product$start - coordinate_tolerance,
    expected_product$start + coordinate_tolerance
  )
  end_range <- expected_product$end_range %||% c(
    expected_product$end - coordinate_tolerance,
    expected_product$end + coordinate_tolerance
  )
  size_range <- expected_product$size_range %||% c(
    expected_product$size - size_tolerance,
    expected_product$size + size_tolerance
  )
  amplicons$reference_id == expected_product$reference_id &
    amplicons$start >= min(start_range) &
    amplicons$start <= max(start_range) &
    amplicons$end >= min(end_range) &
    amplicons$end <= max(end_range) &
    amplicons$product_size >= min(size_range) &
    amplicons$product_size <= max(size_range) &
    !amplicons$invalid_size
}

evaluate_pair_specificity <- function(
  forward_probe,
  reverse_probe,
  references,
  expected_product,
  config = primer_qc_defaults(),
  binding_site_cache = NULL,
  reference_key = NULL
) {
  if (!is.null(binding_site_cache) && is.null(reference_key)) {
    reference_key <- reference_digest(references)
  }
  find_sites <- function(probe, primer_id) {
    if (is.null(binding_site_cache)) {
      return(enumerate_primer_binding_sites(probe, primer_id, references, config))
    }
    key <- paste0("binding_sites:", digest::digest(list(
      reference_key, toupper(as.character(probe)), primer_id,
      config[c("max_mismatches", "critical_3p_bases", "max_3p_mismatches",
               "max_product_size")]
    )))
    if (!exists(key, binding_site_cache, inherits = FALSE)) {
      assign(key, enumerate_primer_binding_sites(probe, primer_id, references, config),
             binding_site_cache)
    }
    get(key, binding_site_cache, inherits = FALSE)
  }
  sites <- bind_rows(
    find_sites(forward_probe, "forward"),
    find_sites(reverse_probe, "reverse")
  )
  amplicons <- enumerate_primer_amplicons(sites, references, config)
  reduced <- match_probe_pair_all_references(
    forward_probe,
    reverse_probe,
    references,
    config
  )
  if (nrow(amplicons)) {
    reduced_key <- paste(
      reduced$reference_id,
      reduced$start,
      reduced$end,
      reduced$product_size,
      sep = ":"
    )
    amplicon_key <- paste(
      amplicons$reference_id,
      amplicons$start,
      amplicons$end,
      amplicons$product_size,
      sep = ":"
    )
    amplicons$in_match_probe_pair <- amplicon_key %in% reduced_key
    amplicons$intended <- expected_product_match(
      amplicons,
      expected_product,
      config
    )
    amplicons$off_target <- !amplicons$intended
    amplicons$rejection_reason <- ifelse(
      amplicons$intended,
      "",
      ifelse(
        amplicons$invalid_size,
        "invalid_product_size",
        "off_target_product"
      )
    )
  } else {
    amplicons$in_match_probe_pair <- logical()
    amplicons$intended <- logical()
    amplicons$off_target <- logical()
    amplicons$rejection_reason <- character()
  }
  intended_count <- sum(amplicons$intended)
  allowed_expected <- expected_product$allowed_products %||% 1L
  valid_offtargets <- amplicons$off_target & !amplicons$invalid_size
  offtarget_count <- sum(valid_offtargets)
  high_risk <- valid_offtargets &
    amplicons$forward_mismatches_3p == 0L &
    amplicons$reverse_mismatches_3p == 0L
  offtarget_site_ids <- unique(c(
    amplicons$forward_site_id[valid_offtargets],
    amplicons$reverse_site_id[valid_offtargets]
  ))
  perfect_3p_sites <- sites$site_id %in% offtarget_site_ids &
    sites$mismatches_3p == 0L
  reasons <- character()
  if (intended_count < 1L) {
    reasons <- c(reasons, "expected_product_not_found")
  }
  if (intended_count > allowed_expected) {
    reasons <- c(reasons, "too_many_expected_products")
  }
  if (offtarget_count > config$max_allowed_offtarget_products) {
    reasons <- c(reasons, sprintf("off_target_products=%d", offtarget_count))
  }
  report_site_key <- paste(
    sites$primer_id,
    sites$reference_id,
    sites$start,
    sites$end,
    sites$strand,
    sites$mismatches,
    sites$mismatch_positions,
    sep = ":"
  )
  list(
    passed = !length(reasons),
    rejection_reason = paste(reasons, collapse = ";"),
    binding_sites = sites[!duplicated(report_site_key), , drop = FALSE],
    amplicons = amplicons,
    match_probe_pair_amplicons = reduced,
    n_high_risk_offtarget_products = sum(high_risk),
    n_all_offtarget_products = offtarget_count,
    n_perfect_3p_offtarget_sites = sum(perfect_3p_sites),
    n_expected_products = intended_count,
    n_reduced_matches_missed = sum(!amplicons$in_match_probe_pair)
  )
}

assert_openprimer_constraints <- function(active_constraints, required_constraints) {
  missing <- setdiff(required_constraints, active_constraints)
  if (length(missing)) {
    stop(
      sprintf(
        paste(
          "openPrimeR отключил обязательные constraints: %s.",
          "Проверьте внешние программы до запуска primer QC"
        ),
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

configure_openprimer_environment <- function() {
  current <- Sys.getenv("UNAFOLDDAT", unset = "")
  if (nzchar(current) && file.exists(file.path(current, "stack.DAT"))) {
    return(invisible(current))
  }
  hybrid_min <- Sys.which("hybrid-min")
  candidates <- unique(c(
    if (nzchar(hybrid_min)) {
      file.path(dirname(dirname(hybrid_min)), "share", "oligoarrayaux")
    },
    "/usr/local/share/oligoarrayaux"
  ))
  candidates <- candidates[
    file.exists(file.path(candidates, "stack.DAT")) &
      file.exists(file.path(candidates, "miscloop.DAT"))
  ]
  if (length(candidates)) {
    Sys.setenv(UNAFOLDDAT = candidates[[1]])
    return(invisible(candidates[[1]]))
  }
  invisible("")
}

validate_openprimer_tools <- function(required_constraints) {
  configure_openprimer_environment()
  constraint_tools <- list(
    melting_temp_range = c(MELTING = "melting-batch"),
    melting_temp_diff = c(MELTING = "melting-batch"),
    self_dimerization = c(OligoArrayAux = "hybrid-min"),
    cross_dimerization = c(OligoArrayAux = "hybrid-min"),
    secondary_structure = c(ViennaRNA = "RNAfold"),
    primer_efficiency = c(OligoArrayAux = "hybrid-min"),
    annealing_DeltaG = c(OligoArrayAux = "hybrid-min")
  )
  tools <- unique(unlist(constraint_tools[required_constraints]))
  missing <- tools[!nzchar(Sys.which(unname(tools)))]
  if (length(missing)) {
    stop(
      sprintf(
        "Для обязательного openPrimeR QC недоступны: %s",
        paste(paste(names(missing), unname(missing), sep = "="), collapse = ", ")
      ),
      call. = FALSE
    )
  }
  if ("hybrid-min" %in% unname(tools)) {
    status <- suppressWarnings(system2(
      Sys.which("hybrid-min"),
      c(
        "-n", "DNA", "-t", "50", "-T", "50", "-N", "0.05", "-E",
        "-q", "ACAGGTGCCCACTCCCAGGTGCAG",
        "CTGCACCTGGGAGTGGGCACCTGT"
      ),
      stdout = FALSE,
      stderr = FALSE
    ))
    if (!identical(status, 0L)) {
      stop(
        paste(
          "OligoArrayAux найден, но hybrid-min не работает.",
          "Проверьте UNAFOLDDAT"
        ),
        call. = FALSE
      )
    }
  }
  invisible(tools)
}

load_openprimer_settings <- function(
  config = primer_qc_defaults(),
  buffer = primer3_buffer_parameters()
) {
  configure_openprimer_environment()
  if (!requireNamespace("openPrimeR", quietly = TRUE)) {
    stop("R-пакет openPrimeR не установлен", call. = FALSE)
  }
  profile <- system.file(
    "extdata",
    "settings",
    config$openprimer_profile,
    package = "openPrimeR"
  )
  if (!nzchar(profile)) {
    stop(
      sprintf("Не найден профиль openPrimeR: %s", config$openprimer_profile),
      call. = FALSE
    )
  }
  settings <- openPrimeR::read_settings(profile)
  required <- unique(c(
    config$openprimer_critical_annealing,
    config$openprimer_critical_full
  ))
  active <- names(openPrimeR::constraints(settings))
  assert_openprimer_constraints(active, required)
  validate_openprimer_tools(c(required, "primer_efficiency"))

  coverage_constraints <- list(
    primer_efficiency = c(min = config$primer_efficiency_min)
  )
  settings <- openPrimeR::`cvg_constraints<-`(
    settings,
    coverage_constraints
  )
  options <- openPrimeR::conOptions(settings)
  options$allowed_mismatches <- config$max_mismatches
  settings <- openPrimeR::`conOptions<-`(settings, options)
  pcr <- openPrimeR::PCR(settings)
  pcr$Na_concentration <- 0
  pcr$K_concentration <- buffer[["monovalent_salt_mm"]] / 1000
  pcr$Mg_concentration <- buffer[["divalent_salt_mm"]] / 1000
  pcr$primer_concentration <- buffer[["dna_nm"]] * 1e-9
  settings <- openPrimeR::`PCR<-`(settings, pcr)
  list(
    settings = settings,
    active_constraints = active,
    required_constraints = required,
    profile = profile,
    pcr = pcr,
    unit_notes = c(
      monovalent_salt = "Primer3 mM -> openPrimeR mol/L (K)",
      divalent_salt = "Primer3 mM -> openPrimeR mol/L (Mg)",
      primer_dna = "Primer3 nM -> openPrimeR mol/L",
      dntp = "Primer3 dNTP has no openPrimeR PCR field"
    )
  )
}

make_openprimer_objects <- function(forward, reverse, template, id) {
  primer_path <- tempfile("2pac-openprimer-primers-", fileext = ".fasta")
  template_path <- tempfile("2pac-openprimer-template-", fileext = ".fasta")
  on.exit(unlink(c(primer_path, template_path)), add = TRUE)
  writeXStringSet(
    DNAStringSet(setNames(c(forward, reverse), paste0(id, c("_fw", "_rev")))),
    primer_path
  )
  writeXStringSet(DNAStringSet(setNames(template, paste0(id, "_template"))), template_path)
  list(
    primers = openPrimeR::read_primers(primer_path),
    templates = openPrimeR::read_templates(
      template_path,
      fw.region = c(1L, nchar(template)),
      rev.region = c(1L, nchar(template))
    )
  )
}

openprimer_sequence_forms <- function(
  reaction,
  annealing_forward,
  annealing_reverse,
  full_forward,
  full_reverse
) {
  data.frame(
    reaction = rep(reaction, 2L),
    primer_id = c("forward", "reverse"),
    annealing_sequence = c(annealing_forward, annealing_reverse),
    full_oligo_sequence = c(full_forward, full_reverse),
    stringsAsFactors = FALSE
  )
}

with_openprimer_locale <- function(expression) {
  old_monetary_locale <- Sys.getlocale("LC_MONETARY")
  monetary_locale <- suppressWarnings(Sys.setlocale("LC_MONETARY", "en_US.UTF-8"))
  if (is.na(monetary_locale) || !nzchar(Sys.localeconv()[["mon_decimal_point"]])) {
    stop(
      paste(
        "openPrimeR/MELTING требует locale с monetary decimal point;",
        "не удалось включить en_US.UTF-8"
      ),
      call. = FALSE
    )
  }
  on.exit(Sys.setlocale("LC_MONETARY", old_monetary_locale), add = TRUE)
  force(expression)
}

evaluate_openprimer_form <- function(
  forward,
  reverse,
  template_sequence,
  settings,
  active_constraints,
  identifier
) {
  objects <- make_openprimer_objects(
    forward,
    reverse,
    template_sequence,
    identifier
  )
  qc <- with_openprimer_locale(openPrimeR::check_constraints(
    objects$primers,
    objects$templates,
    settings,
    active.constraints = active_constraints
  ))
  qc_frame <- as.data.frame(qc, stringsAsFactors = FALSE)
  metrics <- data.frame(
    lapply(qc_frame, identity),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  passed <- isTRUE(metrics$constraints_passed[[1]])
  metrics$constraints_passed <- NULL
  score <- openPrimeR::score_primers(
    qc,
    settings,
    active.constraints = active_constraints
  )
  list(
    metrics = metrics,
    passed = passed,
    penalty = score$Penalty[[1]],
    active_constraints = active_constraints
  )
}

combine_openprimer_form_results <- function(
  annealing_result,
  full_result,
  sequence_forms,
  loaded,
  config,
  reaction
) {
  annealing_metrics <- annealing_result$metrics
  full_metrics <- if (is.null(full_result)) data.frame() else full_result$metrics
  input_columns <- c(
    "Identifier", "ID", "Forward", "Reverse", "primer_length_fw",
    "primer_length_rev", "Direction", "Degeneracy_fw", "Degeneracy_rev", "Run"
  )
  if (ncol(full_metrics)) {
    full_metrics <- full_metrics[
      ,
      setdiff(names(full_metrics), input_columns),
      drop = FALSE
    ]
    duplicate_columns <- intersect(names(full_metrics), names(annealing_metrics))
    names(full_metrics)[names(full_metrics) %in% duplicate_columns] <- paste0(
      "full_",
      names(full_metrics)[names(full_metrics) %in% duplicate_columns]
    )
    metrics <- cbind(annealing_metrics, full_metrics)
  } else {
    metrics <- annealing_metrics
  }
  eval_columns <- grep("^EVAL_", names(metrics), value = TRUE)
  soft_eval <- intersect(
    paste0("EVAL_", config$openprimer_soft_constraints),
    eval_columns
  )
  dimer_values <- unlist(metrics[, intersect(
    c("Self_Dimer_DeltaG", "Cross_Dimer_DeltaG", "Structure_deltaG"),
    names(metrics)
  ), drop = FALSE])
  dimer_values <- suppressWarnings(as.numeric(dimer_values))
  dimer_values <- dimer_values[is.finite(dimer_values)]
  evaluated_constraints <- c(
    annealing_result$active_constraints,
    if (!is.null(full_result)) full_result$active_constraints
  )
  metrics$reaction <- reaction
  metrics$annealing_forward <- sequence_forms$annealing_sequence[[1]]
  metrics$annealing_reverse <- sequence_forms$annealing_sequence[[2]]
  metrics$full_forward <- sequence_forms$full_oligo_sequence[[1]]
  metrics$full_reverse <- sequence_forms$full_oligo_sequence[[2]]
  metrics$annealing_constraints_passed <- annealing_result$passed
  metrics$full_constraints_passed <- !is.null(full_result) && full_result$passed
  metrics$constraints_passed <- annealing_result$passed &&
    !is.null(full_result) && full_result$passed
  metrics$penalty <- annealing_result$penalty + if (is.null(full_result)) {
    0
  } else full_result$penalty
  metrics$openprimer_failed_soft_constraints <- if (length(soft_eval)) {
    sum(!unlist(metrics[1, soft_eval, drop = FALSE]))
  } else 0L
  metrics$max_dimer_risk <- if (length(dimer_values)) {
    max(0, -min(dimer_values))
  } else 0
  metrics$abs_tm_diff <- abs(metrics$melting_temp_diff[[1]])
  metrics$unavailable_constraints <- ""
  metrics$skipped_constraints <- paste(
    setdiff(loaded$active_constraints, evaluated_constraints),
    collapse = ","
  )
  failed <- eval_columns[!unlist(metrics[1, eval_columns, drop = FALSE])]
  list(
    passed = isTRUE(metrics$constraints_passed[[1]]),
    rejection_reason = if (isTRUE(metrics$constraints_passed[[1]])) "" else {
      paste0("openprimer_failed:", paste(failed, collapse = ","))
    },
    metrics = metrics,
    penalty = metrics$penalty[[1]],
    failed_soft_constraints = metrics$openprimer_failed_soft_constraints[[1]],
    max_dimer_risk = metrics$max_dimer_risk[[1]],
    abs_tm_diff = metrics$abs_tm_diff[[1]],
    unavailable_constraints = character(),
    skipped_constraints = if (nzchar(metrics$skipped_constraints[[1]])) {
      strsplit(metrics$skipped_constraints[[1]], ",", fixed = TRUE)[[1]]
    } else character()
  )
}

evaluate_openprimer_pair <- function(
  annealing_forward,
  annealing_reverse,
  full_forward,
  full_reverse,
  template_sequence,
  config = primer_qc_defaults(),
  buffer = primer3_buffer_parameters(),
  loaded_settings = NULL,
  reaction = "reaction"
) {
  loaded <- loaded_settings %||% load_openprimer_settings(config, buffer)
  sequence_forms <- openprimer_sequence_forms(
    reaction,
    annealing_forward,
    annealing_reverse,
    full_forward,
    full_reverse
  )
  annealing_constraints <- unique(c(
    config$openprimer_critical_annealing,
    "primer_coverage"
  ))
  full_constraints <- config$openprimer_critical_full
  annealing_result <- evaluate_openprimer_form(
    annealing_forward,
    annealing_reverse,
    template_sequence,
    loaded$settings,
    annealing_constraints,
    paste0(reaction, "_annealing")
  )
  full_result <- NULL
  if (annealing_result$passed) {
    full_template <- paste0(
      full_forward,
      paste(rep("A", 20L), collapse = ""),
      as.character(reverseComplement(DNAString(full_reverse)))
    )
    full_result <- evaluate_openprimer_form(
      full_forward,
      full_reverse,
      full_template,
      loaded$settings,
      full_constraints,
      paste0(reaction, "_full")
    )
  }
  combine_openprimer_form_results(
    annealing_result,
    full_result,
    sequence_forms,
    loaded,
    config,
    reaction
  )
}

unavailable_openprimer_result <- function(
  reaction,
  reason,
  abs_tm_diff = Inf,
  tm_source = "unavailable"
) {
  reason <- if (inherits(reason, "condition")) conditionMessage(reason) else {
    as.character(reason)
  }
  reason <- gsub("[\r\n\t]+", " ", reason)
  list(
    passed = FALSE,
    rejection_reason = paste0("openprimer_unavailable:", reason),
    metrics = data.frame(
      reaction = reaction,
      constraints_passed = FALSE,
      penalty = Inf,
      melting_temp_diff = abs_tm_diff,
      tm_source = tm_source,
      unavailable_constraints = reason,
      stringsAsFactors = FALSE
    ),
    penalty = Inf,
    failed_soft_constraints = 0L,
    max_dimer_risk = Inf,
    abs_tm_diff = abs_tm_diff,
    unavailable_constraints = reason,
    skipped_constraints = character()
  )
}

evaluate_filtering_policy <- function(
  specificity,
  openprimer,
  filtering_level,
  config = primer_qc_defaults()
) {
  mode <- filtering_level_name(filtering_level)
  expected_ok <- identical(as.integer(specificity$n_expected_products), 1L)
  high_risk_ok <- specificity$n_high_risk_offtarget_products == 0L
  tm_diff <- if (is.null(openprimer)) Inf else openprimer$abs_tm_diff
  tm_ok <- length(tm_diff) == 1L && is.finite(tm_diff) &&
    tm_diff <= config$hard_max_tm_diff
  strict_passed <- isTRUE(specificity$passed) &&
    !is.null(openprimer) && isTRUE(openprimer$passed)

  blocking_reasons <- switch(
    as.character(filtering_level),
    `1` = character(),
    `2` = if (!expected_ok) {
      paste0("expected_product_count=", specificity$n_expected_products)
    } else character(),
    `3` = c(
      if (!expected_ok) {
        paste0("expected_product_count=", specificity$n_expected_products)
      },
      if (!high_risk_ok) {
        paste0(
          "high_risk_off_target_products=",
          specificity$n_high_risk_offtarget_products
        )
      },
      if (!tm_ok) {
        if (is.finite(tm_diff)) {
          paste0("tm_difference=", tm_diff, ">", config$hard_max_tm_diff)
        } else "tm_difference_unavailable"
      }
    )
  )
  warnings <- unique(c(
    if (!isTRUE(specificity$passed)) specificity$rejection_reason,
    if (specificity$n_high_risk_offtarget_products > 0L) {
      paste0(
        "high_risk_off_target_products=",
        specificity$n_high_risk_offtarget_products
      )
    },
    if (specificity$n_perfect_3p_offtarget_sites > 0L) {
      paste0(
        "perfect_3p_off_target_sites=",
        specificity$n_perfect_3p_offtarget_sites
      )
    },
    if (is.null(openprimer)) "openprimer_not_evaluated" else if (
      !isTRUE(openprimer$passed)
    ) openprimer$rejection_reason
  ))
  warnings <- warnings[!is.na(warnings) & nzchar(warnings)]
  list(
    filtering_level = as.integer(filtering_level),
    filtering_mode = mode,
    strict_passed = strict_passed,
    blocking_passed = !length(blocking_reasons),
    blocking_reasons = paste(blocking_reasons, collapse = ";"),
    warnings = paste(warnings, collapse = ";")
  )
}

select_best_primer_pair <- function(candidates) {
  required_gates <- c(
    "structure_passed",
    "blocking_passed"
  )
  missing_gates <- setdiff(c(required_gates, "strict_qc_passed"), names(candidates))
  if (length(missing_gates)) {
    stop(
      sprintf(
        "Таблица ranking не содержит hard gates: %s",
        paste(missing_gates, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  rank_columns <- c(
    "n_expected_product_deviations",
    "n_high_risk_offtarget_products",
    "n_all_offtarget_products",
    "n_perfect_3p_offtarget_sites",
    "openprimer_failed_soft_constraints",
    "openprimer_penalty",
    "max_dimer_risk",
    "abs_tm_diff",
    "primer3_pair_penalty",
    "deleted_nt",
    "primer3_index"
  )
  for (column in setdiff(rank_columns, names(candidates))) {
    candidates[[column]] <- Inf
  }
  gate_values <- lapply(candidates[required_gates], function(value) {
    !is.na(value) & as.logical(value)
  })
  eligible <- Reduce(`&`, gate_values)
  strict <- eligible & !is.na(candidates$strict_qc_passed) &
    as.logical(candidates$strict_qc_passed)
  selection_pool <- if (any(strict)) strict else eligible
  candidates$selected <- FALSE
  candidates$fallback_selected <- FALSE
  gate_reasons <- apply(
    as.data.frame(gate_values),
    1L,
    function(status) paste(
      sub("_passed$", "", required_gates[!status]),
      collapse = ","
    )
  )
  if (!"rejection_reason" %in% names(candidates)) {
    candidates$rejection_reason <- ""
  }
  candidates$rejection_reason[!eligible] <- ifelse(
    nzchar(candidates$rejection_reason[!eligible]),
    candidates$rejection_reason[!eligible],
    paste0("hard_gate:", gate_reasons[!eligible])
  )
  if (!any(selection_pool)) {
    return(list(pair = NULL, ranking = candidates))
  }
  eligible_indices <- which(selection_pool)
  order_args <- lapply(
    candidates[eligible_indices, rank_columns, drop = FALSE],
    function(values) {
      values <- suppressWarnings(as.numeric(values))
      values[is.na(values)] <- Inf
      values
    }
  )
  selected_index <- eligible_indices[do.call(order, order_args)[[1]]]
  candidates$selected[[selected_index]] <- TRUE
  candidates$fallback_selected[[selected_index]] <- !strict[[selected_index]]
  candidates$rejection_reason[eligible & !selection_pool] <-
    "strict_qc_candidate_available"
  candidates$rejection_reason[selection_pool & !candidates$selected] <-
    "lower_deterministic_rank"
  list(
    pair = candidates[selected_index, , drop = FALSE],
    ranking = candidates
  )
}

output_layout <- function(output_dir) {
  list(
    wet_lab = file.path(output_dir, "WetLab"),
    tech_report = file.path(output_dir, "TechReport")
  )
}

write_run_parameters <- function(input, targets, path) {
  openprimer_key <- ".openprimer_settings"
  loaded_openprimer <- if (exists(
    openprimer_key,
    input$primer_qc_cache,
    inherits = FALSE
  )) {
    get(openprimer_key, input$primer_qc_cache, inherits = FALSE)
  } else NULL
  openprimer_limits <- character()
  if (!is.null(loaded_openprimer)) {
    selected_limits <- openPrimeR::constraints(loaded_openprimer$settings)[
      loaded_openprimer$required_constraints
    ]
    openprimer_limits <- unlist(lapply(
      names(selected_limits),
      function(constraint) {
        setNames(
          selected_limits[[constraint]],
          paste0(
            "openprimer_constraint_",
            constraint,
            "_",
            names(selected_limits[[constraint]])
          )
        )
      }
    ))
  }
  parameters <- c(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    genome_file = input$genome_path,
    genome_annotation_file = input$annotation_path,
    annotation_format = input$annotation_format,
    target_plasmid_file = input$target_plasmid,
    ptarget_site1 = input$parameters$site1,
    ptarget_site2 = input$parameters$site2,
    cas_plasmid_file = if (is.null(input$cas_plasmid)) NA_character_ else {
      input$cas_plasmid
    },
    cds_targets = paste(targets$gene[targets$class == "cds"], collapse = ","),
    ncrna_targets = paste(
      targets$gene[targets$class == "ncrna"],
      collapse = ","
    ),
    output_directory = input$output_dir,
    genome_index_directory = file.path(
      output_layout(input$output_dir)$tech_report,
      "genome_indexes"
    ),
    chopchop_script = input$tools$chopchop_script,
    chopchop_python = input$tools$chopchop_python,
    chopchop_config_snapshot = file.path(
      output_layout(input$output_dir)$tech_report,
      "chopchop_config.json"
    ),
    primer3_executable = input$tools$primer3,
    primer3_thermodynamic_parameters = input$tools$primer3_config,
    filtering_level = input$parameters$filtering_level,
    filtering_mode = filtering_level_name(input$parameters$filtering_level),
    n20_count = input$parameters$n20_mn,
    n20_strands = input$parameters$n20_strands,
    n20_offtarget_thresholds = paste(
      input$parameters$n20_offtarget,
      collapse = ","
    ),
    cds_frame_restriction = input$parameters$cds_fs,
    ncrna_frame_restriction = input$parameters$ncrna_fs,
    setNames(
      input$parameters$left_arm,
      paste0("left_arm_", names(input$parameters$left_arm), "_nt")
    ),
    setNames(
      input$parameters$right_arm,
      paste0("right_arm_", names(input$parameters$right_arm), "_nt")
    ),
    n20_arm_min_distance_nt = input$parameters$n20_arm_min_distance,
    setNames(
      input$parameters$primer3_buffer,
      paste0("primer3_buffer_", names(input$parameters$primer3_buffer))
    ),
    primer_qc_max_mismatches = input$parameters$primer_qc$max_mismatches,
    primer_qc_critical_3p_bases = input$parameters$primer_qc$critical_3p_bases,
    primer_qc_max_3p_mismatches =
      input$parameters$primer_qc$max_3p_mismatches,
    primer_qc_min_product_size = input$parameters$primer_qc$min_product_size,
    primer_qc_max_product_size = input$parameters$primer_qc$max_product_size,
    primer_qc_max_allowed_offtarget_products =
      input$parameters$primer_qc$max_allowed_offtarget_products,
    primer_qc_hard_max_tm_diff =
      input$parameters$primer_qc$hard_max_tm_diff,
    primer_qc_primer_efficiency_min =
      input$parameters$primer_qc$primer_efficiency_min,
    openprimer_profile = input$parameters$primer_qc$openprimer_profile,
    openprimer_active_constraints = if (is.null(loaded_openprimer)) NA else {
      paste(loaded_openprimer$active_constraints, collapse = ",")
    },
    openprimer_required_constraints = if (is.null(loaded_openprimer)) NA else {
      paste(loaded_openprimer$required_constraints, collapse = ",")
    },
    openprimer_load_error = if (exists(
      ".openprimer_load_error",
      input$primer_qc_cache,
      inherits = FALSE
    )) get(
      ".openprimer_load_error",
      input$primer_qc_cache,
      inherits = FALSE
    ) else NA_character_,
    openprimer_limits,
    openprimer_version = tryCatch(
      as.character(utils::packageVersion("openPrimeR")),
      error = function(e) NA_character_
    ),
    decipher_version = tryCatch(
      as.character(utils::packageVersion("DECIPHER")),
      error = function(e) NA_character_
    ),
    biostrings_version = as.character(packageVersion("Biostrings")),
    melting_executable = Sys.which("melting-batch"),
    viennarna_executable = Sys.which("RNAfold"),
    oligoarrayaux_executable = Sys.which("hybrid-min"),
    mafft_executable = Sys.which("mafft")
  )
  write_tsv(
    data.frame(
      parameter = names(parameters),
      value = unname(as.character(parameters)),
      stringsAsFactors = FALSE
    ),
    path
  )
  invisible(path)
}

make_design_input <- function(cli) {
  project_tool <- function(path, default) {
    if (identical(path, default)) file.path(.dhole_project_dir, path) else path
  }
  # TODO: only a complete single-contig bacterial genome is supported for now.
  # Multi-contig FASTA/GFF input needs contig-aware sequence extraction.
  genome_set <- readDNAStringSet(cli$genome[[1]], nrec = 1)
  references <- make_specificity_references(
    cli$genome[[1]],
    cli$target_plasmid[[1]],
    if (length(cli$cas_plasmid)) cli$cas_plasmid[[1]] else NULL
  )
  target_records <- which(references$reference_type == "target_plasmid")
  if (length(target_records) != 1L) {
    stop(
      "Для моделирования требуется ровно одна запись pTarget в FASTA",
      call. = FALSE
    )
  }
  find_oriented_restriction_pair(
    references$sequence[[target_records[[1]]]],
    cli$site1,
    cli$site2
  )
  input <- list(
    genome_path = cli$genome[[1]],
    annotation_path = cli$genome_annotation[[1]],
    annotation_format = cli$annotation_format[[1]],
    genome = genome_set[[1]],
    genome_contig = names(genome_set)[[1]],
    annotation = read_genome_annotation(
      cli$genome_annotation[[1]],
      cli$annotation_format[[1]]
    ),
    target_plasmid = cli$target_plasmid[[1]],
    target_plasmid_sequence = references$sequence[[target_records[[1]]]],
    target_plasmid_name = references$source_name[[target_records[[1]]]],
    cas_plasmid = if (length(cli$cas_plasmid)) cli$cas_plasmid[[1]] else NULL,
    output_dir = cli$output_dir[[1]],
    tools = list(
      chopchop_script = project_tool(cli$chopchop_script[[1]], "chopchop/chopchop.py"),
      chopchop_python = cli$chopchop_python[[1]],
      primer3 = project_tool(cli$primer3[[1]], "primer3/src/primer3_core"),
      primer3_config = file.path(.dhole_project_dir, "primer3/src/primer3_config")
    ),
    parameters = list(
      filtering_level = cli$filtering_level,
      n20_mn = cli$n20_mn,
      n20_strands = cli$n20_strands,
      n20_offtarget = cli$n20_offtarget,
      site1 = cli$site1,
      site2 = cli$site2,
      cds_fs = cli$cds_fs,
      ncrna_fs = cli$ncrna_fs,
      left_arm = cli$left_arm,
      right_arm = cli$right_arm,
      n20_arm_min_distance = cli$n20_arm_min_distance,
      primer3_buffer = primer3_buffer_parameters(),
      primer_qc = cli$primer_qc
    ),
    specificity_references = references,
    specificity_reference_digest = reference_digest(references),
    primer_qc_cache = new.env(parent = emptyenv())
  )
  genome_ids <- references$reference_id[
    references$reference_type == "genome" &
      references$source_name == input$genome_contig
  ]
  input$genome_reference_id <- if (length(genome_ids)) {
    genome_ids[[1]]
  } else {
    references$reference_id[references$reference_type == "genome"][[1]]
  }
  input
}

feature_record <- function(input, gene_name) {
  idx <- find_target_feature(input$annotation, gene_name)
  row <- input$annotation[idx, , drop = FALSE]
  contig <- as.character(row$seqid[[1]])
  candidates <- input$annotation
  if (!is.na(contig)) {
    candidates <- filter(candidates, seqid == contig)
  }
  previous_stop <- suppressWarnings(max(
    candidates$stop[candidates$stop < row$start[[1]]],
    na.rm = TRUE
  ))
  next_start <- suppressWarnings(min(
    candidates$start[candidates$start > row$stop[[1]]],
    na.rm = TRUE
  ))
  if (!is.finite(previous_stop)) {
    previous_stop <- 0L
  }
  if (!is.finite(next_start)) {
    next_start <- length(input$genome) + 1L
  }
  list(
    index = idx,
    query_name = gene_name,
    display_name = if (!is.na(row$gene[[1]]) && nzchar(row$gene[[1]])) {
      row$gene[[1]]
    } else {
      gene_name
    },
    contig = if (!is.na(contig)) contig else input$genome_contig,
    start = as.integer(row$start[[1]]),
    end = as.integer(row$stop[[1]]),
    strand = as.character(row$strand[[1]]),
    length = as.integer(row$stop[[1]] - row$start[[1]] + 1),
    left_bound = as.integer(previous_stop + 1L),
    right_bound = as.integer(next_start - 1L)
  )
}

cut_interval <- function(feature, design_class, genome_length) {
  if (feature$length <= 500) {
    if (design_class == "ncrna") {
      extension <- round(feature$length / 10, -1) * 2
      interval <- c(feature$start - extension, feature$end + extension)
    } else {
      interval <- c(feature$start, feature$end)
    }
  } else if (feature$length < 1500) {
    midpoint <- feature$start + feature$length %/% 2
    interval <- c(midpoint - 250, midpoint + 250)
  } else {
    flank <- floor((feature$end - feature$start) / 3)
    interval <- c(feature$start + flank, feature$end - flank)
  }
  interval <- pmax(1L, pmin(as.integer(interval), genome_length))
  if (design_class == "ncrna") {
    interval[[1]] <- max(interval[[1]], feature$left_bound)
    interval[[2]] <- min(interval[[2]], feature$right_bound)
  }
  interval
}

prepare_chopchop_assets <- function(input) {
  genome_name <- tools::file_path_sans_ext(basename(input$genome_path))
  index_dir <- file.path(
    output_layout(input$output_dir)$tech_report,
    "genome_indexes"
  )
  dir.create(index_dir, recursive = TRUE, showWarnings = FALSE)
  index_dir <- normalizePath(index_dir)
  two_bit <- file.path(index_dir, paste0(genome_name, ".2bit"))
  bowtie_prefix <- file.path(index_dir, genome_name)
  index_files <- c(two_bit, paste0(bowtie_prefix, c(
    ".1.ebwt", ".2.ebwt", ".3.ebwt", ".4.ebwt", ".rev.1.ebwt", ".rev.2.ebwt"
  )))
  checksum_file <- paste0(bowtie_prefix, ".fasta.md5")
  checksum <- unname(tools::md5sum(input$genome_path))
  if (is.na(checksum)) stop("Не удалось прочитать FASTA для индексации", call. = FALSE)
  reusable <- file.exists(checksum_file) &&
    identical(readLines(checksum_file, warn = FALSE), checksum) &&
    all(file.exists(index_files)) && all(file.size(index_files) > 0L)
  if (!reusable) {
    # A marker is written only after both tools completed and every file exists.
    unlink(c(checksum_file, index_files))
    run_tool("faToTwoBit", c(input$genome_path, two_bit))
    run_tool("bowtie-build", c(input$genome_path, bowtie_prefix))
    if (!all(file.exists(index_files)) || any(file.size(index_files) == 0L)) {
      stop("Индексация генома не создала полный комплект файлов", call. = FALSE)
    }
    writeLines(checksum, checksum_file)
  }
  list(name = genome_name, directory = index_dir)
}

configure_chopchop <- function(input, genome_assets) {
  script_path <- normalizePath(input$tools$chopchop_script)
  chopchop_dir <- dirname(script_path)
  two_bit_to_fa <- Sys.which("twoBitToFa")
  bowtie <- Sys.which("bowtie")
  if (!nzchar(two_bit_to_fa) || !nzchar(bowtie)) {
    stop("twoBitToFa или bowtie не найдены в PATH", call. = FALSE)
  }
  config <- file.path(chopchop_dir, "config_local.json")
  writeLines(
    c(
      "{",
      "  \"PATH\": {",
      paste0("    \"PRIMER3\": \"", normalizePath(input$tools$primer3), "\","),
      paste0("    \"BOWTIE\": \"", bowtie, "\","),
      paste0("    \"TWOBITTOFA\": \"", two_bit_to_fa, "\","),
      paste0("    \"TWOBIT_INDEX_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"BOWTIE_INDEX_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"ISOFORMS_INDEX_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"ISOFORMS_MT_DIR\": \"", genome_assets$directory, "\","),
      paste0("    \"GENE_TABLE_INDEX_DIR\": \"", genome_assets$directory, "\""),
      "  },",
      "  \"THREADS\": 1",
      "}"
    ),
    config
  )
  invisible(config)
}

run_chopchop <- function(
  input,
  genome_name,
  feature,
  design_class,
  target_dir
) {
  interval <- cut_interval(feature, design_class, length(input$genome))
  target <- paste0(feature$contig, ":", interval[[1]], "-", interval[[2]])
  table_path <- file.path(target_dir, "n20_table.tsv")
  run_tool(
    input$tools$chopchop_python,
    c(
      input$tools$chopchop_script,
      "-Target",
      target,
      "-G",
      genome_name,
      "-M",
      "NGG",
      "-T",
      "1",
      "-g",
      "20",
      "-m",
      "3",
      "--padSize",
      "0",
      "-O",
      "50",
      "--scoringMethod",
      "DOENCH_2016",
      "-o",
      target_dir
    ),
    stdout = table_path,
    stderr = file.path(target_dir, "chopchop.stderr.log")
  )
  if (!file.exists(table_path) || file.size(table_path) == 0) {
    stop("ChopChop did not produce n20_table.tsv", call. = FALSE)
  }
  offtarget_files <- list.files(
    target_dir,
    pattern = "\\.offtargets$",
    full.names = TRUE
  )
  if (length(offtarget_files)) {
    offtarget_dir <- file.path(target_dir, "chopchop_offtargets")
    dir.create(offtarget_dir, showWarnings = FALSE)
    file.rename(
      offtarget_files,
      file.path(offtarget_dir, basename(offtarget_files))
    )
  }
  unlink(file.path(
    target_dir,
    c(
      "sequence.fa",
      "gene_file.fa",
      "output.sam",
      "bowtie.err",
      "twoBitToFa.err"
    )
  ))
  interval
}

filter_grnas <- function(
  table_path,
  feature,
  design_class,
  offtarget_thresholds,
  genome = NULL
) {
  grnas <- read_tsv(table_path, show_col_types = FALSE) |>
    janitor::clean_names()
  required <- c(
    "genomic_location",
    "target_sequence",
    "strand",
    "self_complementarity"
  )
  missing <- setdiff(required, names(grnas))
  if (length(missing)) {
    stop(
      sprintf(
        "Таблица CHOPCHOP не содержит обязательные колонки: %s",
        paste(missing, collapse = ", ")
      ),
      call. = FALSE
    )
  }
  mm_columns <- grep("^mm[0-9]+$", names(grnas), value = TRUE)
  mm_columns <- mm_columns[order(as.integer(sub("^mm", "", mm_columns)))]
  if (length(offtarget_thresholds) > length(mm_columns)) {
    stop(
      sprintf(
        paste(
          "--n20-offtarget содержит %d значений,",
          "но CHOPCHOP вывел только %d колонок MM"
        ),
        length(offtarget_thresholds),
        length(mm_columns)
      ),
      call. = FALSE
    )
  }
  applied_mm_columns <- head(mm_columns, length(offtarget_thresholds))
  keep_offtarget <- rep(TRUE, nrow(grnas))
  for (i in seq_along(applied_mm_columns)) {
    values <- suppressWarnings(as.numeric(grnas[[applied_mm_columns[[i]]]]))
    keep_offtarget <- keep_offtarget &
      !is.na(values) &
      values <= offtarget_thresholds[[i]]
  }
  grnas <- grnas[keep_offtarget, , drop = FALSE] |>
    mutate(
      genomic_start = as.numeric(sub("^.*:", "", genomic_location)),
      genomic_end = genomic_start + 22L,
      n20_start = genomic_start + ifelse(strand == "-", 3L, 0L),
      n20_end = n20_start + 19L,
      pam_start = genomic_start + ifelse(strand == "-", 0L, 20L),
      pam_end = pam_start + 2L,
      mid_closeness = abs(
        (feature$start + feature$end) %/%
          2 -
          (n20_start + n20_end) %/% 2
      ) /
        feature$length
    ) |>
    filter(self_complementarity == 0)
  if (anyNA(grnas$genomic_start)) {
    stop(
      "Не удалось разобрать genomic_location в таблице CHOPCHOP",
      call. = FALSE
    )
  }
  if (any(!grnas$strand %in% c("+", "-"))) {
    stop("Некорректная цепь N20 в таблице CHOPCHOP", call. = FALSE)
  }
  if (!is.null(genome) && nrow(grnas)) {
    if (any(grnas$genomic_start < 1L | grnas$genomic_end > length(genome))) {
      stop("Участок N20/PAM выходит за границы генома", call. = FALSE)
    }
    matches <- vapply(seq_len(nrow(grnas)), function(i) {
      sequence <- genome[grnas$genomic_start[[i]]:grnas$genomic_end[[i]]]
      if (grnas$strand[[i]] == "-") sequence <- reverseComplement(sequence)
      identical(as.character(sequence), toupper(grnas$target_sequence[[i]]))
    }, logical(1))
    if (!all(matches)) {
      stop("Последовательность N20/PAM CHOPCHOP не совпадает с геномом", call. = FALSE)
    }
  }
  if (design_class == "cds") {
    grnas <- filter(grnas, mid_closeness <= 0.18)
  }
  arrange(grnas, mid_closeness)
}

prepare_grna_pool <- function(grnas, n20_mn, strand_mode) {
  if (strand_mode == "plus") {
    grnas <- filter(grnas, strand == "+")
  } else if (strand_mode == "minus") {
    grnas <- filter(grnas, strand == "-")
  }
  if (nrow(grnas) < n20_mn) {
    stop(
      sprintf(
        "Недостаточно подходящих N20: найдено %d, требуется не менее %d",
        nrow(grnas),
        n20_mn
      ),
      call. = FALSE
    )
  }
  if (
    n20_mn > 1L &&
      strand_mode == "both" &&
      (!any(grnas$strand == "+") || !any(grnas$strand == "-"))
  ) {
    stop(
      "Для режима both нужны подходящие N20 на плюс- и минус-цепях",
      call. = FALSE
    )
  }
  grnas
}

visit_grna_sets <- function(grnas, n20_mn, strand_mode, visitor) {
  selected_indices <- integer(n20_mn)
  visit <- function(depth, first_index) {
    if (depth > n20_mn) {
      candidates <- grnas[selected_indices, , drop = FALSE]
      if (
        n20_mn > 1L &&
          strand_mode == "both" &&
          !all(c("+", "-") %in% candidates$strand)
      ) {
        return(NULL)
      }
      return(visitor(candidates))
    }
    remaining <- n20_mn - depth
    last_index <- nrow(grnas) - remaining
    for (index in seq.int(first_index, last_index)) {
      selected_indices[[depth]] <<- index
      result <- visit(depth + 1L, index + 1L)
      if (!is.null(result)) {
        return(result)
      }
    }
    NULL
  }
  visit(1L, 1L)
}

write_selected_grnas <- function(selected, feature, target_dir) {
  write_tsv(
    selected$table,
    file.path(target_dir, "selected_n20_table.tsv")
  )
  n20 <- DNAStringSet(substr(selected$table$target_sequence, 1, 20))
  names(n20) <- paste0(feature$display_name, "_n20_", seq_along(n20))
  writeXStringSet(n20, file.path(target_dir, "selected_n20.fasta"))
  invisible(selected)
}

write_primer3_settings <- function(
  path,
  left_length,
  right_length,
  buffer = primer3_buffer_parameters()
) {
  product <- round(c(
    min(left_length, right_length),
    max(left_length, right_length) * 1.5
  ))
  writeLines(
    c(
      "Primer3 File - http://primer3.org",
      "P3_FILE_TYPE=settings",
      "",
      "PRIMER_TASK=generic",
      "PRIMER_PICK_LEFT_PRIMER=1",
      "PRIMER_PICK_RIGHT_PRIMER=1",
      "PRIMER_NUM_RETURN=10",
      "PRIMER_MIN_SIZE=18",
      "PRIMER_OPT_SIZE=21",
      "PRIMER_MAX_SIZE=27",
      "PRIMER_MIN_TM=59.0",
      "PRIMER_OPT_TM=60.0",
      "PRIMER_MAX_TM=61.0",
      paste0("PRIMER_SALT_MONOVALENT=", buffer[["monovalent_salt_mm"]]),
      paste0("PRIMER_SALT_DIVALENT=", buffer[["divalent_salt_mm"]]),
      paste0("PRIMER_DNTP_CONC=", buffer[["dntp_mm"]]),
      paste0("PRIMER_DNA_CONC=", buffer[["dna_nm"]]),
      "PRIMER_PAIR_MAX_DIFF_TM=8.0",
      "PRIMER_MIN_GC=40.0",
      "PRIMER_MAX_GC=60.0",
      "PRIMER_MAX_SELF_ANY=12.0",
      "PRIMER_MAX_SELF_END=8.0",
      "PRIMER_PAIR_MAX_COMPL_ANY=12.0",
      "PRIMER_PAIR_MAX_COMPL_END=8.0",
      "PRIMER_MAX_POLY_X=5",
      paste0("PRIMER_PRODUCT_SIZE_RANGE=", product[[1]], "-", product[[2]]),
      "PRIMER_EXPLAIN_FLAG=1",
      "PRIMER_FIRST_BASE_INDEX=1",
      "="
    ),
    path
  )
}

add_genome_positions <- function(left, right, arm, plus_strand) {
  if (plus_strand) {
    left <- mutate(
      left,
      genome_start = PRIMER_LEFT_pos + arm$left_start - 1,
      genome_end = PRIMER_RIGHT_pos + arm$left_start - 1
    )
    right <- mutate(
      right,
      genome_start = PRIMER_LEFT_pos + arm$right_start - 1,
      genome_end = PRIMER_RIGHT_pos + arm$right_start - 1
    )
  } else {
    left <- mutate(
      left,
      genome_start = arm$left_end - PRIMER_RIGHT_pos + 1,
      genome_end = arm$left_end - PRIMER_LEFT_pos + 1
    )
    right <- mutate(
      right,
      genome_start = arm$right_end - PRIMER_RIGHT_pos + 1,
      genome_end = arm$right_end - PRIMER_LEFT_pos + 1
    )
  }
  list(left = left, right = right)
}

structural_pair_candidates <- function(
  primers,
  plus_strand,
  feature,
  selected,
  minimum_n20_distance,
  restrict_frame_shift
) {
  left <- primers$left
  right <- primers$right
  candidates <- list()
  if (!nrow(left) || !nrow(right)) {
    return(data.frame(
      left_index = integer(), right_index = integer(),
      deleted_nt = integer(), left_distance = integer(),
      right_distance = integer(), stringsAsFactors = FALSE
    ))
  }
  for (left_index in seq_len(nrow(left))) {
    for (right_index in seq_len(nrow(right))) {
      if (plus_strand) {
        lower_boundary <- left$genome_end[[left_index]]
        upper_boundary <- right$genome_start[[right_index]]
      } else {
        lower_boundary <- right$genome_end[[right_index]]
        upper_boundary <- left$genome_start[[left_index]]
      }
      deleted_nt <- upper_boundary - lower_boundary - 1L
      left_distance <- selected$n20_range[[1]] - lower_boundary - 1L
      right_distance <- upper_boundary - selected$n20_range[[2]] - 1L
      eligible <- deleted_nt >= 0L &&
        deleted_nt < feature$length &&
        left_distance >= minimum_n20_distance &&
        right_distance >= minimum_n20_distance &&
        (!restrict_frame_shift || deleted_nt %% 3L == 0L)
      if (eligible) {
        candidates[[length(candidates) + 1L]] <- data.frame(
          left_index = left_index,
          right_index = right_index,
          deleted_nt = deleted_nt,
          left_distance = left_distance,
          right_distance = right_distance,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  if (!length(candidates)) {
    return(data.frame(
      left_index = integer(), right_index = integer(),
      deleted_nt = integer(), left_distance = integer(),
      right_distance = integer(), stringsAsFactors = FALSE
    ))
  }
  bind_rows(candidates)
}

choose_pair <- function(
  primers,
  plus_strand,
  feature,
  selected,
  minimum_n20_distance,
  restrict_frame_shift
) {
  left <- primers$left
  right <- primers$right
  candidates <- structural_pair_candidates(
    primers,
    plus_strand,
    feature,
    selected,
    minimum_n20_distance,
    restrict_frame_shift
  )
  if (!nrow(candidates)) {
    return(NULL)
  }
  selected_index <- order(
    candidates$deleted_nt,
    candidates$left_index,
    candidates$right_index
  )[[1]]
  candidate <- candidates[selected_index, , drop = FALSE]
  bind_rows(
    left[candidate$left_index[[1]], , drop = FALSE],
    right[candidate$right_index[[1]], , drop = FALSE]
  )
}

new_primer_qc_trace <- function() {
  trace <- new.env(parent = emptyenv())
  trace$binding_sites <- list()
  trace$amplicons <- list()
  trace$openprimer <- list()
  trace$ranking <- list()
  trace
}

append_primer_qc_trace <- function(
  trace,
  reaction,
  pair_id,
  specificity = NULL,
  openprimer = NULL,
  ranking = NULL
) {
  if (!is.null(specificity)) {
    sites <- specificity$binding_sites
    sites$reaction <- rep(reaction, nrow(sites))
    sites$pair_id <- rep(pair_id, nrow(sites))
    trace$binding_sites[[length(trace$binding_sites) + 1L]] <- sites
    amplicons <- specificity$amplicons
    amplicons$reaction <- rep(reaction, nrow(amplicons))
    amplicons$pair_id <- rep(pair_id, nrow(amplicons))
    trace$amplicons[[length(trace$amplicons) + 1L]] <- amplicons
  }
  if (!is.null(openprimer)) {
    metrics <- openprimer$metrics
    metrics$pair_id <- rep(pair_id, nrow(metrics))
    trace$openprimer[[length(trace$openprimer) + 1L]] <- metrics
  }
  if (!is.null(ranking)) {
    if (!"selected_for_attempt" %in% names(ranking)) {
      ranking$selected_for_attempt <- ranking$selected
    }
    trace$ranking[[length(trace$ranking) + 1L]] <- ranking
  }
  invisible(trace)
}

write_primer_qc_trace <- function(trace, target_dir) {
  bind_or_empty <- function(records, empty, label) {
    if (!length(records)) {
      return(empty)
    }
    tryCatch(
      bind_rows(records),
      error = function(e) {
        stop(
          sprintf("Не удалось собрать primer QC trace '%s': %s", label, e$message),
          call. = FALSE
        )
      }
    )
  }
  sites <- bind_or_empty(
    trace$binding_sites,
    data.frame(
      reaction = character(), pair_id = character(), primer_id = character(),
      reference_id = character(), reference_type = character(),
      start = integer(), end = integer(), strand = character(),
      mismatches = integer(), mismatch_positions = character(),
      mismatches_3p = integer(), stringsAsFactors = FALSE
    ),
    "binding_sites"
  )
  site_columns <- c(
    "reaction", "pair_id", "primer_id", "reference_id", "reference_type",
    "start", "end", "strand", "mismatches", "mismatch_positions",
    "mismatches_3p", "passes_3p", "circular_wrap"
  )
  for (column in setdiff(site_columns, names(sites))) {
    sites[[column]] <- rep(NA, nrow(sites))
  }
  write_tsv(
    sites[, site_columns, drop = FALSE],
    file.path(target_dir, "primer_binding_sites.tsv")
  )

  amplicons <- bind_or_empty(
    trace$amplicons,
    data.frame(
      reaction = character(), pair_id = character(), reference_id = character(),
      start = integer(), end = integer(), product_size = integer(),
      intended = logical(), off_target = logical(), circular_wrap = logical(),
      rejection_reason = character(), stringsAsFactors = FALSE
    ),
    "amplicons"
  )
  amplicon_columns <- c(
    "reaction", "pair_id", "reference_id", "reference_type", "start", "end",
    "product_size", "intended", "off_target", "invalid_size",
    "circular_wrap", "in_match_probe_pair", "rejection_reason"
  )
  for (column in setdiff(amplicon_columns, names(amplicons))) {
    amplicons[[column]] <- rep(NA, nrow(amplicons))
  }
  write_tsv(
    amplicons[, amplicon_columns, drop = FALSE],
    file.path(target_dir, "primer_amplicons.tsv")
  )

  openprimer <- bind_or_empty(
    trace$openprimer,
    data.frame(
      reaction = character(), pair_id = character(),
      constraints_passed = logical(), penalty = numeric(),
      stringsAsFactors = FALSE
    ),
    "openprimer"
  )
  list_columns <- names(openprimer)[vapply(openprimer, is.list, logical(1))]
  for (column in list_columns) {
    openprimer[[column]] <- vapply(
      openprimer[[column]],
      function(value) paste(value, collapse = ","),
      character(1)
    )
  }
  text_columns <- names(openprimer)[vapply(
    openprimer,
    function(column) is.character(column) || is.factor(column),
    logical(1)
  )]
  for (column in text_columns) {
    openprimer[[column]] <- gsub(
      "[\r\n\t]+",
      " ",
      as.character(openprimer[[column]])
    )
  }
  write_tsv(
    openprimer,
    file.path(target_dir, "primer_openprimer_qc.tsv")
  )

  ranking <- bind_or_empty(
    trace$ranking,
    data.frame(
      pair_id = character(), reaction = character(),
      primer3_index = integer(), structure_passed = logical(),
      specificity_passed = logical(), openprimer_passed = logical(),
      selected = logical(), rejection_reason = character(),
      stringsAsFactors = FALSE
    ),
    "ranking"
  )
  write_tsv(ranking, file.path(target_dir, "primer_pair_ranking.tsv"))
  invisible(target_dir)
}

cached_pair_specificity <- function(input, forward, reverse, expected_product) {
  config <- input$parameters$primer_qc
  key <- paste0(
    "specificity:",
    digest::digest(list(
      input$specificity_reference_digest,
      toupper(forward),
      toupper(reverse),
      expected_product,
      config[c(
        "max_mismatches",
        "critical_3p_bases",
        "max_3p_mismatches",
        "min_product_size",
        "max_product_size",
        "max_allowed_offtarget_products",
        "expected_coordinate_tolerance",
        "expected_size_tolerance"
      )]
    ))
  )
  if (!exists(key, input$primer_qc_cache, inherits = FALSE)) {
    assign(
      key,
      evaluate_pair_specificity(
        forward,
        reverse,
        input$specificity_references,
        expected_product,
        config,
        input$primer_qc_cache,
        input$specificity_reference_digest
      ),
      input$primer_qc_cache
    )
  }
  get(key, input$primer_qc_cache, inherits = FALSE)
}

cached_openprimer_pair <- function(
  input,
  annealing_forward,
  annealing_reverse,
  full_forward,
  full_reverse,
  template_sequence,
  reaction
) {
  config <- input$parameters$primer_qc
  settings_key <- ".openprimer_settings"
  if (!exists(settings_key, input$primer_qc_cache, inherits = FALSE)) {
    if (exists(".openprimer_load_error", input$primer_qc_cache, inherits = FALSE)) {
      stop(get(
        ".openprimer_load_error",
        input$primer_qc_cache,
        inherits = FALSE
      ), call. = FALSE)
    }
    loaded <- tryCatch(
      load_openprimer_settings(config, input$parameters$primer3_buffer),
      error = function(e) {
        stop(structure(
          list(
            message = conditionMessage(e),
            call = NULL,
            stage = "primer_qc",
            parent = e
          ),
          class = c("primer_qc_error", "error", "condition")
        ))
      }
    )
    assign(settings_key, loaded, input$primer_qc_cache)
  }
  loaded <- get(settings_key, input$primer_qc_cache, inherits = FALSE)
  sequence_forms <- openprimer_sequence_forms(
    reaction,
    annealing_forward,
    annealing_reverse,
    full_forward,
    full_reverse
  )
  annealing_constraints <- unique(c(
    config$openprimer_critical_annealing,
    "primer_coverage"
  ))
  annealing_key <- paste0(
    "openprimer_annealing:",
    digest::digest(list(
      toupper(c(annealing_forward, annealing_reverse, template_sequence)),
      annealing_constraints,
      input$parameters$primer3_buffer
    ))
  )
  if (!exists(annealing_key, input$primer_qc_cache, inherits = FALSE)) {
    assign(
      annealing_key,
      evaluate_openprimer_form(
        annealing_forward,
        annealing_reverse,
        template_sequence,
        loaded$settings,
        annealing_constraints,
        paste0(reaction, "_annealing")
      ),
      input$primer_qc_cache
    )
  }
  annealing_result <- get(
    annealing_key,
    input$primer_qc_cache,
    inherits = FALSE
  )
  full_result <- NULL
  if (annealing_result$passed) {
    full_constraints <- config$openprimer_critical_full
    full_key <- paste0(
      "openprimer_full:",
      digest::digest(list(
        toupper(c(full_forward, full_reverse)),
        full_constraints,
        input$parameters$primer3_buffer
      ))
    )
    if (!exists(full_key, input$primer_qc_cache, inherits = FALSE)) {
      full_template <- paste0(
        full_forward,
        paste(rep("A", 20L), collapse = ""),
        as.character(reverseComplement(DNAString(full_reverse)))
      )
      assign(
        full_key,
        evaluate_openprimer_form(
          full_forward,
          full_reverse,
          full_template,
          loaded$settings,
          full_constraints,
          paste0(reaction, "_full")
        ),
        input$primer_qc_cache
      )
    }
    full_result <- get(full_key, input$primer_qc_cache, inherits = FALSE)
  }
  combine_openprimer_form_results(
    annealing_result,
    full_result,
    sequence_forms,
    loaded,
    config,
    reaction
  )
}

expected_product_from_primer3 <- function(input, primer_row) {
  size <- primer_row$genome_end[[1]] - primer_row$genome_start[[1]] + 1L
  list(
    reference_id = input$genome_reference_id,
    start = as.integer(primer_row$genome_start[[1]]),
    end = as.integer(primer_row$genome_end[[1]]),
    size = as.integer(size),
    allowed_products = 1L
  )
}

evaluate_candidate_reaction <- function(
  input,
  primer_row,
  full_forward,
  full_reverse,
  reaction,
  pair_id,
  trace
) {
  forward <- primer_row$PRIMER_LEFT_SEQUENCE[[1]]
  reverse <- primer_row$PRIMER_RIGHT_SEQUENCE[[1]]
  expected <- expected_product_from_primer3(input, primer_row)
  specificity <- cached_pair_specificity(input, forward, reverse, expected)
  append_primer_qc_trace(
    trace,
    reaction,
    pair_id,
    specificity = specificity
  )
  intended <- specificity$amplicons[specificity$amplicons$intended, , drop = FALSE]
  primer_tm <- suppressWarnings(as.numeric(c(
    primer_row$PRIMER_LEFT_TM[[1]],
    primer_row$PRIMER_RIGHT_TM[[1]]
  )))
  primer3_tm_diff <- if (length(primer_tm) == 2L && all(is.finite(primer_tm))) {
    abs(diff(primer_tm))
  } else Inf
  openprimer <- if (nrow(intended) == 1L && !is.na(intended$sequence[[1]])) {
    tryCatch(
      cached_openprimer_pair(
        input,
        forward,
        reverse,
        full_forward,
        full_reverse,
        intended$sequence[[1]],
        reaction
      ),
      error = function(e) unavailable_openprimer_result(
        reaction,
        e,
        primer3_tm_diff,
        "Primer3 fallback"
      )
    )
  } else {
    unavailable_openprimer_result(
      reaction,
      "expected product is unavailable or ambiguous",
      primer3_tm_diff,
      "Primer3 fallback"
    )
  }
  append_primer_qc_trace(
    trace,
    reaction,
    pair_id,
    openprimer = openprimer
  )
  policy <- evaluate_filtering_policy(
    specificity,
    openprimer,
    input$parameters$filtering_level,
    input$parameters$primer_qc
  )
  list(
    passed = policy$blocking_passed,
    specificity = specificity,
    openprimer = openprimer,
    policy = policy,
    rejection_reason = policy$blocking_reasons,
    risk_warnings = policy$warnings
  )
}

design_homology_arms <- function(
  input,
  feature,
  selected,
  design_class,
  target_dir,
  trace = new_primer_qc_trace(),
  log_path = NULL,
  n20_attempt = 1L
) {
  interval <- cut_interval(feature, design_class, length(input$genome))
  parameters <- input$parameters
  left_limits <- parameters$left_arm
  right_limits <- parameters$right_arm
  minimum_n20_distance <- parameters$n20_arm_min_distance
  restrict_frame_shift <- if (design_class == "cds") {
    parameters$cds_fs
  } else {
    parameters$ncrna_fs
  }
  settings <- file.path(target_dir, "primer3_settings.txt")
  write_primer3_settings(
    settings,
    left_limits[["max"]],
    right_limits[["max"]],
    input$parameters$primer3_buffer
  )
  if (!file.exists(input$tools$primer3)) {
    stop("primer3_core не найден", call. = FALSE)
  }

  attempt <- 0L
  deferred_fallbacks <- list()
  repeat {
    attempt <- attempt + 1L
    left_length <- min(
      left_limits[["opt"]] + attempt - 1L,
      left_limits[["max"]]
    )
    right_length <- min(
      right_limits[["opt"]] + attempt - 1L,
      right_limits[["max"]]
    )
    if (feature$strand == "+") {
      left_end <- interval[[1]] - 1L
      left_start <- max(1L, left_end - left_length + 1L)
      right_start <- interval[[2]] + 1L
      right_end <- min(
        length(input$genome),
        right_start + right_length - 1L
      )
      left_seq <- input$genome[left_start:left_end]
      right_seq <- input$genome[right_start:right_end]
    } else {
      right_end <- interval[[1]] - 1L
      left_start <- interval[[2]] + 1L
      left_end <- min(
        length(input$genome),
        left_start + left_length - 1L
      )
      right_start <- max(1L, right_end - right_length + 1L)
      left_seq <- complement(reverse(input$genome[left_start:left_end]))
      right_seq <- complement(reverse(input$genome[right_start:right_end]))
    }
    if (
      length(left_seq) < left_limits[["min"]] ||
        length(right_seq) < right_limits[["min"]]
    ) {
      stop(
        "Граница генома не позволяет выдержать минимальную длину плеч",
        call. = FALSE
      )
    }
    writeXStringSet(
      DNAStringSet(c(left_arm = left_seq, right_arm = right_seq)),
      file.path(target_dir, "homology_arms_before_primer_search.fasta")
    )
    tm <- if (design_class == "cds") c(62.5, 63, 63.5) else c(60.5, 61, 62.5)
    left <- callPrimer3(
      as.character(left_seq),
      paste0(left_limits[["min"]], "-", length(left_seq)),
      tm,
      2,
      "left_arm",
      primer_num = 10,
      primer3 = input$tools$primer3,
      thermo.param = input$tools$primer3_config,
      settings = settings,
      report = file.path(target_dir, "left_arm_report.txt")
    )
    right <- callPrimer3(
      as.character(right_seq),
      paste0(right_limits[["min"]], "-", length(right_seq)),
      tm,
      2,
      "right_arm",
      primer_num = 10,
      primer3 = input$tools$primer3,
      thermo.param = input$tools$primer3_config,
      settings = settings,
      report = file.path(target_dir, "right_arm_report.txt")
    )
    if (is.data.frame(left) && nrow(left) && is.data.frame(right) && nrow(right)) {
      arm <- list(
        left_start = left_start,
        left_end = left_end,
        right_start = right_start,
        right_end = right_end
      )
      positions <- add_genome_positions(left, right, arm, feature$strand == "+")
      combinations <- structural_pair_candidates(
        positions,
        feature$strand == "+",
        feature,
        selected,
        minimum_n20_distance,
        restrict_frame_shift
      )
      ranking_rows <- list()
      if (nrow(combinations)) {
        combinations$bridge_mod <- combinations$deleted_nt %% 3L
        combinations$structure_passed <- vapply(
          seq_len(nrow(combinations)),
          function(i) {
            combination <- combinations[i, , drop = FALSE]
            pair <- bind_rows(
              positions$left[combination$left_index[[1]], , drop = FALSE],
              positions$right[combination$right_index[[1]], , drop = FALSE]
            )
            pair_arm_lengths <- pair$PRIMER_RIGHT_pos -
              pair$PRIMER_LEFT_pos +
              1L
            ticks <- sort(unlist(pair[, c("genome_start", "genome_end")]))
            pair_arm_lengths[[1]] >= left_limits[["min"]] &&
              pair_arm_lengths[[1]] <= left_limits[["max"]] &&
              pair_arm_lengths[[2]] >= right_limits[["min"]] &&
              pair_arm_lengths[[2]] <= right_limits[["max"]] &&
              ticks[[2]] >= feature$start &&
              ticks[[3]] <= feature$end
          },
          logical(1)
        )
        reaction_results <- new.env(parent = emptyenv())
        evaluate_physical_pair <- function(side, primer_index, bridge_mod) {
          key <- paste(side, primer_index, bridge_mod, sep = ":")
          if (exists(key, reaction_results, inherits = FALSE)) {
            return(invisible(NULL))
          }
          primer_row <- positions[[side]][primer_index, , drop = FALSE]
          bridge <- design_bridge(design_class, bridge_mod)
          reaction <- if (side == "left") "LF_LR" else "RF_RR"
          physical_pair_id <- sprintf(
            "n20_%03d_attempt_%03d_%s_%02d_B%d",
            n20_attempt,
            attempt,
            reaction,
            primer_index,
            bridge_mod
          )
          if (!is.null(log_path)) {
            append_design_log(
              log_path,
              "primer_qc",
              "TRY",
              sprintf("pair_id=%s;reaction=%s", physical_pair_id, reaction)
            )
          }
          full_primers <- homology_primer_pair(
            primer_row, side, bridge, input$parameters$site2
          )
          result <- evaluate_candidate_reaction(
            input, primer_row, full_primers[[1]], full_primers[[2]],
            reaction, physical_pair_id, trace
          )
          assign(
            key,
            list(result = result, pair_id = physical_pair_id),
            reaction_results
          )
          invisible(NULL)
        }
        eligible <- combinations[combinations$structure_passed, , drop = FALSE]
        if (nrow(eligible)) {
          left_physical <- unique(eligible[, c("left_index", "bridge_mod")])
          right_physical <- unique(eligible[, c("right_index", "bridge_mod")])
          for (i in seq_len(nrow(left_physical))) {
            evaluate_physical_pair(
              "left",
              left_physical$left_index[[i]],
              left_physical$bridge_mod[[i]]
            )
          }
          for (i in seq_len(nrow(right_physical))) {
            evaluate_physical_pair(
              "right",
              right_physical$right_index[[i]],
              right_physical$bridge_mod[[i]]
            )
          }
        }
        for (combination_index in seq_len(nrow(combinations))) {
          combination <- combinations[combination_index, , drop = FALSE]
          left_row <- positions$left[
            combination$left_index[[1]],
            ,
            drop = FALSE
          ]
          right_row <- positions$right[
            combination$right_index[[1]],
            ,
            drop = FALSE
          ]
          pair <- bind_rows(left_row, right_row)
          structure_passed <- combinations$structure_passed[[combination_index]]
          pair_id <- sprintf(
            "n20_%03d_attempt_%03d_L%02d_R%02d",
            n20_attempt,
            attempt,
            combination$left_index[[1]],
            combination$right_index[[1]]
          )
          left_result <- NULL
          right_result <- NULL
          left_pair_id <- NA_character_
          right_pair_id <- NA_character_
          rejection_reason <- ""
          if (structure_passed) {
            bridge_mod <- combination$bridge_mod[[1]]
            left_evaluation <- get(
              paste("left", combination$left_index[[1]], bridge_mod, sep = ":"),
              reaction_results,
              inherits = FALSE
            )
            right_evaluation <- get(
              paste("right", combination$right_index[[1]], bridge_mod, sep = ":"),
              reaction_results,
              inherits = FALSE
            )
            left_result <- left_evaluation$result
            right_result <- right_evaluation$result
            left_pair_id <- left_evaluation$pair_id
            right_pair_id <- right_evaluation$pair_id
            rejection_reason <- paste(
              c(
                if (!left_result$passed) {
                  paste0("LF_LR:", left_result$rejection_reason)
                },
                if (!right_result$passed) {
                  paste0("RF_RR:", right_result$rejection_reason)
                }
              ),
              collapse = ";"
            )
          } else {
            rejection_reason <- "structural_constraints"
          }
          specificity_results <- Filter(
            Negate(is.null),
            list(
              if (!is.null(left_result)) left_result$specificity else NULL,
              if (!is.null(right_result)) right_result$specificity else NULL
            )
          )
          openprimer_results <- Filter(
            Negate(is.null),
            list(
              if (!is.null(left_result)) left_result$openprimer else NULL,
              if (!is.null(right_result)) right_result$openprimer else NULL
            )
          )
          policy_results <- Filter(
            Negate(is.null),
            list(
              if (!is.null(left_result)) left_result$policy else NULL,
              if (!is.null(right_result)) right_result$policy else NULL
            )
          )
          risk_warnings <- paste(c(
            if (!is.null(left_result) && nzchar(left_result$risk_warnings)) {
              paste0("LF_LR:", left_result$risk_warnings)
            },
            if (!is.null(right_result) && nzchar(right_result$risk_warnings)) {
              paste0("RF_RR:", right_result$risk_warnings)
            }
          ), collapse = ";")
          primer3_penalties <- suppressWarnings(as.numeric(c(
            if ("PRIMER_PAIR_PENALTY" %in% names(left_row)) {
              left_row$PRIMER_PAIR_PENALTY[[1]]
            } else Inf,
            if ("PRIMER_PAIR_PENALTY" %in% names(right_row)) {
              right_row$PRIMER_PAIR_PENALTY[[1]]
            } else Inf
          )))
          ranking_rows[[length(ranking_rows) + 1L]] <- data.frame(
            pair_id = pair_id,
            reaction = "homology_arms",
            primer3_index = (attempt - 1L) * 100L + combination_index,
            left_primer3_index = combination$left_index[[1]],
            right_primer3_index = combination$right_index[[1]],
            left_pair_id = left_pair_id,
            right_pair_id = right_pair_id,
            filtering_level = input$parameters$filtering_level,
            filtering_mode = filtering_level_name(
              input$parameters$filtering_level
            ),
            structure_passed = structure_passed,
            specificity_passed = length(specificity_results) == 2L &&
              all(vapply(specificity_results, `[[`, logical(1), "passed")),
            openprimer_passed = length(openprimer_results) == 2L &&
              all(vapply(openprimer_results, `[[`, logical(1), "passed")),
            strict_qc_passed = length(policy_results) == 2L &&
              all(vapply(policy_results, `[[`, logical(1), "strict_passed")),
            blocking_passed = length(policy_results) == 2L &&
              all(vapply(policy_results, `[[`, logical(1), "blocking_passed")),
            n_expected_product_deviations = sum(abs(vapply(
              specificity_results,
              `[[`,
              numeric(1),
              "n_expected_products"
            ) - 1L)),
            n_high_risk_offtarget_products = sum(vapply(
              specificity_results,
              `[[`,
              numeric(1),
              "n_high_risk_offtarget_products"
            )),
            n_all_offtarget_products = sum(vapply(
              specificity_results,
              `[[`,
              numeric(1),
              "n_all_offtarget_products"
            )),
            n_perfect_3p_offtarget_sites = sum(vapply(
              specificity_results,
              `[[`,
              numeric(1),
              "n_perfect_3p_offtarget_sites"
            )),
            openprimer_failed_soft_constraints = sum(vapply(
              openprimer_results,
              `[[`,
              numeric(1),
              "failed_soft_constraints"
            )),
            openprimer_penalty = sum(vapply(
              openprimer_results,
              `[[`,
              numeric(1),
              "penalty"
            )),
            max_dimer_risk = if (length(openprimer_results)) {
              max(vapply(
                openprimer_results,
                `[[`,
                numeric(1),
                "max_dimer_risk"
              ))
            } else Inf,
            abs_tm_diff = if (length(openprimer_results)) {
              max(vapply(
                openprimer_results,
                `[[`,
                numeric(1),
                "abs_tm_diff"
              ))
            } else Inf,
            primer3_pair_penalty = sum(primer3_penalties),
            deleted_nt = combination$deleted_nt[[1]],
            risk_warnings = risk_warnings,
            rejection_reason = rejection_reason,
            stringsAsFactors = FALSE
          )
        }
      }
      if (length(ranking_rows)) {
        selection <- select_best_primer_pair(bind_rows(ranking_rows))
        deferred <- !is.null(selection$pair) &&
          isTRUE(selection$pair$fallback_selected[[1]])
        trace_ranking <- selection$ranking
        if (deferred) {
          trace_ranking$selected_for_attempt <- trace_ranking$selected
          provisional <- trace_ranking$selected
          trace_ranking$selected[provisional] <- FALSE
          trace_ranking$rejection_reason[provisional] <-
            "deferred_fallback_candidate"
        }
        append_primer_qc_trace(
          trace,
          "homology_arms",
          paste0("n20_", n20_attempt, "_attempt_", attempt),
          ranking = trace_ranking
        )
        trace_index <- length(trace$ranking)
        if (!is.null(log_path)) {
          for (i in seq_len(nrow(selection$ranking))) {
            candidate <- selection$ranking[i, , drop = FALSE]
            status <- if (candidate$selected[[1]] && deferred) {
              "DEFERRED"
            } else if (candidate$selected[[1]]) "OK" else "REJECTED"
            append_design_log(
              log_path,
              "primer_qc",
              status,
              sprintf(
                "pair_id=%s;reason=%s",
                candidate$pair_id[[1]],
                if (candidate$selected[[1]] &&
                    isTRUE(candidate$fallback_selected[[1]])) {
                  candidate$risk_warnings[[1]]
                } else candidate$rejection_reason[[1]]
              )
            )
          }
        }
        if (!is.null(selection$pair)) {
          selected_rank <- selection$pair
          pair <- bind_rows(
            positions$left[
              selected_rank$left_primer3_index[[1]],
              ,
              drop = FALSE
            ],
            positions$right[
              selected_rank$right_primer3_index[[1]],
              ,
              drop = FALSE
            ]
          )
          ticks <- sort(unlist(pair[, c("genome_start", "genome_end")]))
          position_table <- bind_rows(positions$left, positions$right)
          design <- list(
            pair = pair,
            left = left_seq,
            right = right_seq,
            arm = arm,
            ticks = ticks,
            selected_pair_id = selected_rank$pair_id[[1]],
            primer_qc_trace = trace
          )
          if (!deferred) {
            write_tsv(
              position_table,
              file.path(target_dir, "primer3_table.tsv")
            )
            return(design)
          }
          design$ranking <- selected_rank
          design$positions <- position_table
          design$primer3_reports <- lapply(
            c("left_arm_report.txt", "right_arm_report.txt"),
            function(report_name) {
              report_path <- file.path(target_dir, report_name)
              if (file.exists(report_path)) readLines(report_path, warn = FALSE) else {
                character()
              }
            }
          )
          design$trace_index <- trace_index
          deferred_fallbacks[[length(deferred_fallbacks) + 1L]] <- design
        }
      }
    }
    interval <- interval + c(1L, -1L)
    required_deletion <- selected$n20_range +
      c(-minimum_n20_distance - 1L, minimum_n20_distance + 1L)
    if (
      interval[[1]] > required_deletion[[1]] ||
        interval[[2]] < required_deletion[[2]] ||
        diff(interval) + 1 < feature$length * .3
    ) {
      break
    }
  }
  if (length(deferred_fallbacks)) {
    fallback_ranking <- bind_rows(lapply(
      deferred_fallbacks,
      `[[`,
      "ranking"
    ))
    final_selection <- select_best_primer_pair(fallback_ranking)$pair
    selected_id <- final_selection$pair_id[[1]]
    fallback_index <- which(vapply(
      deferred_fallbacks,
      function(candidate) identical(candidate$selected_pair_id, selected_id),
      logical(1)
    ))[[1]]
    design <- deferred_fallbacks[[fallback_index]]
    trace_ranking <- trace$ranking[[design$trace_index]]
    selected_row <- trace_ranking$pair_id == selected_id
    trace_ranking$selected[selected_row] <- TRUE
    trace_ranking$rejection_reason[selected_row] <- ""
    trace$ranking[[design$trace_index]] <- trace_ranking
    if (!is.null(log_path)) {
      append_design_log(
        log_path,
        "primer_qc",
        "OK",
        sprintf(
          "pair_id=%s;reason=%s",
          selected_id,
          design$ranking$risk_warnings[[1]]
        )
      )
    }
    write_tsv(design$positions, file.path(target_dir, "primer3_table.tsv"))
    writeXStringSet(
      DNAStringSet(c(left_arm = design$left, right_arm = design$right)),
      file.path(target_dir, "homology_arms_before_primer_search.fasta")
    )
    for (i in seq_along(design$primer3_reports)) {
      if (length(design$primer3_reports[[i]])) {
        writeLines(
          design$primer3_reports[[i]],
          file.path(
            target_dir,
            c("left_arm_report.txt", "right_arm_report.txt")[[i]]
          )
        )
      }
    }
    design$ranking <- NULL
    design$positions <- NULL
    design$primer3_reports <- NULL
    design$trace_index <- NULL
    return(design)
  }
  NULL
}

calculate_screening_product_sizes <- function(
  screening,
  original_genome,
  edited_genome
) {
  reference_size <- if (
    "PRIMER_PAIR_PRODUCT_SIZE" %in% names(screening)
  ) {
    suppressWarnings(as.integer(screening$PRIMER_PAIR_PRODUCT_SIZE[[1]]))
  } else {
    as.integer(
      screening$genome_end[[1]] - screening$genome_start[[1]] + 1L
    )
  }
  edited_genome_size <- if (inherits(edited_genome, "XStringSet")) {
    width(edited_genome)[[1]]
  } else {
    length(edited_genome)
  }
  edited_size <- reference_size + edited_genome_size - length(original_genome)
  if (
    length(reference_size) != 1L ||
      is.na(reference_size) ||
      reference_size < 1L ||
      length(edited_size) != 1L ||
      is.na(edited_size) ||
      edited_size < 1L
  ) {
    stop("Не удалось рассчитать размеры скрининговых ПЦР-продуктов", call. = FALSE)
  }
  c(
    unsuccessful_insertion_bp = reference_size,
    successful_insertion_bp = as.integer(edited_size)
  )
}

calculate_n20_arm_distances <- function(selected, arm_ticks, target_strand) {
  required <- c("target_sequence", "strand", "n20_start", "n20_end")
  if (
    !all(required %in% names(selected)) ||
      length(arm_ticks) != 4L ||
      !target_strand %in% c("+", "-")
  ) {
    stop("Не удалось рассчитать расстояния N20 до плеч гомологии", call. = FALSE)
  }
  left_boundary <- sort(arm_ticks)[[2]]
  right_boundary <- sort(arm_ticks)[[3]]
  genomic_left <- selected$n20_start - left_boundary - 1L
  genomic_right <- right_boundary - selected$n20_end - 1L
  distances <- data.frame(
    name = paste0("N20_", seq_len(nrow(selected))),
    sequence = substr(selected$target_sequence, 1L, 20L),
    strand = selected$strand,
    genomic_coordinates = paste0(selected$n20_start, "-", selected$n20_end),
    left_arm_distance_bp = if (target_strand == "+") {
      genomic_left
    } else genomic_right,
    right_arm_distance_bp = if (target_strand == "+") {
      genomic_right
    } else genomic_left,
    stringsAsFactors = FALSE
  )
  if (!nrow(distances) || anyNA(distances) || any(
    distances$left_arm_distance_bp < 0L |
      distances$right_arm_distance_bp < 0L
  )) {
    stop("N20 находится вне промежутка между плечами гомологии", call. = FALSE)
  }
  distances
}

format_openprimer_report_metrics <- function(metrics) {
  labels <- c(
    constraints_passed = "Все обязательные ограничения пройдены",
    primer_length_fw = "Длина forward-праймера, нт",
    primer_length_rev = "Длина reverse-праймера, нт",
    gc_ratio_fw = "GC-состав forward-праймера, %",
    gc_ratio_rev = "GC-состав reverse-праймера, %",
    gc_clamp_fw = "GC-clamp forward-праймера, нт",
    gc_clamp_rev = "GC-clamp reverse-праймера, нт",
    no_runs_fw = "Максимальный гомополимерный участок forward, нт",
    no_runs_rev = "Максимальный гомополимерный участок reverse, нт",
    no_repeats_fw = "Максимальное число повторов forward",
    no_repeats_rev = "Максимальное число повторов reverse",
    Tm_C_fw = "Tm forward по openPrimeR, °C",
    Tm_C_rev = "Tm reverse по openPrimeR, °C",
    melting_temp_diff = "Разница Tm пары, °C",
    tm_source = "Источник разницы Tm",
    Basic_primer_coverage = "Число покрытых целевых шаблонов",
    Basic_Coverage_Ratio = "Доля покрытых целевых шаблонов, %",
    primer_specificity = "Специфичность праймеров, %",
    primer_efficiency = "Эффективность праймеров",
    mean_primer_efficiency = "Средняя эффективность праймеров",
    Self_Dimer_DeltaG = "Худший self-dimer ΔG, ккал/моль",
    Cross_Dimer_DeltaG = "Худший cross-dimer ΔG, ккал/моль",
    Structure_deltaG_fw = "Вторичная структура forward ΔG, ккал/моль",
    Structure_deltaG_rev = "Вторичная структура reverse ΔG, ккал/моль",
    Structure_deltaG = "Худшая вторичная структура ΔG, ккал/моль",
    penalty = "Суммарный штраф openPrimeR",
    unavailable_constraints = "Недоступные проверки openPrimeR",
    EVAL_primer_length = "Проверка длины праймеров",
    EVAL_gc_ratio = "Проверка GC-состава",
    EVAL_gc_clamp = "Проверка GC-clamp",
    EVAL_no_runs = "Проверка гомополимерных участков",
    EVAL_no_repeats = "Проверка повторов",
    EVAL_melting_temp_range = "Проверка диапазона Tm",
    EVAL_melting_temp_diff = "Проверка разницы Tm",
    EVAL_primer_coverage = "Проверка покрытия шаблона",
    EVAL_primer_efficiency = "Проверка эффективности праймеров",
    EVAL_primer_specificity = "Проверка специфичности праймеров",
    EVAL_self_dimerization = "Проверка self-dimer",
    EVAL_cross_dimerization = "Проверка cross-dimer",
    EVAL_secondary_structure = "Проверка вторичной структуры"
  )
  present <- intersect(names(labels), names(metrics))
  if (!length(present)) {
    stop("Для выбранной screening-пары отсутствуют метрики openPrimeR", call. = FALSE)
  }
  percentage_metrics <- c(
    "gc_ratio_fw", "gc_ratio_rev", "Basic_Coverage_Ratio",
    "primer_specificity"
  )
  format_value <- function(value, metric) {
    value <- unlist(value, use.names = FALSE)
    if (is.logical(value)) {
      return(paste(ifelse(value, "пройдено", "не пройдено"), collapse = ", "))
    }
    numeric_value <- suppressWarnings(as.numeric(value))
    if (metric %in% percentage_metrics && all(!is.na(numeric_value))) {
      numeric_value <- numeric_value * 100
    }
    if (all(!is.na(numeric_value))) {
      return(paste(format(round(numeric_value, 3), trim = TRUE), collapse = ", "))
    }
    paste(gsub("[\r\n\t]+", " ", as.character(value)), collapse = ", ")
  }
  data.frame(
    metric = unname(labels[present]),
    value = vapply(
      present,
      function(metric) format_value(metrics[[metric]], metric),
      character(1)
    ),
    stringsAsFactors = FALSE
  )
}

write_wet_lab_outputs <- function(
  wet_lab_dir,
  feature,
  design_class,
  sequences,
  sequence_purposes,
  primer_metrics,
  screening_product_sizes,
  n20_distances,
  screening_qc,
  edited_genome,
  edited_ptargets,
  ptarget_site_pair,
  pcr_products
) {
  required_metrics <- c("name", "purpose", "annealing_sequence", "tm_c")
  required_distances <- c(
    "name", "sequence", "strand", "genomic_coordinates",
    "left_arm_distance_bp", "right_arm_distance_bp"
  )
  required_screening_qc <- c(
    "pair_id", "filtering_level", "filtering_mode", "selection_status",
    "fallback_used", "warnings", "offtarget_products",
    "high_risk_offtarget_products", "perfect_3p_offtarget_sites",
    "openprimer_metrics"
  )
  required_site_pair <- c(
    "orientation", "site1_start", "site2_start"
  )
  required_pcr_products <- c(
    "name", "description", "location", "length_bp", "primers",
    "pcr_conditions", "sequence"
  )
  if (
    !length(sequences) ||
      length(sequence_purposes) != length(sequences) ||
      is.null(names(sequences)) ||
      anyNA(names(sequences)) ||
      any(!nzchar(names(sequences))) ||
      !all(required_metrics %in% names(primer_metrics)) ||
      !all(required_distances %in% names(n20_distances)) ||
      !all(required_screening_qc %in% names(screening_qc)) ||
      length(edited_genome) != 1L ||
      length(edited_ptargets) != nrow(n20_distances) ||
      is.null(names(edited_genome)) ||
      is.null(names(edited_ptargets)) ||
      !all(required_site_pair %in% names(ptarget_site_pair)) ||
      !nrow(pcr_products) ||
      !all(required_pcr_products %in% names(pcr_products)) ||
      !all(c("metric", "value") %in% names(screening_qc$openprimer_metrics)) ||
      !all(
        c("unsuccessful_insertion_bp", "successful_insertion_bp") %in%
          names(screening_product_sizes)
      )
  ) {
    stop("Неполный набор результатов для WetLab", call. = FALSE)
  }
  dir.create(wet_lab_dir, recursive = TRUE, showWarnings = FALSE)
  writeXStringSet(
    sequences,
    file.path(wet_lab_dir, "final_sequences.fasta")
  )
  sequence_table <- data.frame(
    name = names(sequences),
    sequence = as.character(sequences),
    purpose = unname(sequence_purposes),
    stringsAsFactors = FALSE
  )
  write_tsv(
    sequence_table,
    file.path(wet_lab_dir, "final_sequences.txt")
  )
  writeXStringSet(
    edited_genome,
    file.path(wet_lab_dir, "edited_genome.fasta")
  )
  writeXStringSet(
    edited_ptargets,
    file.path(wet_lab_dir, "edited_pTargets.fasta")
  )
  pcr_product_sequences <- DNAStringSet(pcr_products$sequence)
  names(pcr_product_sequences) <- pcr_products$name
  writeXStringSet(
    pcr_product_sequences,
    file.path(wet_lab_dir, "pcr_products.fasta")
  )
  write_tsv(pcr_products, file.path(wet_lab_dir, "pcr_products.tsv"))

  sequence_lines <- apply(
    sequence_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  tm_table <- primer_metrics
  tm_table$tm_c <- format(
    round(tm_table$tm_c, 1),
    nsmall = 1,
    trim = TRUE
  )
  tm_lines <- apply(
    tm_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  n20_table <- n20_distances[, required_distances, drop = FALSE]
  names(n20_table) <- c(
    "N20", "Последовательность N20 (5'-3')", "Цепь",
    "Геномные координаты", "До левого плеча, п.н.",
    "До правого плеча, п.н."
  )
  n20_lines <- apply(
    n20_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  openprimer_table <- screening_qc$openprimer_metrics
  names(openprimer_table) <- c("Метрика openPrimeR", "Значение")
  openprimer_lines <- apply(
    openprimer_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  construction_table <- data.frame(
    file = c(
      "edited_genome.fasta",
      rep("edited_pTargets.fasta", length(edited_ptargets))
    ),
    record = c(names(edited_genome), names(edited_ptargets)),
    type = c("edited_genome", rep("edited_pTarget", length(edited_ptargets))),
    length_bp = c(width(edited_genome), width(edited_ptargets)),
    stringsAsFactors = FALSE
  )
  construction_lines <- apply(
    construction_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  pcr_table <- pcr_products[, required_pcr_products, drop = FALSE]
  names(pcr_table) <- c(
    "ПЦР-продукт", "Описание", "Локация", "Длина, п.н.",
    "Полные праймеры", "Условия ПЦР", "Последовательность (5'-3')"
  )
  pcr_lines <- apply(
    pcr_table,
    1L,
    function(row) paste(row, collapse = "\t")
  )
  report <- c(
    "2PAC: отчёт для мокрой лаборатории",
    paste("Цель", feature$query_name, sep = "\t"),
    paste("Класс", design_class, sep = "\t"),
    paste(
      "Режим фильтрации праймеров",
      sprintf(
        "%d (%s)",
        screening_qc$filtering_level,
        screening_qc$filtering_mode
      ),
      sep = "\t"
    ),
    paste("Статус primer QC", screening_qc$selection_status, sep = "\t"),
    paste("QC fallback", screening_qc$fallback_used, sep = "\t"),
    paste("Предупреждения primer QC", screening_qc$warnings, sep = "\t"),
    "",
    "Итоговый набор последовательностей",
    paste(names(sequence_table), collapse = "\t"),
    sequence_lines,
    "",
    "Моделированные конструкции",
    paste(names(construction_table), collapse = "\t"),
    construction_lines,
    paste(
      "Ориентация пары сайтов рестрикции в исходной pTarget",
      ptarget_site_pair$orientation,
      sep = "\t"
    ),
    "",
    "PCR-продукты, смоделированные DECIPHER::AmplifyDNA",
    "Последовательности включают полные праймеры со служебными 5'-хвостами.",
    paste(names(pcr_table), collapse = "\t"),
    pcr_lines,
    "",
    paste(
      "Температуры отжига праймеров",
      "Tm рассчитана для участка отжига без сервисных последовательностей",
      sep = "\n"
    ),
    paste(names(tm_table), collapse = "\t"),
    tm_lines,
    "",
    "Размеры скрининговых ПЦР-продуктов",
    paste(
      "Без успешного нокаута (исходный аллель), п.н.",
      screening_product_sizes[["unsuccessful_insertion_bp"]],
      sep = "\t"
    ),
    paste(
      "С успешным нокаутом (редактированный аллель), п.н.",
      screening_product_sizes[["successful_insertion_bp"]],
      sep = "\t"
    ),
    "",
    "Расстояния выбранных N20 до плеч гомологии",
    paste(
      "Левое и правое плечи указаны в ориентации target;",
      "расстояние измерено между ближайшими границами N20 и плеча."
    ),
    paste(names(n20_table), collapse = "\t"),
    n20_lines,
    "",
    "QC скрининговых праймеров",
    paste("Выбранная пара", screening_qc$pair_id, sep = "\t"),
    paste(
      "Оффтаргетные ПЦР-продукты, всего",
      screening_qc$offtarget_products,
      sep = "\t"
    ),
    paste(
      "Высокорисковые оффтаргетные ПЦР-продукты",
      screening_qc$high_risk_offtarget_products,
      sep = "\t"
    ),
    paste(
      "Оффтаргетные сайты с идеальным совпадением 3'-концов",
      screening_qc$perfect_3p_offtarget_sites,
      sep = "\t"
    ),
    "",
    "Метрики качества openPrimeR для выбранной screening-пары",
    paste(names(openprimer_table), collapse = "\t"),
    openprimer_lines
  )
  writeLines(
    enc2utf8(report),
    file.path(wet_lab_dir, "wet_lab_report.txt"),
    useBytes = TRUE
  )
  invisible(wet_lab_dir)
}

write_design_outputs <- function(
  input,
  feature,
  selected,
  arms,
  design_class,
  target_dir,
  log_path = NULL
) {
  trace <- arms$primer_qc_trace
  trace$output_attempts <- (trace$output_attempts %||% 0L) + 1L
  output_id <- paste0(arms$selected_pair_id, "_outputs_", trace$output_attempts)
  completed <- FALSE
  on.exit({
    if (!completed) {
      for (i in seq_along(trace$ranking)) trace$ranking[[i]]$selected <- FALSE
    }
  }, add = TRUE)
  pair <- arms$pair
  gap <- arms$ticks[[3]] - arms$ticks[[2]] - 1L
  bridge <- design_bridge(design_class, gap)
  restrict_frame_shift <- if (design_class == "cds") {
    input$parameters$cds_fs
  } else {
    input$parameters$ncrna_fs
  }
  frame_status <- if (restrict_frame_shift) {
    "divisible_by_three"
  } else {
    "not_restricted"
  }
  if (design_class == "cds") {
    if (!restrict_frame_shift && gap %% 3 != 0) {
      frame_status <- sprintf("bridge shortened to %d nt", nchar(bridge))
    }
  }
  bridge_rc <- as.character(complement(reverse(DNAString(bridge))))
  final_arms <- DNAStringSet(c(
    left_homology_arm = arms$left[
      pair$PRIMER_LEFT_pos[[1]]:pair$PRIMER_RIGHT_pos[[1]]
    ],
    right_homology_arm = arms$right[
      pair$PRIMER_LEFT_pos[[2]]:pair$PRIMER_RIGHT_pos[[2]]
    ]
  ))
  writeXStringSet(final_arms, file.path(target_dir, "homology_arms.fasta"))
  sgrnas <- DNAStringSet(paste0(
    "ACG",
    input$parameters$site1,
    substr(selected$table$target_sequence, 1, 20),
    "GTTTTAGAGCTAGAAATAGCAAGTTaaaataaggct"
  ))
  names(sgrnas) <- paste0(feature$display_name, "_sgF", seq_along(sgrnas))
  arm_primers <- DNAStringSet(c(
    homology_primer_pair(pair[1, , drop = FALSE], "left", bridge, input$parameters$site2),
    homology_primer_pair(pair[2, , drop = FALSE], "right", bridge, input$parameters$site2)
  ))
  names(arm_primers) <- paste0(
    feature$display_name,
    c("_LF", "_LR", "_RF", "_RR")
  )
  sgrna_reverse <- DNAStringSet("AGTTGACGCTAAAAAAAGCACCGACTCGGTGCC")
  names(sgrna_reverse) <- paste0(feature$display_name, "_sgR")
  all_primers <- c(sgrnas, sgrna_reverse, arm_primers)
  writeXStringSet(all_primers, file.path(target_dir, "all_primers.fasta"))
  plain <- DNAStringSet(c(
    substr(selected$table$target_sequence, 1, 20),
    pair$PRIMER_LEFT_SEQUENCE[[1]],
    pair$PRIMER_RIGHT_SEQUENCE[[1]],
    pair$PRIMER_LEFT_SEQUENCE[[2]],
    pair$PRIMER_RIGHT_SEQUENCE[[2]]
  ))
  names(plain) <- c(names(sgrnas), names(arm_primers))
  plain_path <- file.path(target_dir, "primers_without_service_sequences.fasta")
  writeXStringSet(plain, plain_path)

  offtarget_range <- range(unlist(pair[, c("genome_start", "genome_end")]))
  screening_range <- pmax(
    1L,
    pmin(as.integer(offtarget_range + c(-200L, 200L)), length(input$genome))
  )
  screening_seq <- input$genome[screening_range[[1]]:screening_range[[2]]]
  screening <- callPrimer3(
    as.character(screening_seq),
    paste0(length(screening_seq) - 100L, "-", length(screening_seq)),
    c(62.5, 63, 63.5),
    2,
    "genome_screening",
    primer_num = 5,
    primer3 = input$tools$primer3,
    thermo.param = input$tools$primer3_config,
    settings = file.path(target_dir, "primer3_settings.txt"),
    report = file.path(target_dir, "genome_screening_report.txt")
  )
  if (!is.data.frame(screening) || !nrow(screening)) {
    stop("Не удалось подобрать скрининговые праймеры", call. = FALSE)
  }
  screening <- mutate(
    screening,
    genome_start = PRIMER_LEFT_pos + screening_range[[1]] - 1L,
    genome_end = PRIMER_RIGHT_pos + screening_range[[1]] - 1L
  )
  write_tsv(screening, file.path(target_dir, "screening_primer3_table.tsv"))
  screening_ranking <- list()
  for (screening_index in seq_len(nrow(screening))) {
    screening_row <- screening[screening_index, , drop = FALSE]
    pair_id <- sprintf("%s_screening_%02d", output_id, screening_index)
    if (!is.null(log_path)) {
      append_design_log(
        log_path,
        "primer_qc",
        "TRY",
        sprintf("pair_id=%s;reaction=scrF_scrR", pair_id)
      )
    }
    result <- evaluate_candidate_reaction(
      input,
      screening_row,
      screening_row$PRIMER_LEFT_SEQUENCE[[1]],
      screening_row$PRIMER_RIGHT_SEQUENCE[[1]],
      "scrF_scrR",
      pair_id,
      arms$primer_qc_trace
    )
    specificity <- result$specificity
    openprimer <- result$openprimer
    primer3_penalty <- if ("PRIMER_PAIR_PENALTY" %in% names(screening_row)) {
      suppressWarnings(as.numeric(screening_row$PRIMER_PAIR_PENALTY[[1]]))
    } else Inf
    screening_ranking[[length(screening_ranking) + 1L]] <- data.frame(
      pair_id = pair_id,
      reaction = "scrF_scrR",
      primer3_index = screening_index,
      filtering_level = input$parameters$filtering_level,
      filtering_mode = filtering_level_name(input$parameters$filtering_level),
      structure_passed = TRUE,
      specificity_passed = specificity$passed,
      openprimer_passed = !is.null(openprimer) && openprimer$passed,
      strict_qc_passed = result$policy$strict_passed,
      blocking_passed = result$policy$blocking_passed,
      n_expected_product_deviations = abs(
        specificity$n_expected_products - 1L
      ),
      n_high_risk_offtarget_products =
        specificity$n_high_risk_offtarget_products,
      n_all_offtarget_products = specificity$n_all_offtarget_products,
      n_perfect_3p_offtarget_sites =
        specificity$n_perfect_3p_offtarget_sites,
      openprimer_failed_soft_constraints = if (is.null(openprimer)) {
        0L
      } else openprimer$failed_soft_constraints,
      openprimer_penalty = if (is.null(openprimer)) Inf else openprimer$penalty,
      max_dimer_risk = if (is.null(openprimer)) Inf else {
        openprimer$max_dimer_risk
      },
      abs_tm_diff = if (is.null(openprimer)) Inf else openprimer$abs_tm_diff,
      primer3_pair_penalty = primer3_penalty,
      deleted_nt = gap,
      risk_warnings = if (nzchar(result$risk_warnings)) {
        paste0("scrF_scrR:", result$risk_warnings)
      } else "",
      rejection_reason = result$rejection_reason,
      stringsAsFactors = FALSE
    )
  }
  screening_selection <- select_best_primer_pair(bind_rows(screening_ranking))
  append_primer_qc_trace(
    arms$primer_qc_trace,
    "scrF_scrR",
    "screening_candidates",
    ranking = screening_selection$ranking
  )
  if (!is.null(log_path)) {
    for (i in seq_len(nrow(screening_selection$ranking))) {
      candidate <- screening_selection$ranking[i, , drop = FALSE]
      append_design_log(
        log_path,
        "primer_qc",
        if (candidate$selected[[1]]) "OK" else "REJECTED",
        sprintf(
          "pair_id=%s;reason=%s",
          candidate$pair_id[[1]],
          if (candidate$selected[[1]] &&
              isTRUE(candidate$fallback_selected[[1]])) {
            candidate$risk_warnings[[1]]
          } else candidate$rejection_reason[[1]]
        )
      )
    }
  }
  if (is.null(screening_selection$pair)) {
    stop("Все screening-пары отклонены primer QC", call. = FALSE)
  }
  screening <- screening[
    screening_selection$pair$primer3_index[[1]],
    ,
    drop = FALSE
  ]
  selected_screening_pair_id <- screening_selection$pair$pair_id[[1]]
  screening_primers <- DNAStringSet(c(
    screening$PRIMER_LEFT_SEQUENCE[[1]],
    screening$PRIMER_RIGHT_SEQUENCE[[1]]
  ))
  names(screening_primers) <- paste0(feature$display_name, "_scr", c("F", "R"))
  screening_path <- file.path(target_dir, "screening_primers.fasta")
  writeXStringSet(screening_primers, screening_path)

  # ticks are sorted genomic coordinates of the complete selected arms.
  # Preserve both arms, including primer sites, in the original FASTA orientation.
  edited_genome <- DNAStringSet(c(
    edited_genome = c(
      input$genome[1:arms$ticks[[2]]],
      DNAString(if (feature$strand == "-") bridge_rc else bridge),
      input$genome[arms$ticks[[3]]:length(input$genome)]
    )
  ))
  writeXStringSet(edited_genome, file.path(target_dir, "edited_genome.fasta"))
  screening_product_sizes <- calculate_screening_product_sizes(
    screening,
    input$genome,
    edited_genome
  )
  pcr_products <- model_design_pcr_products(
    input,
    feature,
    pair,
    sgrnas,
    sgrna_reverse,
    arm_primers,
    screening,
    screening_primers,
    edited_genome,
    screening_product_sizes
  )
  pcr_product_sequences <- DNAStringSet(pcr_products$sequence)
  names(pcr_product_sequences) <- pcr_products$name
  writeXStringSet(
    pcr_product_sequences,
    file.path(target_dir, "pcr_products.fasta")
  )
  write_tsv(pcr_products, file.path(target_dir, "pcr_products.tsv"))
  edited_ptargets <- model_edited_ptargets(
    input$target_plasmid_sequence,
    pcr_products$sequence[startsWith(pcr_products$name, "sgRNA_N20_")],
    pcr_products$sequence[pcr_products$name == "left_homology_arm"],
    bridge,
    pcr_products$sequence[pcr_products$name == "right_homology_arm"],
    input$parameters$site1,
    input$parameters$site2,
    feature$display_name
  )
  writeXStringSet(
    edited_ptargets$sequences,
    file.path(target_dir, "edited_pTargets.fasta")
  )
  final_sequences <- c(all_primers, screening_primers)
  primer_purposes <- c(
    "left_arm_forward_primer",
    "left_arm_reverse_primer",
    "right_arm_forward_primer",
    "right_arm_reverse_primer",
    "screening_forward_primer",
    "screening_reverse_primer"
  )
  sequence_purposes <- c(
    rep("sgRNA_forward_oligo", length(sgrnas)),
    "sgRNA_reverse_oligo",
    primer_purposes
  )
  primer_metrics <- data.frame(
    name = c(names(arm_primers), names(screening_primers)),
    purpose = primer_purposes,
    annealing_sequence = c(
      pair$PRIMER_LEFT_SEQUENCE[[1]],
      pair$PRIMER_RIGHT_SEQUENCE[[1]],
      pair$PRIMER_LEFT_SEQUENCE[[2]],
      pair$PRIMER_RIGHT_SEQUENCE[[2]],
      screening$PRIMER_LEFT_SEQUENCE[[1]],
      screening$PRIMER_RIGHT_SEQUENCE[[1]]
    ),
    tm_c = c(
      pair$PRIMER_LEFT_TM[[1]],
      pair$PRIMER_RIGHT_TM[[1]],
      pair$PRIMER_LEFT_TM[[2]],
      pair$PRIMER_RIGHT_TM[[2]],
      screening$PRIMER_LEFT_TM[[1]],
      screening$PRIMER_RIGHT_TM[[1]]
    ),
    stringsAsFactors = FALSE
  )
  selected_ranking <- bind_rows(arms$primer_qc_trace$ranking)
  selected_ranking <- selected_ranking[selected_ranking$selected, , drop = FALSE]
  homology_rank <- selected_ranking[
    selected_ranking$pair_id == arms$selected_pair_id,
    ,
    drop = FALSE
  ]
  screening_rank <- selected_ranking[
    selected_ranking$pair_id == selected_screening_pair_id,
    ,
    drop = FALSE
  ]
  selected_openprimer <- bind_rows(Filter(
    function(metrics) {
      "pair_id" %in% names(metrics) &&
        any(metrics$pair_id == selected_screening_pair_id)
    },
    arms$primer_qc_trace$openprimer
  ))
  selected_openprimer <- selected_openprimer[
    selected_openprimer$pair_id == selected_screening_pair_id,
    ,
    drop = FALSE
  ]
  if (
    nrow(homology_rank) != 1L ||
      nrow(screening_rank) != 1L ||
      nrow(selected_openprimer) != 1L
  ) {
    stop("Не удалось собрать QC выбранной screening-пары", call. = FALSE)
  }
  qc_warnings <- unique(trimws(c(
    homology_rank$risk_warnings[[1]],
    screening_rank$risk_warnings[[1]]
  )))
  qc_warnings <- qc_warnings[nzchar(qc_warnings)]
  fallback_used <- any(c(
    homology_rank$fallback_selected[[1]],
    screening_rank$fallback_selected[[1]]
  ))
  warning_text <- if (length(qc_warnings)) {
    paste(qc_warnings, collapse = "; ")
  } else "none"
  has_qc_warnings <- length(qc_warnings) > 0L
  qc_status <- if (has_qc_warnings) "selected_with_warnings" else "strict_pass"
  if (has_qc_warnings) {
    detail <- sprintf(
      "filtering_level=%d;target=%s;warnings=%s",
      input$parameters$filtering_level,
      feature$query_name,
      warning_text
    )
    if (!is.null(log_path)) {
      append_design_log(log_path, "primer_qc", "WARNING", detail)
    }
    warning(
      sprintf(
        "[%s] Выбраны праймеры с QC-рисками (%s): %s",
        feature$query_name,
        filtering_level_name(input$parameters$filtering_level),
        warning_text
      ),
      call. = FALSE,
      immediate. = TRUE
    )
  }
  n20_distances <- calculate_n20_arm_distances(
    selected$table,
    arms$ticks,
    feature$strand
  )
  screening_qc <- list(
    pair_id = selected_screening_pair_id,
    filtering_level = input$parameters$filtering_level,
    filtering_mode = filtering_level_name(input$parameters$filtering_level),
    selection_status = qc_status,
    fallback_used = fallback_used,
    warnings = warning_text,
    offtarget_products = screening_rank$n_all_offtarget_products[[1]],
    high_risk_offtarget_products =
      screening_rank$n_high_risk_offtarget_products[[1]],
    perfect_3p_offtarget_sites =
      screening_rank$n_perfect_3p_offtarget_sites[[1]],
    openprimer_metrics = format_openprimer_report_metrics(selected_openprimer)
  )

  writeLines(
    c(
      paste("target", feature$query_name, sep = "\t"),
      paste("class", design_class, sep = "\t"),
      paste("filtering_level", input$parameters$filtering_level, sep = "\t"),
      paste(
        "filtering_mode",
        filtering_level_name(input$parameters$filtering_level),
        sep = "\t"
      ),
      paste("primer_qc_selection_status", qc_status, sep = "\t"),
      paste("primer_qc_fallback_used", fallback_used, sep = "\t"),
      paste("primer_qc_warnings", warning_text, sep = "\t"),
      paste("n20_count", nrow(selected$table), sep = "\t"),
      paste("ptarget_site1", input$parameters$site1, sep = "\t"),
      paste("ptarget_site2", input$parameters$site2, sep = "\t"),
      paste(
        "ptarget_site_pair_orientation",
        edited_ptargets$restriction_pair$orientation,
        sep = "\t"
      ),
      paste(
        "edited_ptarget_count",
        length(edited_ptargets$sequences),
        sep = "\t"
      ),
      paste(
        "n20_strands",
        paste(sort(unique(selected$table$strand)), collapse = ","),
        sep = "\t"
      ),
      paste(
        "n20_offtarget_thresholds",
        paste(input$parameters$n20_offtarget, collapse = ","),
        sep = "\t"
      ),
      paste(
        "n20_arm_min_distance",
        input$parameters$n20_arm_min_distance,
        sep = "\t"
      ),
      paste("deleted_nt", gap, sep = "\t"),
      paste("frame_status", frame_status, sep = "\t"),
      paste("left_arm_nt", length(final_arms[[1]]), sep = "\t"),
      paste("right_arm_nt", length(final_arms[[2]]), sep = "\t"),
      paste(
        "left_forward_primer_tm",
        round(pair$PRIMER_LEFT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "left_reverse_primer_tm",
        round(pair$PRIMER_RIGHT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "right_forward_primer_tm",
        round(pair$PRIMER_LEFT_TM[[2]], 1),
        sep = "\t"
      ),
      paste(
        "right_reverse_primer_tm",
        round(pair$PRIMER_RIGHT_TM[[2]], 1),
        sep = "\t"
      ),
      paste(
        "screening_forward_tm",
        round(screening$PRIMER_LEFT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "screening_reverse_tm",
        round(screening$PRIMER_RIGHT_TM[[1]], 1),
        sep = "\t"
      ),
      paste(
        "screening_unsuccessful_insertion_bp",
        screening_product_sizes[["unsuccessful_insertion_bp"]],
        sep = "\t"
      ),
      paste(
        "screening_successful_insertion_bp",
        screening_product_sizes[["successful_insertion_bp"]],
        sep = "\t"
      ),
      paste("homology_pair_id", arms$selected_pair_id, sep = "\t"),
      paste("screening_pair_id", selected_screening_pair_id, sep = "\t"),
      paste(
        "primer_qc_homology_offtarget_products",
        homology_rank$n_all_offtarget_products[[1]],
        sep = "\t"
      ),
      paste(
        "primer_qc_screening_offtarget_products",
        screening_rank$n_all_offtarget_products[[1]],
        sep = "\t"
      ),
      paste(
        "primer_qc_homology_openprimer_penalty",
        homology_rank$openprimer_penalty[[1]],
        sep = "\t"
      ),
      paste(
        "primer_qc_screening_openprimer_penalty",
        screening_rank$openprimer_penalty[[1]],
        sep = "\t"
      )
    ),
    file.path(target_dir, "report.tsv")
  )
  completed <- TRUE
  list(
    all_primers_path = file.path(target_dir, "all_primers.fasta"),
    plain_path = plain_path,
    screening_path = screening_path,
    wet_lab = list(
      sequences = final_sequences,
      sequence_purposes = sequence_purposes,
      primer_metrics = primer_metrics,
      screening_product_sizes = screening_product_sizes,
      n20_distances = n20_distances,
      screening_qc = screening_qc,
      edited_genome = edited_genome,
      edited_ptargets = edited_ptargets$sequences,
      ptarget_site_pair = edited_ptargets$restriction_pair,
      pcr_products = pcr_products
    ),
    primer_qc_trace = arms$primer_qc_trace,
    homology_pair_id = arms$selected_pair_id,
    screening_pair_id = selected_screening_pair_id
  )
}

append_design_log <- function(log_path, stage, status, detail = "") {
  line <- paste(
    format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    stage,
    status,
    detail,
    sep = "\t"
  )
  cat(line, "\n", file = log_path, append = TRUE, sep = "")
}

run_design_stage <- function(stage, log_path, expression) {
  append_design_log(log_path, stage, "START")
  tryCatch(
    {
      value <- force(expression)
      append_design_log(log_path, stage, "OK")
      value
    },
    error = function(e) {
      reason <- conditionMessage(e)
      error_stage <- if (!is.null(e$stage)) e$stage else stage
      append_design_log(log_path, error_stage, "ERROR", reason)
      stop(structure(
        list(
          message = reason,
          call = NULL,
          stage = error_stage,
          parent = e
        ),
        class = c("design_stage_error", "error", "condition")
      ))
    }
  )
}

design_from_grna_pool <- function(
  input,
  feature,
  grnas,
  design_class,
  target_dir,
  log_path
) {
  attempted_ranges <- new.env(parent = emptyenv())
  primer_qc_trace <- new_primer_qc_trace()
  attempts <- 0L
  last_failure_stage <- "homology_arms"
  result <- visit_grna_sets(
    grnas,
    input$parameters$n20_mn,
    input$parameters$n20_strands,
    function(candidates) {
      selected <- list(
        table = candidates,
        n20_range = range(c(candidates$n20_start, candidates$n20_end))
      )
      range_key <- paste(selected$n20_range, collapse = ":")
      if (exists(range_key, envir = attempted_ranges, inherits = FALSE)) {
        return(NULL)
      }
      assign(range_key, TRUE, envir = attempted_ranges)
      attempts <<- attempts + 1L
      detail <- sprintf(
        "set=%d;n20=%s;range=%s",
        attempts,
        paste(candidates$genomic_location, collapse = ","),
        range_key
      )
      append_design_log(log_path, "homology_arms", "TRY", detail)
      arms <- design_homology_arms(
        input,
        feature,
        selected,
        design_class,
        target_dir,
        primer_qc_trace,
        log_path,
        attempts
      )
      if (is.null(arms)) {
        last_failure_stage <<- if (length(primer_qc_trace$ranking)) {
          "primer_qc"
        } else {
          "homology_arms"
        }
        append_design_log(
          log_path,
          "homology_arms",
          "REJECTED",
          sprintf("set=%d", attempts)
        )
        return(NULL)
      }
      append_design_log(
        log_path,
        "design_outputs",
        "TRY",
        sprintf("set=%d", attempts)
      )
      output_attempt <- tryCatch(
        list(
          ok = TRUE,
          value = write_design_outputs(
            input,
            feature,
            selected,
            arms,
            design_class,
            target_dir,
            log_path
          )
        ),
        error = function(e) list(ok = FALSE, error = e)
      )
      if (!output_attempt$ok) {
        if (inherits(output_attempt$error, "primer3_error")) stop(output_attempt$error)
        last_failure_stage <<- if (grepl(
          "primer QC|openPrimeR|constraint|specificity",
          conditionMessage(output_attempt$error),
          ignore.case = TRUE
        )) {
          "primer_qc"
        } else {
          "design_outputs"
        }
        append_design_log(
          log_path,
          "design_outputs",
          "REJECTED",
          sprintf(
            "set=%d;reason=%s",
            attempts,
            conditionMessage(output_attempt$error)
          )
        )
        return(NULL)
      }
      append_design_log(
        log_path,
        "design_outputs",
        "OK",
        sprintf("set=%d", attempts)
      )
      write_primer_qc_trace(primer_qc_trace, target_dir)
      list(
        selected = selected,
        arms = arms,
        primer_paths = output_attempt$value,
        attempts = attempts
      )
    }
  )
  if (is.null(result)) {
    write_primer_qc_trace(primer_qc_trace, target_dir)
    reason <- sprintf(
      paste(
        "Не удалось завершить дизайн ни для одного допустимого набора N20",
        "(проверено уникальных диапазонов: %d)"
      ),
      attempts
    )
    stop(structure(
      list(
        message = reason,
        call = NULL,
        stage = last_failure_stage
      ),
      class = c("grna_sets_exhausted", "error", "condition")
    ))
  }
  result
}

write_target_error <- function(
  path,
  gene_name,
  design_class,
  stage,
  reason
) {
  writeLines(
    c(
      paste("timestamp", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), sep = "\t"),
      paste("gene", gene_name, sep = "\t"),
      paste("class", design_class, sep = "\t"),
      paste("stage", stage, sep = "\t"),
      paste("reason", reason, sep = "\t")
    ),
    path
  )
}

design_target <- function(input, genome_name, gene_name, design_class) {
  safe_name <- tolower(gsub("[^A-Za-z0-9_.-]", "_", gene_name))
  layout <- output_layout(input$output_dir)
  target_dir <- file.path(
    layout$tech_report,
    paste0(safe_name, "_results")
  )
  wet_lab_dir <- file.path(
    layout$wet_lab,
    paste0(safe_name, "_results")
  )
  unlink(file.path(
    wet_lab_dir,
    c(
      "final_sequences.fasta",
      "final_sequences.txt",
      "wet_lab_report.txt",
      "edited_genome.fasta",
      "edited_pTargets.fasta",
      "pcr_products.fasta",
      "pcr_products.tsv"
    )
  ))
  dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
  log_path <- file.path(target_dir, "design.log")
  error_path <- file.path(target_dir, "error.txt")
  writeLines(
    paste(
      format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
      "target",
      "START",
      sprintf("gene=%s;class=%s", gene_name, design_class),
      sep = "\t"
    ),
    log_path
  )
  if (file.exists(error_path)) {
    unlink(error_path)
  }
  message(sprintf("[%s] %s", design_class, gene_name))
  tryCatch(
    {
      feature <- run_design_stage(
        "feature_lookup",
        log_path,
        feature_record(input, gene_name)
      )
      run_design_stage(
        "chopchop",
        log_path,
        run_chopchop(
          input,
          genome_name,
          feature,
          design_class,
          target_dir
        )
      )
      grnas <- run_design_stage(
        "n20_filter",
        log_path,
        filter_grnas(
          file.path(target_dir, "n20_table.tsv"),
          feature,
          design_class,
          input$parameters$n20_offtarget,
          input$genome
        )
      )
      grnas <- run_design_stage(
        "n20_pool",
        log_path,
        prepare_grna_pool(
          grnas,
          input$parameters$n20_mn,
          input$parameters$n20_strands
        )
      )
      design <- run_design_stage(
        "homology_arms",
        log_path,
        design_from_grna_pool(
          input,
          feature,
          grnas,
          design_class,
          target_dir,
          log_path
        )
      )
      run_design_stage(
        "n20_output",
        log_path,
        write_selected_grnas(design$selected, feature, target_dir)
      )
      run_design_stage(
        "wet_lab_output",
        log_path,
        write_wet_lab_outputs(
          wet_lab_dir,
          feature,
          design_class,
          design$primer_paths$wet_lab$sequences,
          design$primer_paths$wet_lab$sequence_purposes,
          design$primer_paths$wet_lab$primer_metrics,
          design$primer_paths$wet_lab$screening_product_sizes,
          design$primer_paths$wet_lab$n20_distances,
          design$primer_paths$wet_lab$screening_qc,
          design$primer_paths$wet_lab$edited_genome,
          design$primer_paths$wet_lab$edited_ptargets,
          design$primer_paths$wet_lab$ptarget_site_pair,
          design$primer_paths$wet_lab$pcr_products
        )
      )
      append_design_log(log_path, "target", "OK")
      data.frame(
        gene = gene_name,
        class = design_class,
        output_dir = target_dir,
        status = "ok",
        stage = NA_character_,
        reason = NA_character_,
        wet_lab_dir = wet_lab_dir
      )
    },
    error = function(e) {
      stage <- if (inherits(e, "design_stage_error")) {
        e$stage
      } else {
        "unknown"
      }
      reason <- conditionMessage(e)
      write_target_error(
        error_path,
        gene_name,
        design_class,
        stage,
        reason
      )
      append_design_log(log_path, "target", "ERROR", paste(stage, reason))
      stop(structure(
        list(
          message = reason,
          call = NULL,
          gene = gene_name,
          design_class = design_class,
          stage = stage,
          target_dir = target_dir
        ),
        class = c("target_design_error", "error", "condition")
      ))
    }
  )
}

main <- function(args = commandArgs(trailingOnly = TRUE)) {
  cli <- parse_designer_args(args)
  configure_openprimer_environment()
  input <- make_design_input(cli)
  message(sprintf(
    "[QC] filtering_level=%d (%s)",
    input$parameters$filtering_level,
    filtering_level_name(input$parameters$filtering_level)
  ))
  loaded_openprimer <- tryCatch(
    load_openprimer_settings(
      input$parameters$primer_qc,
      input$parameters$primer3_buffer
    ),
    error = identity
  )
  if (inherits(loaded_openprimer, "error")) {
    assign(
      ".openprimer_load_error",
      conditionMessage(loaded_openprimer),
      input$primer_qc_cache
    )
    warning(
      sprintf("openPrimeR QC недоступен: %s", conditionMessage(loaded_openprimer)),
      call. = FALSE,
      immediate. = TRUE
    )
  } else {
    assign(".openprimer_settings", loaded_openprimer, input$primer_qc_cache)
  }
  dir.create(input$output_dir, recursive = TRUE, showWarnings = FALSE)
  layout <- output_layout(input$output_dir)
  dir.create(layout$wet_lab, recursive = TRUE, showWarnings = FALSE)
  dir.create(layout$tech_report, recursive = TRUE, showWarnings = FALSE)
  targets <- bind_rows(
    data.frame(gene = cli$cds, class = rep("cds", length(cli$cds))),
    data.frame(gene = cli$ncrna, class = rep("ncrna", length(cli$ncrna)))
  )
  write_run_parameters(
    input,
    targets,
    file.path(layout$tech_report, "run_parameters.tsv")
  )
  genome_assets <- prepare_chopchop_assets(input)
  chopchop_config <- configure_chopchop(input, genome_assets)
  if (!file.copy(
    chopchop_config,
    file.path(layout$tech_report, "chopchop_config.json"),
    overwrite = TRUE
  )) {
    stop("Не удалось сохранить конфигурацию CHOPCHOP в TechReport", call. = FALSE)
  }
  results <- lapply(seq_len(nrow(targets)), function(i) {
    tryCatch(
      design_target(
        input,
        genome_assets$name,
        targets$gene[[i]],
        targets$class[[i]]
      ),
      error = function(e) {
        stage <- if (inherits(e, "target_design_error")) {
          e$stage
        } else {
          "unknown"
        }
        output_dir <- if (inherits(e, "target_design_error")) {
          e$target_dir
        } else {
          NA_character_
        }
        message(sprintf(
          "[ERROR] gene=%s class=%s stage=%s: %s",
          targets$gene[[i]],
          targets$class[[i]],
          stage,
          conditionMessage(e)
        ))
        data.frame(
          gene = targets$gene[[i]],
          class = targets$class[[i]],
          output_dir = output_dir,
          status = "error",
          stage = stage,
          reason = conditionMessage(e),
          wet_lab_dir = NA_character_
        )
      }
    )
  })
  write_tsv(
    bind_rows(results),
    file.path(layout$tech_report, "design_summary.tsv")
  )
}

if (sys.nframe() == 0) {
  main()
}
