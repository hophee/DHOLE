#!/usr/bin/env Rscript
source("oligo_designer.R")

local({
  check <- function(value, message) {
    if (!isTRUE(value)) stop(message, call. = FALSE)
  }
  expect_error <- function(args, pattern) {
    error <- tryCatch(parse_designer_args(args), error = identity)
    check(inherits(error, "error") && grepl(pattern, conditionMessage(error)),
          paste("Expected argument error:", pattern))
  }
  common <- c("--target-plasmid", "target.fa", "--output-dir", "output",
              "--gene", "shared,only_first,missing")
  args <- c("--genome", "first.fa", "second.fa",
            "--genome-annotation", "first.gff", "second.tsv",
            "--annotation-format", "gff", "bakta", common)
  parsed <- parse_designer_args(args)
  check(identical(parsed$genome, c("first.fa", "second.fa")) &&
          identical(parsed$genome_annotation, c("first.gff", "second.tsv")) &&
          identical(parsed$annotation_format, c("gff", "bakta")),
        "Genome/annotation/format order must be preserved")
  comma <- parse_designer_args(c(
    "--genome", "first.fa,second.fa", "--genome-annotation", "same.gff,same.gff", common
  ))
  check(identical(comma$genome, parsed$genome) &&
          identical(comma$genome_annotation, c("same.gff", "same.gff")) &&
          identical(comma$annotation_format, "gff"),
        "Comma lists must preserve repeated annotation paths and the shared format")
  expect_error(c("--genome", "first.fa", "second.fa",
                 "--genome-annotation", "first.gff", common), "one annotation per genome")
  expect_error(c("--genome", "first.fa", "--genome-annotation", "first.gff",
                 "--annotation-format", "gff", "bakta", common), "one format")
  expect_error(c("--genome", "one/strain.fa", "two/strain.fna",
                 "--genome-annotation", "a.gff", "b.gff", common), "distinct")
  expect_error(c("--genome", "design_summary.tsv.fa", "second.fa",
                 "--genome-annotation", "a.gff", "b.gff", common), "conflict with root reports")
  expect_error(c("--genome", "first.fa,,second.fa", "--genome-annotation", "a.gff", "b.gff", common),
               "empty entries")
  expect_error(c("--genome", "first.fa", "--genome-annotation", "a.gff",
                 "--annotation-format", "unknown", common), "Invalid --annotation-format")

  directory <- tempfile("multi-genome-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  previous <- file.path(directory, "previous.tsv")
  write_tsv(data.frame(sequence = "AAAA", name = "old_primer", gene = "shared",
                      genome = "old", source_gene = "shared"), previous)
  output <- file.path(directory, "batch")
  args[args == "output"] <- output
  args <- c(args, "--primer-table", previous)

  # Isolate the expensive per-genome design, retaining real orchestration and merging.
  batch <- main
  environment(batch) <- new.env(parent = environment(main))
  assign("run_genome", function(cli, previous_primers) {
    check(identical(cli$gene, c("shared", "only_first", "missing")),
          "Every genome must receive the complete requested gene set")
    expected_annotation <- if (cli$genome == "first.fa") "first.gff" else "second.tsv"
    expected_format <- if (cli$genome == "first.fa") "gff" else "bakta"
    check(cli$genome_annotation == expected_annotation &&
            cli$annotation_format == expected_format, "Annotations paired with wrong genomes")
    wet_lab <- file.path(cli$output_dir, "WetLab", "shared_results")
    dir.create(wet_lab, recursive = TRUE)
    write_tsv(data.frame(name = "shared_LF", sequence = "ACGT"),
              file.path(wet_lab, "final_sequences.txt"))
    data.frame(gene = cli$gene, class = "cds",
               status = c("ok", "error", "missing"),
               wet_lab_dir = c(wet_lab, NA, NA))
  }, environment(batch))
  batch(args)
  summary <- read.delim(file.path(output, "design_summary.tsv"))
  merged <- read.delim(file.path(output, "unique_primers.tsv"))
  check(nrow(summary) == 6L &&
          identical(summary$genome, rep(c("first", "second"), each = 3L)) &&
          identical(summary$status, rep(c("ok", "error", "missing"), 2L)),
        "The root summary must retain every genome's successes, errors and missing genes")
  check(nrow(merged) == 2L && merged$gene[merged$sequence == "AAAA"] == "shared_old" &&
          merged$gene[merged$sequence == "ACGT"] == "shared_first;shared_second" &&
          merged$genome[merged$sequence == "ACGT"] == "first;second",
        "Root collection must merge all new genomes and the imported collection")
  check(file.exists(file.path(output, "first", "WetLab", "shared_results", "final_sequences.txt")) &&
          file.exists(file.path(output, "second", "WetLab", "shared_results", "final_sequences.txt")),
        "Genome outputs must use separate directories named after FASTA files")

  # Real pipeline with no matching genes: no external design tools should be needed.
  first <- file.path(directory, "first.fa")
  second <- file.path(directory, "second.fa")
  gff <- file.path(directory, "first.gff")
  bakta <- file.path(directory, "second.tsv")
  plasmid <- file.path(directory, "target.fa")
  writeLines(c(">chr", strrep("ACGT", 100L)), first)
  file.copy(first, second)
  writeLines(c("##gff-version 3",
    "chr\ttest\tCDS\t10\t90\t.\t+\t0\tID=LOC1;locus_tag=LOC1;gene=other"), gff)
  write_tsv(data.frame(seqid = "chr", type = "CDS", start = 10L, stop = 90L,
                      strand = "+", locus_tag = "LOC2", gene = "other"), bakta)
  writeLines(c(">target", paste0("ACTAGT", strrep("C", 20L), SGRNA_SCAFFOLD,
                               "CTGCAG", strrep("A", 200L))), plasmid)
  output <- file.path(directory, "missing_batch")
  stale <- file.path(output, "first", "TechReport", "missing_results")
  stale_wet <- file.path(output, "first", "WetLab", "missing_results")
  dir.create(stale, recursive = TRUE)
  dir.create(stale_wet, recursive = TRUE)
  writeLines("old", file.path(stale, "report.tsv"))
  writeLines("old", file.path(stale_wet, "final_sequences.txt"))
  real_args <- c("--genome", first, second, "--genome-annotation", gff, bakta,
                 "--annotation-format", "gff", "bakta", "--target-plasmid", plasmid,
                 "--output-dir", output, "--gene", "missing", "also_missing",
                 "--primer-table", previous)
  messages <- character()
  withCallingHandlers(main(real_args), message = function(m) {
    messages <<- c(messages, conditionMessage(m))
    invokeRestart("muffleMessage")
  })
  summary <- read.delim(file.path(output, "design_summary.tsv"))
  check(nrow(summary) == 4L && all(summary$status == "missing") &&
          sum(grepl("\\[MISSING\\]", messages)) == 4L,
        "Every absent gene must be reported for every genome in the summary and console")
  check(identical(list.files(stale), "missing_gene.txt") && !dir.exists(stale_wet),
        "A missing gene must leave exactly one TXT file and no stale successful outputs")
  check(all(file.exists(file.path(output, c("first", "second"), "TechReport",
                                  "also_missing_results", "missing_gene.txt"))),
        "Missing targets need their own file in each genome directory")
  note <- readLines(file.path(stale, "missing_gene.txt"))
  check(grepl("missing", note) && grepl(gff, note, fixed = TRUE) &&
          grepl("first", note), "Missing-gene note must identify gene, genome and annotation")
  check(identical(read_primer_table(file.path(output, "unique_primers.tsv")),
                  read_primer_table(previous)),
        "An all-missing batch must preserve imported primers")

  single_output <- file.path(directory, "single")
  single_args <- c("--genome", first, "--genome-annotation", gff,
                   "--target-plasmid", plasmid, "--output-dir", single_output,
                   "--cds", "missing")
  main(single_args)
  check(file.exists(file.path(single_output, "TechReport", "design_summary.tsv")) &&
          file.exists(file.path(single_output, "WetLab", "unique_primers.tsv")) &&
          !file.exists(file.path(single_output, "design_summary.tsv")) &&
          !dir.exists(file.path(single_output, "first")) &&
          nrow(read_primer_table(file.path(single_output, "WetLab", "unique_primers.tsv"))) == 0L,
        "Single-genome output layout and readable empty primer table must remain unchanged")

  # A later annotation contains the target, but its deletion overlaps a neighbor.
  input <- make_design_input(parse_designer_args(single_args))
  input$annotation$gene <- "missing"
  input$annotation <- rbind(input$annotation,
    transform(input$annotation, gene = "neighbor", locus_tag = "LOC_NEIGHBOR"))
  error <- tryCatch(design_target(input, "first", "missing", "ncrna"), error = identity)
  target_dir <- file.path(single_output, "TechReport", "missing_results")
  check(inherits(error, "target_design_error") && error$stage == "deletion_bounds" &&
          !file.exists(file.path(target_dir, "missing_gene.txt")) &&
          file.exists(file.path(target_dir, "error.txt")),
        "A newly found target must clear its missing note and retain normal design error reporting")
})

message("Multi-genome tests passed")
