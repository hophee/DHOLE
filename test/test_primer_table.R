#!/usr/bin/env Rscript
source("oligo_designer.R")

assert_true <- function(value, message) {
  if (!isTRUE(value)) stop(message, call. = FALSE)
}
assert_error <- function(expression, pattern) {
  error <- tryCatch({ force(expression); NULL }, error = identity)
  assert_true(inherits(error, "error") && grepl(pattern, conditionMessage(error)),
              paste("Expected error:", pattern))
}
base_args <- c("--genome", "genome.fasta", "--genome-annotation", "genome.gff",
               "--target-plasmid", "target.fasta", "--output-dir", "output", "--cds", "geneA")
defaults <- parse_designer_args(base_args)

local({
  directory <- tempfile("unique-primers-")
  dir.create(directory)
  on.exit(unlink(directory, recursive = TRUE))
  output <- file.path(directory, "unique_primers.tsv")
  fixture <- function(gene, sequences, names) {
    target <- file.path(directory, gene)
    dir.create(target, showWarnings = FALSE)
    write_tsv(data.frame(name = names, sequence = sequences),
              file.path(target, "final_sequences.txt"))
    data.frame(gene = gene, status = "ok", wet_lab_dir = target)
  }
  gene_a <- fixture("gene_A", c("acgt", "CCCC"), c("gene_A_sgR", "gene_A_LF"))
  gene_b <- fixture("geneB", c("ACGT", "GGGG"), c("geneB_sgR", "geneB_scrF"))
  failed <- data.frame(gene = "failed", status = "error", wet_lab_dir = "nonexistent")
  results <- bind_rows(gene_a, gene_b, failed)
  first <- write_unique_primers(results, "path/strain_one.fna", read_primer_table(), output)
  assert_true(nrow(first) == 3L && first$gene[[1]] == "gene_A;geneB" &&
    first$name[[1]] == "gene_A_sgR;geneB_sgR" &&
    first$genome[[1]] == "strain_one;strain_one",
    "Shared primers must have one row with both genes, names and paired genome sources")
  assert_true(parse_designer_args(c(base_args, "--primer-table", output))$primer_table == output &&
    !length(defaults$primer_table), "Previous primer table option was not parsed")
  previous <- read_primer_table(output)
  again <- write_unique_primers(results, "strain_one.fna", previous, output)
  assert_true(identical(first, again), "Reimporting the same run must not duplicate sources")

  gene_a_new <- fixture("gene_A", c("ACGT", "TTTT"), c("gene_A_sgR", "gene_A_RR"))
  merged <- write_unique_primers(gene_a_new, "strain_two.fasta", read_primer_table(output), output)
  assert_true(nrow(merged) == 4L &&
    merged$gene[[1]] == "gene_A_strain_one;geneB;gene_A_strain_two" &&
    merged$genome[[1]] == "strain_one;strain_one;strain_two" &&
    merged$gene[merged$sequence == "CCCC"] == "gene_A_strain_one" &&
    merged$gene[merged$sequence == "TTTT"] == "gene_A_strain_two" &&
    merged$source_gene[[1]] == "gene_A;geneB;gene_A",
    "Same gene in different genomes must be qualified across every primer, preserving raw names")
  roundtrip <- write_unique_primers(gene_a_new, "strain_two.fasta", read_primer_table(output), output)
  assert_true(identical(merged, roundtrip), "Qualified names with underscores must survive reimport")
  kept <- write_unique_primers(failed, "strain_three.fna", read_primer_table(output), output)
  assert_true(identical(merged, kept), "An unsuccessful run must preserve the previous collection")
  empty <- write_unique_primers(failed, "strain_three.fna", read_primer_table(), output)
  assert_true(nrow(empty) == 0L && nrow(read_primer_table(output)) == 0L,
              "No successful targets must still produce a readable header-only table")

  invalid <- merged
  invalid$genome[[1]] <- "strain_one"
  write_tsv(invalid, output)
  assert_error(read_primer_table(output), "inconsistent gene/genome")
  invalid <- merged
  invalid$sequence[[1]] <- "invalid!"
  write_tsv(invalid, output)
  assert_error(read_primer_table(output), "invalid DNA")
  write_tsv(merged[, "sequence", drop = FALSE], output)
  assert_error(read_primer_table(output), "must contain")
  writeLines(c("sequence\tname\tgene\tgenome\tsource_gene", "ACGT\tprimer\tgene"), output)
  suppressWarnings(assert_error(read_primer_table(output), "malformed TSV"))
})

message("Primer table tests passed")
