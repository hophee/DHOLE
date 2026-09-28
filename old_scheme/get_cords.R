suppressPackageStartupMessages(library(Biostrings))

# --genome
# --genome_annotation
# --gene_name
# --annotation_format (optional: bakta, gff; default: bakta)

args <- commandArgs(trailingOnly = TRUE)
script_path <- sub("^--file=", "", commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))][1])
source(file.path(dirname(normalizePath(script_path)), "pipeline_utils.R"))
input <- parse_named_args(
  args,
  required = c("genome", "genome_annotation", "gene_name"),
  optional = "annotation_format"
)
annotation_format <- if (is.null(input$annotation_format)) "bakta" else input$annotation_format

genome_name <- readDNAStringSet(input$genome)@ranges@NAMES[1]
genome_table <- read_genome_annotation(input$genome_annotation, annotation_format)
g_idx <- find_target_feature(genome_table, input$gene_name)
gene_start <- genome_table$start[g_idx]
gene_end <- genome_table$stop[g_idx]

gene_length <- (gene_end - gene_start + 1)

# ranges_get <- function(x){
#   if (x<=500){
#     return(x)
#   } else if(x>500 & x< 1500){
#     return(500+((x-500)/10))
#   } else{
#     return(x/3)
#   }
# }
#
# plot(x=1:3000, y=sapply(1:3000, ranges_get))

if (gene_length <= 500) {
  #  gap <- round(gene_length/10, -1)
  cut_start <- gene_start #+ gap
  cut_end <- gene_end #- gap
} else if (gene_length > 500 & gene_length < 1500) {
  side <- 250
  mid <- gene_start + gene_length %/% 2
  cut_start <- mid - side
  cut_end <- mid + side
} else {
  cut_start <- floor((gene_end - gene_start) / 3) + gene_start
  cut_end <- gene_end - floor((gene_end - gene_start) / 3)
}
cat(paste0(genome_name, ":", cut_start, "-", cut_end))
