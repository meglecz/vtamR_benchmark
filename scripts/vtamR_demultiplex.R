args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
vsearch_path <- as.character(args[2])
cutadapt_path <- as.character(args[3])
pigz_path <- as.character(args[4])
outdir <- as.character(args[5])

library(vtamR)

# setwd("~/vtamR_benchmark")

fastainfo <- file.path(outdir, dataset, "vtamR", "merge", "fastainfo.csv")
merged_dir <- file.path(outdir, dataset, "vtamR", "merge")
demultiplexed_dir <- file.path(outdir, dataset, "vtamR", "demultiplex_trim")

sampleinfo_df <- demultiplex_and_trim(
  fastainfo = fastainfo, 
  fasta_dir=merged_dir, 
  outdir=demultiplexed_dir, 
  check_reverse=TRUE, # to match vtam
  tag_to_end = FALSE, # to match vtam
  primer_to_end = FALSE,  # to match vtam
  cutadapt_path=cutadapt_path,
  vsearch_path = vsearch_path,
  compress_method = "pigz",
  pigz_path = pigz_path,
  compress=TRUE
)
