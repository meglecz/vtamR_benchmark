args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
vsearch_path <- as.character(args[2])
cutadapt_path <- as.character(args[3])
pigz_path <- as.character(args[4])
log_file <- as.character(args[5])

library(vtamR)
library(dplyr)


# setwd("~/vtamR_benchmark")

fastainfo <- "data/fastq/bat/merge/fastainfo.csv"
merged_dir <- "data/fastq/bat/merge"
demultiplexed_dir <- "data/fastq/bat/demultiplex_trim"
  
## Trim primers
sampleinfo_df <- trim_primers(
  fastainfo = fastainfo, 
  fasta_dir=merged_dir, 
  outdir=demultiplexed_dir, 
  cutadapt_path=cutadapt_path, 
  vsearch_path=vsearch_path,
  primer_to_end=FALSE, # important !!!,
  compress=TRUE,
  compress_method="pigz",
  pigz_path=pigz_path,
  quiet = TRUE,
  log_file = log_file
)



sample_info <- read.csv(file.path(demultiplexed_dir, "sampleinfo.csv"))

sample_info <- sample_info %>%
  mutate(
    run = "bat",
    marker = "coi",
    ) %>%
  select(
    run, marker, sample, replicate, sortedfasta=fasta
  )

write.table(sample_info, file.path(demultiplexed_dir, "sortedinfo.tsv"), sep="\t", row.names = FALSE)



