library(vtamR)

setwd("~/vtamR_benchmark")

## Concatenate files of the same name; 
# These are fastq files of the same samples, same amplicons run twice in Illumina
fastq_dirs = c("data/fastq/shark/200309", "data/fastq/shark/200513")
outdir <- "data/fastq/shark"
concatenate_files(
  dirs=fastq_dirs, 
  outdir=outdir, 
  pattern = "\\.fastq\\.gz"
  )
