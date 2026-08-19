args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
vsearch_path <- as.character(args[2])
pigz_path <- as.character(args[3])
outdir <-  as.character(args[4])

library(vtamR)

# setwd("~/vtamR_benchmark")

fastqinfo <- paste0("vtamR_fastqinfo_", dataset, ".csv")
fastqinfo <- file.path("data", "metadata", fastqinfo)
if(dataset == "bat"){
  outdir = file.path("data", "fastq", "bat", "merge")
}else{
  outdir = file.path(outdir, dataset, "vtamR", "merge")
}
fastq_dir <- file.path("data", "fastq", dataset)
log_file <- file.path(outdir, dataset, "vtamR", "vtamR_log.csv")

fastainfo_df <- merge_fastq_pairs(
  fastqinfo, 
  fastq_dir=fastq_dir, 
  outdir=outdir,
  fastq_maxee=1,
  fastq_maxns=0,
  fastq_maxmergelen=500, # Adjusted to vtam
  fastq_minmergelen= 100, # Adjusted to vtam
  fastq_maxlen = 10000, # no default value in vtam
  vsearch_path = vsearch_path,
  compress=TRUE,
  compress_method="pigz",
  pigz_path=pigz_path,
  log_file = log_file
)

