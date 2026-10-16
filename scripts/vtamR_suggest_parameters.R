args <- commandArgs(trailingOnly = TRUE)

dataset <- as.character(args[1])
vsearch_path  <- as.character(args[2])
by_replicate <- as.logical(args[3])
outdir_pipe <- as.character(args[4])
log_file <- as.character(args[5])


library(vtamR)
library(dplyr)

# setwd("~/vtamR_benchmark")

outdir = file.path(outdir_pipe, dataset, "vtamR","optimize_suggest")
mock_composition <- paste("vtamR_mock_composition_", dataset, ".csv", sep="")
mock_composition <- file.path("data", "metadata", mock_composition)

if(dataset == "bat"){
  sampleinfo <- file.path("data", "fastq", "bat", "demultiplex_trim", "sampleinfo.csv")
}else{
  sampleinfo <- file.path(outdir_pipe, dataset, "vtamR", "demultiplex_trim", "sampleinfo.csv")
}

# The input table in VTAM is identical in content to 2_LFNglobalReadCount.csv
# The final filtering with the suggested parameters is done on this file
read_count_raw_file <- file.path(outdir_pipe, dataset, "vtamR", "filter_default", "2_LFNglobalReadCount.csv")
read_count_raw <- read.csv(read_count_raw_file)
# To establish known_occurrences df used in the suggest functions, use the filtered dataset
read_count_prefilter_file <- file.path(outdir_pipe, dataset, "vtamR", "filter_default", "12_pool_replicates.csv")

run_time <-  data.frame(
  step = character(),
  elapsed = numeric(),
  user = numeric(),
  system = numeric()
)


############### suggest 
# The aim was to benchmark vtamR against VTAM. Therefore, to ensure that the two approaches were directly comparable:
#  - All optimize/suggest steps were run on the unfiltered dataset, following the VTAM approach.
#  - The known occurrences (TP, FP, and FN) from the control samples were calculated based on the prefiltered data.
#  - The filter_lfn_sample_replicate and filter_pcr_error steps were deactivated in VTAM when searching for lfn_variant_cutoff and lfn_read_count_cutoff.

# 1. Run classify_control_occurrences on read_count_prefilter_file to make known_occurrences

start_time <- proc.time()

known_occurrences = file.path(outdir, "known_occurrences.csv")
false_negatives = file.path(outdir, "false_negatives.csv")
performance_metrics = file.path(outdir, "performance_metrics.csv")

results <- classify_control_occurrences(
  read_count = read_count_prefilter_file,
  sampleinfo = sampleinfo,
  mock_composition = mock_composition,
  known_occurrences = known_occurrences,
  false_negatives = false_negatives,
  performance_metrics = performance_metrics,
  log_file = log_file
)

known_occurrences <- results[[1]]

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "classify_control_occurrences", elapsed = timing[3], user = timing[1], system =timing[2])

# 2. run suggest_xxx functions using the mock composition and known_occurrences file on input data (read_count_raw)

# 2.a suggest_sample_cutoff
start_time <- proc.time()
outfile = file.path(outdir, "suggest_sample_cutoff.csv")

suggest_sample_cutoff_df <- suggest_sample_cutoff(
  read_count=read_count_raw,
  mock_composition=mock_composition,
  outfile=outfile,
  log_file = log_file
)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "suggest_sample_cutoff", elapsed = timing[3], user = timing[1], system =timing[2])

# 2.b suggest_pcr_error_cutoff
start_time <- proc.time()
outfile = file.path(outdir, "suggest_pcr_error_cutoff.csv")

suggest_pcr_error_cutoff_df <- suggest_pcr_error_cutoff(
  read_count = read_count_raw, 
  mock_composition = mock_composition,
  vsearch_path = vsearch_path, 
  outfile = outfile,
  log_file = log_file
)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "suggest_pcr_error_cutoff", elapsed = timing[3], user = timing[1], system =timing[2])

# 2.c suggest_variant_read_count
start_time <- proc.time()

suggest_variant_readcount_cutoffs_df <- suggest_variant_readcount_cutoffs(
  read_count = read_count_raw,
  outdir = outdir,
  known_occurrences = known_occurrences,
  min_read_count_cutoff = 10,
  max_read_count_cutoff = 100,
  increment_read_count_cutoff = 10,
  min_variant_cutoff = 0.001,
  max_variant_cutoff = 0.05,
  increment_variant_cutoff = 0.005,
  by_replicate = by_replicate,
  min_replicate_number = 2,
  log_file = log_file
)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "suggest_variant_readcount_cutoffs", elapsed = timing[3], user = timing[1], system =timing[2])

# Print stats ----
runtime_file = file.path(outdir, "runtimes_optimize.csv")
write.csv(run_time, file=runtime_file, row.names = FALSE)


