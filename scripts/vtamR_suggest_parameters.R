args <- commandArgs(trailingOnly = TRUE)

dataset <- as.character(args[1])
vsearch_path  <- as.character(args[2])
pcr_error_var_prop <- as.numeric(args[3])
sample_cutoff<- as.numeric(args[4])
by_replicate <- as.logical(args[5])
outdir_pipe <- as.character(args[6])

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
# To establish known_occurrences df that is used in the suggest function
read_count_prefilter_file <- file.path(outdir_pipe, dataset, "vtamR", "filter_default", "12_pool_replicates.csv")

run_time <-  data.frame(
  step = character(),
  elapsed = numeric(),
  user = numeric(),
  system = numeric()
)


############### suggest 
# This workflow is intended to match as much as possible the VTAM workflow
# 1. run suggest_sample_cutoff and suggest_pcr_error_cutoff on raw unfiltered data
# Both functions depend only of mock composition, and do not need kwown_occurrences
# 2. run filter_occurrences_cutoff and filter_pcr_error with the suggested parameters on unfiltered data => read_count_tmp
# 3. run classify_control_occurrences to make know_occurrences file (with TP, FP) 
# based on data filtered with default values
# 4. run suggest_variant_read_count on read_count_tmp data using know_occurrences established on default filtered data (point3)


## 1. suggest_sample_cutoff and suggest_pcr_error_cutoff ----
start_time <- proc.time()
outfile = file.path(outdir, "suggest_sample_cutoff.csv")

suggest_sample_cutoff_df <- suggest_sample_cutoff(
  read_count=read_count_raw,
  mock_composition=mock_composition,
  outfile=outfile
)
head(suggest_sample_cutoff_df)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "suggest_sample_cutoff", elapsed = timing[3], user = timing[1], system =timing[2])


start_time <- proc.time()
outfile = file.path(outdir, "suggest_pcr_error_cutoff.csv")

suggest_pcr_error_cutoff_df <- suggest_pcr_error_cutoff(
  read_count = read_count_raw, 
  mock_composition=mock_composition,
  vsearch_path=vsearch_path, 
  outfile=outfile
)
head(suggest_pcr_error_cutoff_df)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "suggest_pcr_error_cutoff", elapsed = timing[3], user = timing[1], system =timing[2])

# 2. run filter_occurrences_cutoff and filter_pcr_error with the suggested parameters on unfiltered data => read_count_tmp
# filter_occurrence_sample
start_time <- proc.time()
read_count_df_sample <- filter_occurrence_sample(
  read_count_raw, 
  cutoff=sample_cutoff
)
# filter_pcr_error
by_sample <- TRUE
sample_prop <- 0.8
read_count_df_pcr <- filter_pcr_error(
  read_count_raw,
  pcr_error_var_prop = pcr_error_var_prop,
  sample_prop = sample_prop,
  vsearch_path = vsearch_path
)
# pool_filters
read_count_df_tmp <- pool_filters(
  read_count_df_sample, 
  read_count_df_pcr
  )
# filter_min_replicate
min_replicate_number <- 2
read_count_df_tmp <- filter_min_replicate(
  read_count_df_tmp, 
  cutoff=min_replicate_number
)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "prefilter_sample_pcr", elapsed = timing[3], user = timing[1], system =timing[2])

# 3. run classify_control_occurrences to make know_occurrences file (with TP, FP) 
# based on data filtered with default values
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
  performance_metrics = performance_metrics
)

known_occurrences <- results[[1]]

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "classify_control_occurrences", elapsed = timing[3], user = timing[1], system =timing[2])


###################"

# 4. run suggest_variant_read_count on read_count_tmp data using know_occurrences established on default filtered data (point3)
start_time <- proc.time()

suggest_variant_readcount_cutoffs_df <- suggest_variant_readcount_cutoffs(
  read_count = read_count_df_tmp,
  outdir = outdir,
  known_occurrences = known_occurrences,
  min_read_count_cutoff = 10,
  max_read_count_cutoff = 100,
  increment_read_count_cutoff = 10,
  min_variant_cutoff = 0.001,
  max_variant_cutoff = 0.05,
  increment_variant_cutoff = 0.005,
  by_replicate = by_replicate,
  min_replicate_number = 2
)

head(suggest_variant_readcount_cutoffs_df)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "suggest_variant_readcount_cutoffs", elapsed = timing[3], user = timing[1], system =timing[2])



# Print stats ----

runtime_file = file.path(outdir, "runtimes_optimize.csv")
write.csv(run_time, file=runtime_file, row.names = FALSE)


