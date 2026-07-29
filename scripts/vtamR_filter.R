args <- commandArgs(trailingOnly = TRUE)

dataset <- as.character(args[1])
filterdir <- as.character(args[2])
vsearch_path <- as.character(args[3])
cutadapt_path <- as.character(args[4])
pcr_error_var_prop <- as.numeric(args[5]) # 0.1
sample_cutoff <- as.numeric(args[6]) # 0.001
read_count_cutoff <-as.numeric(args[7]) # 10
variant_cutoff <- as.numeric(args[8]) # 0.001
by_replicate <- as.logical(args[9]) # TRUE
coding <- as.logical(args[10]) # TRUE
outdir_pipe <- as.character(args[11])

library(vtamR)
library(dplyr)

search()
packageVersion("vtamR")

# setwd("~/vtamR_benchmark")


asv_list = file.path(outdir_pipe, "asv_all_datasets.csv")
asv_list_df = read.csv(asv_list)
colnames(asv_list_df) <- c("asv_id", "asv")

outdir = file.path(outdir_pipe, dataset, "vtamR", filterdir)
mock_composition <- paste("vtamR_mock_composition_", dataset, ".csv", sep="")
mock_composition <- file.path("data", "metadata", mock_composition)


if(dataset == "bat"){
  sampleinfo <- file.path("data", "fastq",  dataset, "demultiplex_trim", "sampleinfo.csv")
  demultiplexed_dir <- file.path("data", "fastq",  dataset, "demultiplex_trim")
}else{
  sampleinfo <- file.path(outdir_pipe, dataset, "vtamR", "demultiplex_trim", "sampleinfo.csv")
  demultiplexed_dir <- file.path(outdir_pipe, dataset, "vtamR", "demultiplex_trim")
}


run_time <-  data.frame(
  step = character(),
  elapsed = numeric(),
  user = numeric(),
  system = numeric()
)

stat_df <- data.frame(
  parameters = character(),
  asv_count = integer(),
  read_count = integer(),
  sample_count = integer(),
  sample_replicate_count = integer()
)


# Filter ##################################

## 1. SampleInformationDereplicate ----
start_time <- proc.time()
outfile <- file.path(outdir, "1_before_filter.csv")
read_count_df <- dereplicate(
  sampleinfo = sampleinfo, 
  outfile=outfile,
  dir=demultiplexed_dir,
  input_asv_list = asv_list_df
  # output_asv_list=updated_asv_list
)
stat_df <- get_stat(read_count_df, stat_df=stat_df, stage="input", params=NA)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "dereplicate", elapsed = timing[3], user = timing[1], system =timing[2])

## 2.VariantReadCount ----
start_time <- proc.time()
global_read_count_cutoff = 2
outfile <- file.path(outdir, "2_LFNglobalReadCount.csv")
read_count_df <- filter_asv_global(
  read_count_df, 
  outfile=outfile,
  cutoff=global_read_count_cutoff
)
param <- paste("global_read_count_cutoff:", global_read_count_cutoff, sep=" ")
stat_df <- get_stat(read_count_df, stat_df, stage="filter_asv_global", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_asv_global", elapsed = timing[3], user = timing[1], system =timing[2])

## 3. FilterLFN ----

### LFN sample ----
start_time <- proc.time()
outfile <- file.path(outdir, "3-1_filter_occurrence_sample.csv")
read_count_df_sample <- filter_occurrence_sample(
  read_count_df, 
  cutoff=sample_cutoff, 
  outfile=outfile
)

param <- paste("cutoff:", sample_cutoff, sep=" ")
stat_df <- get_stat(read_count_df_sample, stat_df, stage="filter_occurrence_sample", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_occurrence_sample", elapsed = timing[3], user = timing[1], system =timing[2])

### LFN variant ----
start_time <- proc.time()
outfile <- file.path(outdir, "3-2_filter_occurrence_variant.csv")
read_count_df_variant <- filter_occurrence_variant(
  read_count_df, 
  cutoff=variant_cutoff, 
  outfile=outfile,
  by_replicate = by_replicate
)

param <- paste("cutoff:", variant_cutoff, "by_replicate:", by_replicate, sep=" ")         
stat_df <- get_stat(read_count_df_variant, stat_df, stage="filter_occurrence_variant", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_occurrence_variant", elapsed = timing[3], user = timing[1], system =timing[2])

### LFN read_count ----
start_time <- proc.time()
outfile <- file.path(outdir, "3-3_filter_occurrence_read_count.csv")
read_count_df_rc <- filter_occurrence_read_count(
  read_count_df, 
  cutoff=read_count_cutoff, 
  outfile=outfile
)

param <- paste("cutoff:", read_count_cutoff, sep=" ")        
stat_df <- get_stat(read_count_df_rc, stat_df, stage="filter_occurrence_read_count", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_occurrence_read_count", elapsed = timing[3], user = timing[1], system =timing[2])

##### Pool LFN
start_time <- proc.time()
outfile <- file.path(outdir, "3-4_pool_lfn_filters.csv")

read_count_df <- pool_filters(
  read_count_df_rc, read_count_df_variant, read_count_df_sample,
  outfile=outfile
)

stat_df <- get_stat(read_count_df, stat_df, stage="pool_filters", params=NA)

# delete temporary data frames
rm(read_count_df_rc)
rm(read_count_df_variant)
rm(read_count_df_sample)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "pool_filters", elapsed = timing[3], user = timing[1], system =timing[2])

## 4. FilterMinReplicateNumber ----
start_time <- proc.time()
min_replicate_number <- 2
outfile <- file.path(outdir, "4_filter_min_replicate_1.csv")

read_count_df <- filter_min_replicate(
  read_count_df, 
  cutoff=min_replicate_number, 
  outfile=outfile
)

param <- paste("cutoff:", min_replicate_number, sep=" ")  
stat_df <- get_stat(read_count_df, stat_df, stage="filter_min_replicate_1", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_min_replicate_1", elapsed = timing[3], user = timing[1], system =timing[2])

## 5. FilterPCRerror ----

start_time <- proc.time()
outfile <- file.path(outdir, "5_filter_pcr_error.csv")
by_sample <- TRUE
sample_prop <- 0.8
filter_occurrence <- TRUE
min_read_count <- 10

read_count_df <- filter_pcr_error(
  read_count_df,
  outfile = outfile,
  pcr_error_var_prop = pcr_error_var_prop,
  sample_prop = sample_prop,
  vsearch_path = vsearch_path,
  filter_occurrence = filter_occurrence,
  min_read_count = min_read_count
)

param <- paste("pcr_error_var_prop:", pcr_error_var_prop, "by_sample:", by_sample, "sample_prop:", sample_prop, sep=" ")
stat_df <- get_stat(read_count_df, stat_df, stage="filter_pcr_error", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_pcr_error", elapsed = timing[3], user = timing[1], system =timing[2])

## 6.FilterChimera ----

start_time <- proc.time()
outfile <- file.path(outdir, "6_filter_chimera.csv")
abskew=16
by_sample = TRUE
sample_prop <- 0.8
filter_occurrence <- TRUE

read_count_df <- filter_chimera(
  read_count_df, 
  outfile=outfile,
  by_sample=by_sample, 
  abskew=abskew,
  vsearch_path = vsearch_path,
  filter_occurrence = filter_occurrence
)

param <- paste("by_sample:", by_sample, "abskew:", abskew, "sample_prop:", sample_prop, sep=" ")
stat_df <- get_stat(read_count_df, stat_df, stage="filter_chimera", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_chimera", elapsed = timing[3], user = timing[1], system =timing[2])


## 7. FilterMinReplicateNumber2 ----

start_time <- proc.time()
min_replicate_number <- 2
outfile <- file.path(outdir, "7_filter_min_replicate_2.csv")

read_count_df <- filter_min_replicate(
  read_count_df, 
  cutoff=min_replicate_number, 
  outfile=outfile
)

param <- paste("cutoff:", min_replicate_number, sep=" ")  
stat_df <- get_stat(read_count_df, stat_df, stage="filter_min_replicate_2", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_min_replicate_2", elapsed = timing[3], user = timing[1], system =timing[2])

## 8. FilterRenkonen ----

start_time <- proc.time()
outfile <- file.path(outdir, "8_filter_replicate.csv")
renkonen_distance_quantile <- 0.9

read_count_df <- filter_replicate(
  read_count_df, 
  renkonen_distance_quantile = renkonen_distance_quantile,
  outfile = outfile
)

param <- paste("renkonen_distance_quantile:", renkonen_distance_quantile, sep=" ")
stat_df <- get_stat(read_count_df, stat_df, stage="filter_replicate", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_replicate", elapsed = timing[3], user = timing[1], system =timing[2])


## 9. FilterMinReplicateNumber3 ----

start_time <- proc.time()
min_replicate_number <- 2
outfile <- file.path(outdir, "9_filter_min_replicate_3.csv")

read_count_df <- filter_min_replicate(
  read_count_df, 
  cutoff=min_replicate_number, 
  outfile=outfile
)

param <- paste("cutoff:", min_replicate_number, sep=" ")  
stat_df <- get_stat(read_count_df, stat_df, stage="filter_min_replicate_3", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "filter_min_replicate_3", elapsed = timing[3], user = timing[1], system =timing[2])

if(coding){
  ## 10. FilterIndel ----
  
  start_time <- proc.time()
  outfile <- file.path(outdir, "10_filter_indel.csv")
  
  read_count_df <- filter_indel(
    read_count_df,
    outfile=outfile
  )
  
  stat_df <- get_stat(read_count_df, stat_df, stage="filter_indel", params=NA)
  
  
  end_time <- proc.time()
  timing <- end_time - start_time
  run_time <- add_row(run_time, step = "filter_indel", elapsed = timing[3], user = timing[1], system =timing[2])
  
  
  ## 11. FilterCodonStop ----
  
  start_time <- proc.time()
  outfile <- file.path(outdir, "11_filter_stop_codon.csv")
  
  genetic_code = 5
  read_count_df <- filter_stop_codon(
    read_count_df, 
    genetic_code=genetic_code,
    outfile=outfile
  )
  
  param <- paste("genetic_code:", genetic_code, sep=" ")
  stat_df <- get_stat(read_count_df, stat_df, stage="filter_stop_codon", params=param)
  
  end_time <- proc.time()
  timing <- end_time - start_time
  run_time <- add_row(run_time, step = "filter_stop_codon", elapsed = timing[3], user = timing[1], system =timing[2])
}

## 12. ReadCountAverageOverReplicates ----
#vtam write sum readcount over replicates and not mean to asv_table

start_time <- proc.time()
outfile <- file.path(outdir, "12_pool_replicates.csv")
method = "sum"

read_count_sample_df <- pool_replicates(
  read_count_df,
  method = method,
  outfile = outfile
)

param <- paste("method:", method, sep=" ")
stat_df <- get_stat(read_count_sample_df, stat_df, stage="pool_replicates", params=param)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "pool_replicates", elapsed = timing[3], user = timing[1], system =timing[2])


# Performance metrics ----

start_time <- proc.time()
outfile=file.path(outdir, "asv_table.csv")

asv_table_df <- write_asv_table(read_count_sample_df, 
                                outfile=outfile, 
                                sampleinfo=sampleinfo, 
                                add_empty_samples=TRUE, 
                                add_sums_by_sample=TRUE, 
                                add_sums_by_asv=TRUE, 
                                add_expected_asv=TRUE, 
                                mock_composition=mock_composition)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "write_asv_table", elapsed = timing[3], user = timing[1], system =timing[2])


start_time <- proc.time()
false_negatives <- file.path(outdir, "performance", "false_negatives.csv")
performance_metrics <- file.path(outdir, "performance", "performance_metrics.csv")
known_occurrences <- file.path(outdir, "performance", "known_occurrences.csv")

results <- classify_control_occurrences(read_count_sample_df,
                                        sampleinfo=sampleinfo, 
                                        mock_composition=mock_composition, 
                                        known_occurrences=known_occurrences, 
                                        false_negatives=false_negatives,
                                        performance_metrics=performance_metrics)

# give explicit names to the 3 output data frames
known_occurrences_df <- results[[1]]
false_negatives_df <- results[[2]]
performance_metrics_df <- results[[3]]

print(performance_metrics_df)

end_time <- proc.time()
timing <- end_time - start_time
run_time <- add_row(run_time, step = "classify_control_occurrences", elapsed = timing[3], user = timing[1], system =timing[2])


# Print stats ----

stat_file = file.path(outdir, "performance", "stat_read_counts.csv")
write.csv(stat_df, file=stat_file, row.names = TRUE)

runtime_file = file.path(outdir, "performance", "runtimes.csv")
write.csv(run_time, file=runtime_file, row.names = FALSE)


