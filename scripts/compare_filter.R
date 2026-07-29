library(vtamR)
library(dplyr)
library(tidyr)

setwd("~/vtamR_benchmark")

datasets =  c("fish_mfzr", "fish_zfzr", "shark", "bat")
# The content of the 10 runs are identical. Check only one of them
dir = "pipeline_output_1" 
filter_dir <- "filter_default"

outdir <- file.path("comparaison_results", "filter")
outdir <- check_dir(path=outdir, is_file=FALSE)
pooled_rc <- data.frame()

for(dataset in datasets){

  # outputs of the final filtering step for VTAM and vtamR
  vtamR_file <- file.path(dir, dataset, "vtamR",  filter_dir, "12_pool_replicates.csv")
  vtam_file <- file.path(dir, dataset, "VTAM",  filter_dir, "asvtable_default.tsv")
  
  ##### pool vtam and vtamR ----
  vtamR_rc_long <- read.csv(vtamR_file) %>%
    rename(read_count_vtamR = read_count,
           asv_id_vtamR = asv_id)
  
  vtam_rc_long <- read.table(vtam_file, sep="\t", header=TRUE) %>%
    select(-run, -marker, -sequence_length, -read_count, -clusterid, -clustersize, -chimera_borderline) %>%
    rename(asv = sequence,
           asv_id_vtam = variant) %>%
    pivot_longer(           
      cols = c(-asv, -asv_id_vtam),       
      names_to = "sample",  
      values_to = "read_count_vtam"  
    ) %>%
    filter(read_count_vtam > 0) %>%
    mutate(sample = sub("^X", "", sample)) %>% # if sample name starts by a number read.table adds an X before
    select(asv_id_vtam, sample, read_count_vtam, asv) 
  
  rc_filter <- full_join(vtam_rc_long, vtamR_rc_long, by=c("asv", "sample")) %>%
    mutate(read_count_vtam = if_else(is.na(read_count_vtam), 0, read_count_vtam)) %>%
    mutate(read_count_vtamR = if_else(is.na(read_count_vtamR), 0, read_count_vtamR)) %>%
    mutate(diff = read_count_vtam-read_count_vtamR) %>%
    mutate("dataset" = dataset) %>%
    select("dataset", asv_id_vtam, asv_id_vtamR, sample, read_count_vtam, read_count_vtamR, diff, asv) %>%
    arrange(diff, asv, sample)
  
  
  pooled_rc <- rbind(pooled_rc, rc_filter)
}


write.csv(pooled_rc, file=file.path(outdir, "filter_read_counts.csv"), row.names = FALSE)

diff <- pooled_rc %>%
  filter(diff !=0)
print("Occurrences with different read counts:")
print(diff)


