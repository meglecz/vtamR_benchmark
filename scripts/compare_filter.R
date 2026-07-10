library(vtamR)
library(dplyr)
library(tidyr)

setwd("~/vtamR_benchmark")

datasets =  c("fish_mfzr", "fish_zfzr", "shark", "bat")
dir = "pipeline_output_1"
filter_dir <- "filter_default"

outdir <- file.path("comparaison_results", "compare_filter")
outdir <- check_dir(path=outdir, is_file=FALSE)

##### read vtam results (extracted from db.sqlite) ----
vtam_sample <- read.csv("comparaison_results/db_tables/samples.csv")
vtam_run <- read.csv("comparaison_results/db_tables/runs.csv")
vtam_marker <- read.csv("comparaison_results/db_tables/markers.csv")
vtam_VariantReadCount <- read.csv("comparaison_results/db_tables/VariantReadCount.csv")
vtam_FilterLFN <- read.csv("comparaison_results/db_tables/FilterLFN.csv")
vtam_FilterMinReplicateNumber <- read.csv("comparaison_results/db_tables/FilterMinReplicateNumber.csv")
vtam_FilterPCRerror <- read.csv("comparaison_results/db_tables/FilterPCRerror.csv")
vtam_FilterChimera <- read.csv("comparaison_results/db_tables/FilterChimera.csv")
vtam_FilterMinReplicateNumber2 <- read.csv("comparaison_results/db_tables/FilterMinReplicateNumber2.csv")
vtam_FilterRenkonen <- read.csv("comparaison_results/db_tables/FilterRenkonen.csv")
vtam_FilterMinReplicateNumber3 <- read.csv("comparaison_results/db_tables/FilterMinReplicateNumber3.csv")
vtam_FilterIndel <- read.csv("comparaison_results/db_tables/FilterIndel.csv")
vtam_FilterCodonStop <- read.csv("comparaison_results/db_tables/FilterCodonStop.csv")

for(dataset in datasets){

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
    select(asv_id_vtam, asv_id_vtamR, sample, read_count_vtam, read_count_vtamR, diff, asv) %>%
    arrange(diff, asv, sample)
  
  ##### missing occurrences in vtam ----
  
  # get run_id and marker_id
  if(dataset %in% c("fish_mfzr", "fish_zfzr")){
    vtam_run_id <- vtam_run[vtam_run$name == "fish", "id"]
    marker <- sub("fish_", "", dataset)
    vtam_marker_id <- vtam_marker[vtam_marker$name == marker, "id"]
  }else {
   vtam_run_id <- vtam_run[vtam_run$name == dataset, "id"]
   if(dataset == "bat"){
     vtam_marker_id <- vtam_marker[vtam_marker$name == "coi", "id"]
   }else{
     vtam_marker_id <- vtam_marker[vtam_marker$name == "16S", "id"]
   }
  }

  
  missing_vtam <- rc_filter %>%
    filter(diff < 0) %>%
    arrange(asv_id_vtamR, sample) %>%
    select(asv_id = asv_id_vtamR, sample, read_count_vtam, read_count_vtamR, asv)
  
  # add sample id
  missing_vtam <- missing_vtam %>%
    left_join(vtam_sample, by=c("sample" = "name")) %>%
    rename(sample_id = id)
  
  # decompose LFN filters

  vtam_FilterLFN_tmp <- vtam_FilterLFN %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0 ) %>% # select the occurrences not filtered out
    select(sample_id, asv_id = variant_id, filter_id) %>%
    distinct() # if one replicate is present it is still OK
  
  # For bat lfn_per_variant, for the others lfn_per_variant_replicate was used
  if(dataset == "bat"){
    lfn_per_variant_id = 2
  } else{
    lfn_per_variant_id = 3
  }
  
  lfn_per_variant <- vtam_FilterLFN_tmp %>%
    filter(filter_id == lfn_per_variant_id) %>%
    rename(lfn_per_variant = filter_id)
  
  lfn_sample_replicate <- vtam_FilterLFN_tmp %>%
    filter(filter_id == 6) %>%
    rename(lfn_sample_replicate = filter_id)
  
  lfn_read_count <- vtam_FilterLFN_tmp %>%
    filter(filter_id == 7) %>%
    rename(lfn_read_count = filter_id)
  
  lfn <- vtam_FilterLFN_tmp %>%
    filter(filter_id == 8) %>%
    rename(lfn = filter_id)
  
  missing_vtam <- missing_vtam %>%
    left_join(lfn_per_variant, by=c("sample_id", "asv_id")) %>%
    left_join(lfn_sample_replicate, by=c("sample_id", "asv_id")) %>%
    left_join(lfn_read_count, by=c("sample_id",  "asv_id")) %>%
    left_join(lfn, by=c("sample_id",  "asv_id"))
    

  vtam_FilterMinReplicateNumber_tmp <- vtam_FilterMinReplicateNumber %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterMinReplicateNumber = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterMinReplicateNumber_tmp, by=c("sample_id", "asv_id"))
  
  
  vtam_FilterPCRerror_tmp <- vtam_FilterPCRerror %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterPCRerror = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterPCRerror_tmp, by=c("sample_id", "asv_id"))
  
  vtam_FilterChimera_tmp <-vtam_FilterChimera %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterChimera = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterChimera_tmp, by=c("sample_id", "asv_id"))
  
  vtam_FilterMinReplicateNumber2_tmp <- vtam_FilterMinReplicateNumber2 %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterMinReplicateNumber2 = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterMinReplicateNumber2_tmp, by=c("sample_id", "asv_id"))
  
  vtam_FilterRenkonen_tmp <- vtam_FilterRenkonen %>%
    filter(run_id == vtam_run_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, vtam_FilterRenkonen = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterRenkonen_tmp, by=c("sample_id", "asv_id"))
  
  
  vtam_FilterMinReplicateNumber3_tmp <- vtam_FilterMinReplicateNumber3 %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterMinReplicateNumber3 = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterMinReplicateNumber3_tmp, by=c("sample_id", "asv_id"))
  
  
  vtam_FilterIndel_tmp <- vtam_FilterIndel %>%
    filter(run_id == vtam_run_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterIndel = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterIndel_tmp, by=c("sample_id", "asv_id"))
  
  vtam_FilterCodonStop_tmp <- vtam_FilterCodonStop %>%
    filter(run_id == vtam_run_id,
           marker_id == vtam_marker_id,
           filter_delete == 0) %>%
    select(sample_id, asv_id = variant_id, FilterCodonStop = filter_delete) %>%
    distinct()
  missing_vtam <- missing_vtam %>%
    left_join(vtam_FilterCodonStop_tmp, by=c("sample_id", "asv_id"))
  
  # change values for TRUE/FALSE for presence of asv in sample after each sample
  missing_vtam <- missing_vtam %>%
    mutate(
      across(lfn_per_variant:FilterCodonStop, ~ !is.na(.x))
    ) %>%
    select(-asv, everything(), asv)
    
  write.csv(missing_vtam, file=file.path(outdir, paste0("occurrences_missing_in_vtam_", dataset, ".csv")), row.names = FALSE)
  write.csv(rc_filter, file=file.path(outdir, paste0("read_count_VTAM_vtamR_", dataset, ".csv")), row.names = FALSE)
  
  
  ###############################################
  ## missing occurrences in vtamR ----
  # select sv_id, sample combination that are missing from vtamR but present in VTAM
  missing_vtamR <- rc_filter %>%
    filter(diff > 0) %>%
    select(asv_id = asv_id_vtam, sample) %>%
    mutate(flag=TRUE)
  
  dir_temp <- file.path(dir, dataset, "vtamR", filter_dir)
  hist_vtamR <- history_by(dir=dir_temp, feature = "asv_id", values= unique(missing_vtamR$asv_id)) %>%
    left_join(missing_vtamR, by = c("asv_id", "sample")) %>%
    filter(flag==TRUE) %>%
    select(-replicate, -read_count, -flag) %>%
    distinct() %>%
    arrange(asv_id, file, sample)
  
  write.csv(hist_vtamR, file=file.path(outdir, paste0("occurrences_missing_in_vtamR_", dataset, ".csv")), row.names = FALSE)
}


