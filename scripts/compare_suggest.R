library(vtamR)
library(dplyr)
library(tidyr)


setwd("~/vtamR_benchmark")

datasets =  c("fish_mfzr", "fish_zfzr", "shark", "bat")
# The content of the 10 runs are identical. Check only one of them
dir = "pipeline_output_1" 
suggest_dir <- "optimize_suggest"
outdir <- file.path("comparaison_results", suggest_dir)
outdir <- check_dir(path=outdir, is_file=FALSE)
asv_all <- file.path(dir, "asv_all_datasets.csv")

########### known_occurrences
known_occ <- data.frame(
  "dataset" = character(),
  "action" = character(), 
  "sample" = character(), 
  "asv" = character(),
  "present_vtam"= character(),
  "present_vtamR"= character()
  )
for(dataset in datasets){
  
  print("#################")
  print(paste0("DATASET:", dataset))
  
  filename <- paste0("known_occurrences_", dataset, ".tsv")
  known_occurrences_vtam <- file.path(dir, dataset, "VTAM", suggest_dir, filename)
  known_occurrences_vtam_df <- read.table(known_occurrences_vtam, sep="\t", header = TRUE) %>%
    rename(asv = sequence) %>%
    select(sample, action, asv) %>%
    arrange(action, sample, asv)
  known_occurrences_vtamR <- file.path(dir, dataset, "vtamR", suggest_dir, "known_occurrences.csv" )
  known_occurrences_vtamR_df <- read.table(known_occurrences_vtamR, sep=",", header = TRUE) %>%
    select(sample, action, asv) %>%
    arrange(action, sample, asv)
  
  # test if dfs are identical
  ident <- identical(known_occurrences_vtam_df, known_occurrences_vtamR_df)
  # add vtam and vtamR columns and make full join to see which lines are different (if any)
  known_occurrences_vtam_df <- known_occurrences_vtam_df  %>%
      mutate(vtam=TRUE)
  known_occurrences_vtamR_df <- known_occurrences_vtamR_df  %>%
      mutate(vtamR=TRUE)
    
  tmp <- full_join(known_occurrences_vtam_df, known_occurrences_vtamR_df, 
                     by = c("action", "sample", "asv")) %>%
    mutate(dataset = dataset) %>%
    select(dataset, everything())
  
  known_occ <- rbind(known_occ, tmp)
  
  # see only the lines that are different
  tmp <- tmp  %>%
      filter(is.na(vtam) | is.na(vtamR))
    
  msg <- paste0("known_occurrences: ", ident)
  print(msg)
  print(tmp)
  
}
write.csv(known_occ, file = file.path(outdir, "known_occurrences.csv"), row.names = FALSE)


################## optimize_lfn_sample_replicate

sample_repl <- data.frame(
  "dataset" = character(),
  "sample" = character(), 
  "replicate" = character(), 
  "asv_id" = character(),
  "sample_cutoff" = numeric(),
  "present_vtam"= character(),
  "present_vtamR"= character()
)
for(dataset in datasets){
  
  print("#################")
  print(paste0("DATASET:", dataset))

  vtam <- file.path(dir, dataset, "VTAM", suggest_dir, "optimize_lfn_sample_replicate.tsv")
  vtam_df <-  read.table(vtam, sep="\t", header = TRUE)  %>%
    rename(asv_id = variant, sample_cutoff = round_down) %>%
    select(sample, replicate, asv_id, sample_cutoff) %>%
    arrange(sample, replicate, asv_id, sample_cutoff)
  
  vtamR <- file.path(dir, dataset, "vtamR", suggest_dir, "suggest_sample_cutoff.csv")
  vtamR_df <- read.table(vtamR, sep=",", header = TRUE) %>%
    select(sample, replicate, asv_id, sample_cutoff) %>%
    arrange(sample, replicate, asv_id, sample_cutoff)
  
  ident <- identical(vtam_df, vtamR_df)


  vtam_df <- vtam_df  %>%
      mutate(vtam=TRUE)
  vtamR_df <- vtamR_df  %>%
      mutate(vtamR=TRUE)
    
  tmp <- full_join(vtam_df, vtamR_df, 
                     by = c("sample", "replicate", "asv_id", "sample_cutoff"))%>%
    mutate(dataset = dataset) %>%
    select(dataset, everything())
  
  sample_repl <- rbind(sample_repl, tmp)
  
  tmp <- tmp %>%
      filter(is.na(vtam) | is.na(vtamR))
    
    msg <- paste0("optimize_lfn_sample_replicate: ", ident)
    print(msg)
    print(tmp)
}
write.csv(sample_repl, file = file.path(outdir, "sample_replicate_cutoffs.csv"), row.names = FALSE)

################## optimize_pcr_error
pcr_error <- data.frame(
  "dataset" = character(),
  "sample" = character(), 
  "expected_asv_id" = character(), 
  "unexpected_asv_id" = character(),
  "pcr_error_var_prop" = numeric(),
  "unexpected_read_count" = numeric(),
  "present_vtam"= character(),
  "present_vtamR"= character()
)
for(dataset in datasets){
  
  print("#################")
  print(paste0("DATASET:", dataset))

  vtam <- file.path(dir, dataset, "VTAM", suggest_dir, "optimize_pcr_error.tsv")
  vtam_df <-  read.table(vtam, sep="\t", header = TRUE)  %>%
    rename(expected_asv_id = variant_expected, 
           unexpected_asv_id = variant_unexpected, 
           pcr_error_var_prop = N_ij_unexpected_to_expected_ratio,
           unexpected_read_count = N_ij_unexpected) %>%
    mutate(pcr_error_var_prop = round(pcr_error_var_prop, digits = 5)) %>%
    select(sample, expected_asv_id, unexpected_asv_id, pcr_error_var_prop, unexpected_read_count) %>%
    arrange(sample, expected_asv_id, unexpected_asv_id, pcr_error_var_prop, unexpected_read_count)
  
  vtamR <- file.path(dir, dataset, "vtamR", suggest_dir, "suggest_pcr_error_cutoff.csv")
  vtamR_df <- read.table(vtamR, sep=",", header = TRUE) %>%
    select(sample, expected_asv_id, unexpected_asv_id, pcr_error_var_prop, unexpected_read_count) %>%
    mutate(pcr_error_var_prop = round(pcr_error_var_prop, digits = 5)) %>%
    arrange(sample, expected_asv_id, unexpected_asv_id, pcr_error_var_prop, unexpected_read_count)
  
  ident <- identical(vtam_df, vtamR_df)
  
  vtam_df <- vtam_df  %>%
      mutate(vtam=TRUE)
  vtamR_df <- vtamR_df  %>%
      mutate(vtamR=TRUE)
    
  tmp <- full_join(vtam_df, vtamR_df, 
                   by = c("sample", "expected_asv_id", "unexpected_asv_id", "pcr_error_var_prop", "unexpected_read_count" )) %>%
    mutate(dataset = dataset) %>%
    select(dataset, everything())
  
  pcr_error <- rbind(pcr_error, tmp)
  
  tmp <- tmp%>%
      filter(is.na(vtam) | is.na(vtamR)) %>%
      filter(unexpected_read_count > 10)
    
    msg <- paste0("optimize_pcr_error: ", ident)
    print(msg)
    print(tmp)
}
write.csv(pcr_error, file = file.path(outdir, "pcr_error_cutoffs.csv"), row.names = FALSE)

### In vtamR, samples occurrences with fewer than 10 reads are ignored to avoid 
# estimating the error proportion from low read counts, which are more susceptible 
# to random fluctuations and would likely be removed by other filters anyway. 
# This improves robustness and reduces runtime.
### vtamR identifies ASV pairs as PCR errors when they differ by one nucleotide 
# in length while being otherwise identical. VTAM does not apply this criterion.

# extract discordant results between vtam and vtamR, where unexpected_read_count >= 10
df <- pcr_error %>%
  filter((is.na(vtamR) & !is.na(vtam)) | (!is.na(vtamR) & is.na(vtam))) %>%
  filter(unexpected_read_count >= 10) %>%
  select(expected_asv_id, unexpected_asv_id) %>%
  distinct()

# add sequence
asv_all_df <- read.csv(asv_all)
df <- left_join(df, asv_all_df, by=c("expected_asv_id" = "id")) %>%
  rename(expected_asv = sequence) %>%
  left_join(asv_all_df, by=c("unexpected_asv_id" = "id")) %>%
  rename(unexpected_asv = sequence)

# write files with sequence pairs and align them
aligned_dir <- file.path(outdir, "alignments")
for(i in 1:nrow(df)){
  outfile <- paste(df[i, "expected_asv_id"], df[i, "unexpected_asv_id"], sep="-")
  outfile <- file.path(aligned_dir, outfile)
  
  con <- file(outfile, "wt")
  txt <- paste0(">", df[i, "expected_asv_id"], "\n", df[i, "expected_asv"], "\n", ">", df[i, "unexpected_asv_id"], "\n", df[i, "unexpected_asv"], "\n")
  writeLines(txt, con)
  close(con)
  
  aligned <- paste0(outfile, "_aligned.fas")
  cmd <- paste0("muscle -in ", outfile," -out ", aligned)
  system(cmd)
}

count <- pcr_error %>%
  group_by(vtam, vtamR) %>%
  summarise(Nb_occurrences = n())
  


################## optimize_lfn_read_count_and_lfn_variant
lfn_read_count_and_lfn_variant <- data.frame(
  "dataset" = character(),
  "read_count_cutoff" = numeric(),
  "variant_cutoff" = numeric(),
  "TP_vtam" = numeric(),
  "FP_vtam" = numeric(),
  "TP_vtamR" = numeric(),
  "FP_vtamR" = numeric(),
  "diff_TP" = numeric(),
  "diff_FP" = numeric()
)
for(dataset in datasets){
  
  print("#################")
  print(paste0("DATASET:", dataset))
    ################## optimize_lfn_read_count_and_lfn_variant
  
  if(dataset == "bat"){
    vtam <- file.path(dir, dataset, "VTAM", suggest_dir, "optimize_lfn_read_count_and_lfn_variant.tsv")
    vtam_df <-  read.table(vtam, sep="\t", header = TRUE)  %>%
      rename(read_count_cutoff =lfn_nijk_cutoff,
             variant_cutoff = lfn_variant_cutoff,
             TP	= occurrence_nb_keep,
             FP = occurrence_nb_delete
      )  %>%
      select(read_count_cutoff, variant_cutoff, TP, FP) %>%
      mutate(variant_cutoff = 
               if_else(variant_cutoff > 0.031, variant_cutoff - 0.001, variant_cutoff)) %>% # correct vtam rounding
      arrange(read_count_cutoff, variant_cutoff, TP, FP) 
  }else{
    vtam <- file.path(dir, dataset, "VTAM", suggest_dir, "optimize_lfn_read_count_and_lfn_variant_replicate.tsv")
    vtam_df <-  read.table(vtam, sep="\t", header = TRUE)  %>%
      rename(read_count_cutoff =lfn_nijk_cutoff,
             variant_cutoff = lfn_variant_replicate_cutoff,
             TP	= occurrence_nb_keep,
             FP = occurrence_nb_delete
      ) %>%
      select(read_count_cutoff, variant_cutoff, TP, FP) %>%
      mutate(variant_cutoff = 
               if_else(variant_cutoff > 0.031, variant_cutoff - 0.001, variant_cutoff)) %>% # correct vtam rounding
      arrange(read_count_cutoff, variant_cutoff, TP, FP) 
  }

  vtamR <- file.path(dir, dataset, "vtamR", suggest_dir, "suggest_variant_readcount_cutoffs.csv")
  vtamR_df <- read.table(vtamR, sep=",", header = TRUE) %>%
    select(read_count_cutoff, variant_cutoff, TP, FP)  %>%
    arrange(read_count_cutoff, variant_cutoff, TP, FP) 
  
  ident <- identical(vtam_df, vtamR_df)

  tmp <- full_join(vtam_df, vtamR_df, 
                   by = c("read_count_cutoff", "variant_cutoff"), 
                   suffix=c("_vtam", "_vtamR")) %>%
    mutate(diff_TP = TP_vtamR - TP_vtam) %>%
    mutate(diff_FP = FP_vtamR - FP_vtam) %>%
    mutate(dataset = dataset) %>%
    select(dataset, everything())
  
  lfn_read_count_and_lfn_variant <- rbind(lfn_read_count_and_lfn_variant, tmp)
  
  tmp <- tmp %>%
      filter(diff_FP > 0)
    
    msg <- paste0("optimize_lfn_read_count_and_lfn_variant: ", ident)
    print(msg)
    print(tmp)
}
write.csv(lfn_read_count_and_lfn_variant, file = file.path(outdir, "lfn_read_count_and_lfn_variant_cutoffs.csv"), row.names = FALSE)

### variant_cutoffs in VTAM are determined by generating 10 values between the 
# minimum and maximum. This can lead to rounding of the cutoff values, which may 
# not exactly match the cutoffs used in vtamR, where values are determined based 
# the minimum, maximum, and increment.
# This can lead to slight differences in the number of false positives and true positives.