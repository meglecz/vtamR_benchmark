
library(vtamR)
library(dplyr)
library(tidyr)

setwd("~/vtamR_benchmark")

######################################################################
compare_rc <- function(step, datasets, dir){
  
rc <- data.frame(
  filename = character(),
  read_count_vtam = numeric(),
  read_count_vtamr = numeric(),
  dataset = character(),
  step = character()
)

for(dataset in datasets){
  
  rc_vtam <- count_reads_in_dir(
    dir=file.path(dir, dataset, "VTAM", step), 
    pattern="\\.fasta.gz$", 
    file_type="fasta"
  )
  
  if(step== "demultiplex_trim"){
    sortedinfo = file.path(dir, dataset, "VTAM", "demultiplex_trim", "sortedinfo.tsv")
    sortedinfo_df <- read.table(sortedinfo, header = TRUE, sep="\t") %>%
      rename(filename = sortedfasta) %>%
      mutate(filename_clean = paste0(sample, "-", replicate, ".fasta.gz")) %>%
      select(filename, filename_clean)
    
    rc_vtam <- rc_vtam %>%
      left_join(sortedinfo_df, by="filename") %>%
      select(-filename) %>%
      rename(filename = filename_clean) %>%
      select(filename, read_count)
  }
  
  rc_vtamr <- count_reads_in_dir(
    dir=file.path(dir, dataset, "vtamR", step), 
    pattern="\\.fasta.gz", 
    file_type="fasta"
  )
  
  tmp <- full_join(
    rc_vtam, 
    rc_vtamr, 
    by="filename", 
    suffix = 
      c("_vtam", "_vtamr")
    ) %>%
    mutate("dataset" = dataset,
           "step" = step)

  rc <- rbind(rc, tmp)
}

rc <- rc %>%
  mutate(difference = read_count_vtam - read_count_vtamr)
return(rc)
}
######################################################################

datasets = c("fish_mfzr", "fish_zfzr", "shark") # bat were preprocessed only by vtamR
# count the number of sequences in each merged fasta file
merge_rc <- compare_rc(step="merge", datasets, "pipeline_output_1")
demultiplex_rc <- compare_rc(step="demultiplex_trim", datasets, "pipeline_output_1")

rc <- rbind(merge_rc, demultiplex_rc)
outdir <- file.path("comparaison_results", "preprocess")
outdir <- check_dir(outdir)
write.csv(rc, file = file.path(outdir, "read_count_merge_demultiplex.csv"), row.names = FALSE)


