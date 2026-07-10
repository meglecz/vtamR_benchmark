args <- commandArgs(trailingOnly = TRUE)
dataset <- as.character(args[1])
blast_path <- as.character(args[2])
taxonomy <- as.character(args[3])
blast_db <- as.character(args[4])
outdir_pipe <- as.character(args[5])

library(vtamR)
library(dplyr)

# setwd("~/vtamR_benchmark")

# use VTAM output to be able to compare runtimes (same input for VTAM and vtamR taxassing)
asv <- file.path(outdir_pipe, dataset, "VTAM", "filter_default", "asvtable_default.tsv")
outfile <- file.path(outdir_pipe, dataset, "vtamR", "taxassign", "asv_default_taxa.csv")

# select only sequences are sequence ids and rename them to be compatible with vtamR
asv_df <- read.table(asv, header = TRUE, sep="\t") %>%
  select(
    asv_id = variant,
    asv = sequence
  )

#parametres that match as much as possible VTAM taxassign
ltg_params_df = data.frame( pid=c(100,99,97,95,90,85,80,75),
                              pcov=c(80,80,80,80,80,80,80,80),
                              phit=c(90,90,90,90,90,90,90,90),
                              taxn=c(1,1,1,3,3,3,3,3),
                              seqn=c(1,1,1,3,3,3,3,3),
                              refres=c(0,0,0,0,0,0,0,0),
                              ltgres=c(8,8,8,8,8,8,8,8)
                              )

if(dataset == "shark"){ # use rdp classifier for 16S
  taxa <- assign_taxonomy_rdp(
    asv=asv_df,
    max_memory = 8,
    confidence = 0.8,
    rm_chloroplast = TRUE,
    outfile = outfile
  )
}else{
  asv_tax <- assign_taxonomy_ltg(
    asv=asv_df,
    taxonomy=taxonomy, 
    blast_db=blast_db, 
    blast_path=blast_path,
    outfile=outfile,
    ltg_params = ltg_params_df
    )
}
