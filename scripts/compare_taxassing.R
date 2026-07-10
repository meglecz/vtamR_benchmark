
library(vtamR)
library(dplyr)
library(tidyr)

setwd("~/vtamR_benchmark")

datasets =  c("fish_mfzr", "fish_zfzr", "bat") # shark is only analyzed by vtamR
dir = "pipeline_output_1"
outdir <- file.path("comparaison_results", "compare_taxassign")
outdir = check_dir(path=outdir, is_file=FALSE)

count <- data.frame(
  dataset = character(),
  total_asv = integer(),
  different = integer()
)

for(dataset in datasets){
  
  vtam_file <- file.path(dir, dataset, "VTAM", "taxassign", "asvtable_default_taxa.tsv")
  vtam_df <- read.table(vtam_file, sep="\t", header = TRUE, na.strings = "") %>%
    select(asv_id=variant, 
           ltg_taxid = ltg_tax_id, 
           ltg_name = ltg_tax_name,
           ltg_rank,
           pid = identity,
           phylum, class, order, family, genus, species, 
           asv = sequence
           )
  
  vtamr_file <- file.path(dir, dataset, "vtamR", "taxassign", "asv_default_taxa.csv")
  vtamr_df <- read.csv(vtamr_file) %>%
    select(asv_id, ltg_taxid, ltg_name, ltg_rank, pid,
           phylum, class, order, family, genus, species, asv)
  
  pool <- full_join(vtam_df, vtamr_df, by=c("asv_id", "asv"), suffix= c("_vtam", "_vtamR")) %>%
    mutate(dataset = dataset) %>%
    select(dataset, asv_id, ltg_taxid_vtam, ltg_taxid_vtamR, ltg_name_vtam, ltg_name_vtamR, ltg_rank_vtam, ltg_rank_vtamR, pid_vtam, pid_vtamR,
           phylum_vtam, phylum_vtamR, class_vtam, class_vtamR, order_vtam, order_vtamR, family_vtam, family_vtamR, genus_vtam, genus_vtamR,
           species_vtam, species_vtamR, asv )
  
  pooled_outfile <- paste0("VTAM_vtamR_taxassign_pooled_", dataset, ".csv")
  pooled_outfile <- file.path(outdir, pooled_outfile)
  write.csv(pool, file=pooled_outfile, row.names = FALSE)
  
  # select ASV where the assignment is different at the phylum, class, order, family, genus or species level
  diff <-  pool %>%
    filter(!(is.na(phylum_vtam) & is.na(phylum_vtamR))) %>%
    filter(phylum_vtam != phylum_vtamR | is.na(phylum_vtam) | is.na(phylum_vtamR))
  
  diff_class <-  pool %>%
    filter( !asv_id %in% diff$asv_id) %>%
    filter(!(is.na(class_vtam) & is.na(class_vtamR))) %>%
    filter(class_vtam != class_vtamR | is.na(class_vtam) | is.na(class_vtamR))
  diff <- rbind(diff, diff_class)
  
  diff_order <-  pool %>%
    filter( !asv_id %in% diff$asv_id) %>%
    filter(!(is.na(order_vtam) & is.na(order_vtamR))) %>%
    filter(order_vtam != order_vtamR | is.na(order_vtam) | is.na(order_vtamR))
  diff <- rbind(diff, diff_order)
  
  diff_family <-  pool %>%
    filter( !asv_id %in% diff$asv_id) %>%
    filter(!(is.na(family_vtam) & is.na(family_vtamR))) %>%
    filter(family_vtam != family_vtamR | is.na(family_vtam) | is.na(family_vtamR))
  diff <- rbind(diff, diff_family)
  
  diff_genus <-  pool %>%
    filter( !asv_id %in% diff$asv_id) %>%
    filter(!(is.na(genus_vtam) & is.na(genus_vtamR))) %>%
    filter(genus_vtam != genus_vtamR | is.na(genus_vtam) | is.na(genus_vtamR))
  diff <- rbind(diff, diff_genus)  
  
  diff_species <-  pool %>%
    filter( !asv_id %in% diff$asv_id) %>%
    filter(!(is.na(species_vtam) & is.na(species_vtamR))) %>%
    filter(species_vtam != species_vtamR | is.na(species_vtam) | is.na(species_vtamR))
  diff <- rbind(diff, diff_species)  
  
  diff_file = paste("VTAM_vtamR_different_",  dataset, ".csv", sep="")
  diff_file = file.path(outdir, diff_file)
  write.csv(diff, file=diff_file, row.names=FALSE)
  
  count_tmp <- data.frame(
    dataset = dataset,
    total_asv = nrow(pool),
    different = nrow(diff)
  )
  count <- rbind(count, count_tmp)
}

count <- count %>%
  mutate(prop =  different / total_asv)
count_file = file.path(outdir, "summary.csv")
write.csv(count, file=count_file, row.names=FALSE)

