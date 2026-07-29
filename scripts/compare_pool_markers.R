library(vtamR)
library(tidyr)
library(dplyr)

setwd("~/vtamR_benchmark")
dir = "pipeline_output_1"

##### Results of vtamR pool_markers ####
vtamr_pooled_file <- file.path(dir, "fish_pooled_markers", "vtamR", "asvtable_pooled_markers.csv")
vtamr_pooled <- read.csv(vtamr_pooled_file, header = TRUE)

##### Results of VTAM pool --readcounts ####
vtam_pooled_file <- file.path(dir, "fish_pooled_markers", "VTAM", "asvtable_pooled_markers.tsv")
vtam_pooled <- read.table(vtam_pooled_file, header = TRUE, sep="\t")

vtam_pooled <- vtam_pooled %>%
  select(-c(run, marker, clusterid, clustersize, pooled_sequences)) %>%
  rename(asv_id = variant_id, asv = sequence) %>%
  pivot_longer(
    cols = -c("asv_id", "asv", "pooled_variants" ),
    names_to = "sample", 
    values_to = "read_count") %>%
  filter(read_count > 0) %>%
  mutate(sample = sub("^fish.", "", sample)) %>%
  select(asv_id, pooled_variants, sample, read_count, asv)

#### Compare results of pool_markers (vtamR) and pool (VTAM) #### 
comparaison <- full_join(vtam_pooled, vtamr_pooled, by=c("asv_id", "asv", "sample"), suffix = c("_vtam", "_vtamR")) %>%
  mutate(diff_read_count = read_count_vtam - read_count_vtamR) %>%
  select(asv_id, pooled_variants, sample, read_count_vtam, read_count_vtamR, diff_read_count, asv)

outdir <- file.path("comparaison_results", "pool_markers")
outdir <- check_dir(outdir)
outfile <- file.path(outdir, "pool_markers_read_counts.csv" )
write.csv(comparaison, file=outfile, row.names = FALSE)

diff <- comparaison %>%
  filter(diff_read_count != 0)
print("Occurrences with different read counts:")
print(diff)

