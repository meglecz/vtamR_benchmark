args <- commandArgs(trailingOnly = TRUE)
vsearch_path <- as.character(args[1])
outdir_pipe <- as.character(args[2])

library(vtamR)
library(dplyr)
library(tidyr)

# setwd("~/vtamR_benchmark")

#### Rational #### 
# The results of filtering are slightly different between VTAM and vtamR =>
# use the filtered results of VTAM for pooling markers both for VTAM pool function, ans vtamR pool_markers function


####  Input: filtered files by VTAM in large format ######
mfzr_filtered_vtam_file <- file.path(outdir_pipe, "fish_mfzr/VTAM/filter_default/asvtable_default.tsv")
mfzr_filtered_vtam <- read.table(mfzr_filtered_vtam_file, header = TRUE, sep="\t")
zfzr_filtered_vtam_file <- file.path(outdir_pipe, "fish_zfzr/VTAM/filter_default/asvtable_default.tsv")
zfzr_filtered_vtam <- read.table(zfzr_filtered_vtam_file, header = TRUE, sep="\t")

# Make long format ----

mfzr_filtered_vtam <- mfzr_filtered_vtam %>%
  select(-c(run, marker, sequence_length, read_count, clusterid, clustersize, chimera_borderline)) %>%
  rename(asv_id = variant, asv = sequence) %>%
  pivot_longer(
    cols = -c("asv_id", "asv"),
    names_to = "sample", 
    values_to = "read_count") %>%
  filter(read_count > 0) %>%
  mutate(sample = sub("^X", "", sample)) %>%
  select(asv_id, sample, read_count, asv)

zfzr_filtered_vtam <- zfzr_filtered_vtam %>%
  select(-c(run, marker, sequence_length, read_count, clusterid, clustersize, chimera_borderline)) %>%
  rename(asv_id = variant, asv = sequence) %>%
  pivot_longer(
    cols = -c("asv_id", "asv"),
    names_to = "sample", 
    values_to = "read_count") %>%
  filter(read_count > 0) %>%
  mutate(sample = sub("^X", "", sample)) %>%
  select(asv_id, sample, read_count, asv)

# Write files, since pool_markers does not accept df ----
mfzr_filtered <- file.path(outdir_pipe, "fish_mfzr", "VTAM", "filter_default","asvtable_default_long.csv" )
write.csv(mfzr_filtered_vtam, file = mfzr_filtered, row.names = FALSE)
zfzr_filtered <- file.path(outdir_pipe, "fish_zfzr", "VTAM", "filter_default","asvtable_default_long.csv" )
write.csv(zfzr_filtered_vtam, file = zfzr_filtered, row.names = FALSE)


#### Use vtamR to pool the results of the 2 markers obtained by VTAM #####
outfile <- file.path(outdir_pipe, "fish_pooled_markers", "vtamR", "asvtable_pooled_markers.csv") 
asv_with_centroids <- file.path(
  "fish_pooled_markers", 
  "vtamR", 
  "asvtable_pooled_markers_with_centroids.csv"
)

vtamr_pooled <- pool_markers(
  files = c(mfzr_filtered, zfzr_filtered),
  outfile = outfile,  # ASVs identical over their overlapping region are merged into a single line
  asv_with_centroids = asv_with_centroids,  # same as input, with added centroid_id and centroid
  method = "sum",  # mean read counts across pooled ASVs
  vsearch_path = vsearch_path,
  quiet = TRUE
)
