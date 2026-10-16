args <- commandArgs(trailingOnly = TRUE)
vsearch_path <- as.character(args[1])
outdir_pipe <- as.character(args[2])
log_file <- as.character(args[3])

library(vtamR)
library(dplyr)
library(tidyr)

# setwd("~/vtamR_benchmark")

#### Rational #### 
# The results of filtering are slightly different between VTAM and vtamR =>
# use the filtered results of VTAM for pooling markers both for VTAM pool function, ans vtamR pool_markers function


####  Input: filtered files by VTAM in large format ######
mfzr_filtered_vtam_file <- file.path(outdir_pipe, "fish_mfzr/vtamR/filter_default/12_pool_replicates.csv")
#mfzr_filtered_vtam <- read.table(mfzr_filtered_vtam_file, header = TRUE, sep=",")
zfzr_filtered_vtam_file <- file.path(outdir_pipe, "fish_zfzr/vtamR/filter_default/12_pool_replicates.csv")
#zfzr_filtered_vtam <- read.table(zfzr_filtered_vtam_file, header = TRUE, sep=",")




#### Use vtamR to pool the results of the 2 markers obtained by VTAM #####
outfile <- file.path(outdir_pipe, "fish_pooled_markers", "vtamR", "asvtable_pooled_markers.csv") 
asv_with_centroids <- file.path(
  "fish_pooled_markers", 
  "vtamR", 
  "asvtable_pooled_markers_with_centroids.csv"
)

vtamr_pooled <- pool_markers(
  files = c(mfzr_filtered_vtam_file, zfzr_filtered_vtam_file),
  outfile = outfile,  # ASVs identical over their overlapping region are merged into a single line
  asv_with_centroids = asv_with_centroids,  # same as input, with added centroid_id and centroid
  method = "sum",  # mean read counts across pooled ASVs
  vsearch_path = vsearch_path,
  quiet = TRUE,
  log_file = log_file
)
