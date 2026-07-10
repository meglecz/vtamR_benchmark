
library(vtamR)
library(dplyr)
library(tidyr)

setwd("~/vtamR_benchmark")

vtam <- read.csv("runtime/time_VTAM.csv", header = FALSE, strip.white = TRUE)
vtamr <- read.csv("runtime/time_vtamR.csv", header = FALSE, strip.white = TRUE)
colnames(vtam) <- c("pipeline", "run",  "step", "dataset", "elapsed", "user", "system")
colnames(vtamr) <- c("pipeline",  "run","step", "dataset", "elapsed", "user", "system")

############ transform elapsed to seconds, calculate total_cpu and average_cpu
to_seconds <- function(x) {
  sapply(strsplit(x, ":"), function(z)
    as.numeric(z[1]) * 60 + as.numeric(z[2]))
}

vtam <- vtam %>%
  mutate(elapsed = to_seconds(elapsed)) %>%
  mutate(total_cpu = user + system) %>%
  mutate(average_cpu = total_cpu / elapsed)

vtamr <- vtamr %>%
  mutate(elapsed = to_seconds(elapsed)) %>%
  mutate(total_cpu = user + system) %>%
  mutate(average_cpu = total_cpu / elapsed) %>%
  filter(!(dataset == "shark" & step == "taxassign"))# delete taxaassign shark, since it was done only by vtamR

############"
# sum optimize steps of vtam to "suggest_parameters" to match steps between vtam and vtamR
runtime_cols <- c("elapsed", "user", "system", "total_cpu", "average_cpu")
suggest_steps <- c( "make_known_occurrences", "OptimizePCRerror", "OptimizeLFNsampleReplicate", "OptimizeLFNreadCountAndLFNvariant")
optimise <- vtam %>%
  filter(step %in% suggest_steps) %>%
  group_by(dataset, run) %>%
  summarise(
    pipeline = first(pipeline),
    step = "suggest_parameters",
    across(all_of(runtime_cols), sum),
    .groups = "drop"
  )

vtam <- vtam %>%
  filter(!step %in% suggest_steps) %>%
  rbind(optimise)

compare <- full_join(vtam, vtamr, by=c("step", "dataset"), suffix= c("_vtam", "_vtamR"))

compare <- compare %>%
  mutate(prop_elapsed = elapsed_vtamR/elapsed_vtam) %>%
  mutate(prop_total_cpu = total_cpu_vtamR / total_cpu_vtam)

sum(compare$elapsed_vtamR)
sum(compare$elapsed_vtam)
sum(compare$elapsed_vtam) / sum(compare$elapsed_vtamR)

outfile <- file.path("runtime", "runtime_comparaison.csv")
write.csv(compare, file=outfile, row.names = FALSE)


