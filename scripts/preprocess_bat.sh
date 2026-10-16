#!/bin/bash

cd "$HOME/vtamR_benchmark/"

vsearch_path="$HOME/miniconda3/envs/vtam/bin/vsearch"
cutadapt_path="$HOME/miniconda3/envs/vtam_benchmark/bin/cutadapt"
pigz_path="$HOME/miniconda3/envs/vtam/bin/pigz"
log_file="$HOME/vtamR_benchmark/data/fastq/preprocess_log.csv"


/usr/bin/time -f \
"vtamR, merge, NA, bat, %E, %U, %S" \
-a -o runtime_preprocess_bat.csv \
Rscript \
    scripts/vtamR_merge.R bat $vsearch_path $pigz_path $log_file\
> /dev/null


/usr/bin/time -f \
"vtamR, trim_primers, NA, bat, %E, %U, %S" \
-a -o runtime_preprocess_bat.csv \
Rscript \
    scripts/vtamR_trim_primers.R bat $vsearch_path $cutadapt_path $pigz_path $log_file\
> /dev/null

