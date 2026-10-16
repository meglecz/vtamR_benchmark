#!/bin/bash

cd "$HOME/vtamR_benchmark/"

log_file="$HOME/vtamR_benchmark/data/fastq/preprocess_log.csv"

/usr/bin/time -f \
"vtamR, concatenate_files, NA, bat, %E, %U, %S" \
-a -o runtime_concatenate_shark.csv \
Rscript scripts/vtamR_concatenate_files.R $log_file\
> /dev/null
