#!/bin/bash

cd "$HOME/vtamR_benchmark/"

/usr/bin/time -f \
"vtamR, concatenate_files, NA, bat, %E, %U, %S" \
-a -o runtime_concatenate_shark.csv \
Rscript scripts/vtamR_concatenate_files.R \
> /dev/null
