#!/bin/bash

set -euo pipefail
set -x

cd "$HOME/vtamR_benchmark/"

db_dir="$HOME/mkCOInr/COInr/COInr_for_vtam_2026_05_03_dbV5"
taxonomy="$HOME/mkCOInr/COInr/COInr_for_vtam_2026_05_03_dbV5/COInr_for_vtam_taxonomy.tsv"

datasets=("bat" "fish_mfzr" "fish_zfzr" "shark")

for i in {9..10}; do
    echo "Iteration $i"

    outdir="pipeline_output_$i"
    mkdir -p "$outdir"
       
    for ds in "${datasets[@]}"; do

        echo "Processing $ds"

    #################################################
    # STEP merge : skip only for bat
    #################################################

        if [[ "$ds" != "bat" ]]; then

            /usr/bin/time -f \
            "VTAM, $i, merge, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam merge \
                --fastqinfo data/metadata/vtam_fastqinfo_${ds}.tsv \
                --fastqdir data/fastq/${ds} \
                --fastainfo $outdir/${ds}/VTAM/merge/fastainfo.tsv \
                --fastadir $outdir/${ds}/VTAM/merge \
                -v \
                > /dev/null

        else
            echo "Skipping merge for $ds"
        fi

    #################################################
    # STEP sortreads : skip only for bat
    #################################################

        if [[ "$ds" != "bat" ]]; then

            /usr/bin/time -f \
            "VTAM, $i, demultiplex_trim, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam sortreads \
                --fastainfo $outdir/${ds}/VTAM/merge/fastainfo.tsv \
                --fastadir $outdir/${ds}/VTAM/merge \
                --sorteddir $outdir/${ds}/VTAM/demultiplex_trim \
                -v \
            > /dev/null

        else
            echo "Skipping sortreads for $ds"
        fi

    #################################################
    # STEP filter_default - lfn_variant_replicate for fish and shark - skip FilterIndel, FilterCodonStop for shark
    #################################################
    
        mkdir -p ${ds}/VTAM/filter_default
        if [[ "$ds" == "bat" ]]; then

            /usr/bin/time -f \
            "VTAM, $i, filter_default, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam filter \
                --db $outdir/db.sqlite \
                --sortedinfo data/fastq/bat/demultiplex_trim/sortedinfo.tsv \
                --sorteddir data/fastq/bat/demultiplex_trim/ \
                --asvtable $outdir/${ds}/VTAM/filter_default/asvtable_default.tsv \
                -v \
            > /dev/null
        
        elif [[ "$ds" == "shark" ]]; then

            /usr/bin/time -f \
            "VTAM, $i, filter_default, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam filter \
                --db $outdir/db.sqlite \
                --sortedinfo $outdir/${ds}/VTAM/demultiplex_trim/sortedinfo.tsv \
                --sorteddir $outdir/${ds}/VTAM/demultiplex_trim \
                --asvtable $outdir/${ds}/VTAM/filter_default/asvtable_default.tsv \
                --lfn_variant_replicate \
                --params data/metadata/params_filter1_${ds}.yml \
                -v \
            > /dev/null
    
        else

            /usr/bin/time -f \
            "VTAM, $i, filter_default, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam filter \
                --db $outdir/db.sqlite \
                --sortedinfo $outdir/${ds}/VTAM/demultiplex_trim/sortedinfo.tsv \
                --sorteddir $outdir/${ds}/VTAM/demultiplex_trim \
                --asvtable $outdir/${ds}/VTAM/filter_default/asvtable_default.tsv \
                --lfn_variant_replicate \
                -v \
            > /dev/null

        fi
    
    #################################################
    # STEP make_known_occurrences
    #################################################
    
         mkdir -p ${ds}/VTAM/optimize_suggest
    
        /usr/bin/time -f \
        "VTAM, $i, make_known_occurrences, ${ds}, %E, %U, %S" \
        -a -o runtime/time_VTAM.csv \
        vtam make_known_occurrences \
            --asvtable $outdir/${ds}/VTAM/filter_default/asvtable_default.tsv \
            --sample_types data/metadata/sample_types_${ds}.tsv \
            --mock_composition data/metadata/vtam_mock_composition_${ds}.tsv \
            --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
            --missing_occurrences $outdir/${ds}/VTAM/optimize_suggest/missing_occurrences_${ds}.tsv \
            -v \
        > /dev/null

    #################################################
    # STEP OptimizePCRerror 
    #################################################
    
        if [[ "$ds" == "bat" ]]; then
    
            /usr/bin/time -f \
            "VTAM, $i, OptimizePCRerror, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam optimize \
                --db $outdir/db.sqlite  \
                --sortedinfo data/fastq/bat/demultiplex_trim/sortedinfo.tsv \
                --sorteddir data/fastq/bat/demultiplex_trim \
                --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
                --outdir $outdir/${ds}/VTAM/optimize_suggest \
                --until OptimizePCRerror \
                -v \
            > /dev/null
    
        else
    
            /usr/bin/time -f \
            "VTAM, $i, OptimizePCRerror, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam optimize \
                --db $outdir/db.sqlite  \
                --sortedinfo $outdir/${ds}/VTAM/demultiplex_trim/sortedinfo.tsv \
                --sorteddir $outdir/${ds}/VTAM/demultiplex_trim \
                --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
                --outdir $outdir/${ds}/VTAM/optimize_suggest \
                --until OptimizePCRerror \
                -v \
            > /dev/null
    
        fi
    


    #################################################
    # STEP OptimizeLFNsampleReplicate
    #################################################
    
        if [[ "$ds" == "bat" ]]; then
    
            /usr/bin/time -f \
            "VTAM, $i, OptimizeLFNsampleReplicate, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam optimize \
                --db $outdir/db.sqlite \
                --sortedinfo data/fastq/bat/demultiplex_trim/sortedinfo.tsv \
                --sorteddir data/fastq/bat/demultiplex_trim \
                --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
                --outdir $outdir/${ds}/VTAM/optimize_suggest \
                --until OptimizeLFNsampleReplicate \
                -v \
            > /dev/null
    
        else
    
        /usr/bin/time -f \
        "VTAM, $i, OptimizeLFNsampleReplicate, ${ds}, %E, %U, %S" \
        -a -o runtime/time_VTAM.csv \
        vtam optimize \
            --db $outdir/db.sqlite \
            --sortedinfo $outdir/${ds}/VTAM/demultiplex_trim/sortedinfo.tsv \
            --sorteddir $outdir/${ds}/VTAM/demultiplex_trim \
            --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
            --outdir $outdir/${ds}/VTAM/optimize_suggest \
            --until OptimizeLFNsampleReplicate \
            -v \
        > /dev/null
    
        fi
 
     #################################################
    # STEP OptimizeLFNreadCountAndLFNvariant lfn_variant_replicate for fish and shark
    #################################################

        if [[ "$ds" == "bat" ]]; then
    
        /usr/bin/time -f \
        "VTAM, $i, OptimizeLFNreadCountAndLFNvariant, ${ds}, %E, %U, %S" \
        -a -o runtime/time_VTAM.csv \
        vtam optimize \
            --db $outdir/db.sqlite \
            --sortedinfo data/fastq/bat/demultiplex_trim/sortedinfo.tsv \
            --sorteddir data/fastq/bat/demultiplex_trim \
            --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
            --outdir $outdir/${ds}/VTAM/optimize_suggest \
            --params data/metadata/params_optimize_${ds}.yml \
            --until OptimizeLFNreadCountAndLFNvariant \
            -v \
        > /dev/null


        else

        /usr/bin/time -f \
        "VTAM, $i, OptimizeLFNreadCountAndLFNvariant, ${ds}, %E, %U, %S" \
        -a -o runtime/time_VTAM.csv \
        vtam optimize \
            --db $outdir/db.sqlite \
            --sortedinfo $outdir/${ds}/VTAM/demultiplex_trim/sortedinfo.tsv \
            --sorteddir $outdir/${ds}/VTAM/demultiplex_trim \
            --known_occurrences $outdir/${ds}/VTAM/optimize_suggest/known_occurrences_${ds}.tsv \
            --outdir $outdir/${ds}/VTAM/optimize_suggest \
            --lfn_variant_replicate \
            --params data/metadata/params_optimize_${ds}.yml \
            --until OptimizeLFNreadCountAndLFNvariant \
            -v \
        > /dev/null

        fi


    #################################################
    # STEP taxassign skip for shark
    #################################################
    
        if [[ "$ds" != "shark" ]]; then
    
         mkdir -p $outdir/${ds}/VTAM/taxassign

            /usr/bin/time -f \
            "VTAM, $i, taxassign, ${ds}, %E, %U, %S" \
            -a -o runtime/time_VTAM.csv \
            vtam taxassign \
                --db $outdir/db.sqlite \
                --asvtable $outdir/${ds}/VTAM/filter_default/asvtable_default.tsv \
                --output $outdir/${ds}/VTAM/taxassign/asvtable_default_taxa.tsv \
                --taxonomy $taxonomy \
                --blastdbdir $db_dir \
                --blastdbname COInr_for_vtam \
                --mode reset \
                -v \
            > /dev/null

        else
            echo "Skipping taxassign for $ds"
        fi
    
    done
    
    #################################################
    # STEP pool fish_mfzr and fish_zfzr
    #################################################
    
     mkdir -p $outdir/fish_pooled_markers/VTAM
     
    /usr/bin/time -f \
    "VTAM, $i, pool_markers, fish, %E, %U, %S" \
        -a -o runtime/time_VTAM.csv \
        vtam pool \
        --db $outdir/db.sqlite \
        --runmarker data/metadata/pool_markers.tsv \
        --asvtable $outdir/fish_pooled_markers/VTAM/asvtable_pooled_markers.tsv \
        --readcounts \
        -v
    > /dev/null

done
