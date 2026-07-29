#!/bin/bash
set -euo pipefail

cd "$HOME/vtamR_benchmark/"

vsearch_path="$HOME/miniconda3/envs/vtam/bin/vsearch"
cutadapt_path="$HOME/miniconda3/envs/vtam_benchmark/bin/cutadapt"
pigz_path="$HOME/miniconda3/envs/vtam/bin/pigz"
blast_path="$HOME/miniconda3/envs/vtam/bin/blastn"
blast_db="$HOME/mkCOInr/COInr/COInr_for_vtam_2026_05_03_dbV5/COInr_for_vtam"
taxonomy="$HOME/mkCOInr/COInr/COInr_for_vtam_2026_05_03_dbV5/COInr_for_vtam_taxonomy.tsv"

datasets=("bat" "fish_mfzr" "fish_zfzr" "shark")

for i in {1..10}; do
    echo "Iteration $i"
    
    outdir="pipeline_output_$i"
    
#################################################
# STEP Extract complete list of asv asv_id from VTAM database
#################################################

    sqlite_db="$outdir/db.sqlite"
    asv_id="$outdir/asv_all_datasets.csv"
    sqlite3 -header -csv $sqlite_db "SELECT * FROM Variant;" > $asv_id


    for ds in "${datasets[@]}"; do

        echo "Processing $ds"

    #################################################
    # STEP merge : skip only for bat
    #################################################

        if [[ "$ds" != "bat" ]]; then

            /usr/bin/time -f \
           "vtamR, $i, merge, ${ds}, %E, %U, %S" \
            -a -o runtime/time_vtamR.csv \
            Rscript scripts/vtamR_merge.R \
            ${ds} $vsearch_path $pigz_path $outdir\
            > /dev/null
        
        else
            echo "Skipping merge for $ds"
       fi
    


    #################################################
    # STEP demultiplex : skip only for bat
    #################################################

        if [[ "$ds" != "bat" ]]; then

            /usr/bin/time -f \
            "vtamR, $i, demultiplex_trim, ${ds}, %E, %U, %S" \
            -a -o runtime/time_vtamR.csv \
            Rscript scripts/vtamR_demultiplex.R \
                ${ds} $vsearch_path $cutadapt_path $pigz_path $outdir\
            > /dev/null
        
        else
            echo "Skipping demultiplex_trim for $ds"
        fi


    
    #################################################
    # STEP filter_default
    #################################################

        pcr_error_var_prop="0.1"
        sample_cutoff="0.001"
        read_count_cutoff="10"
        variant_cutoff="0.001"
    
        if [[ "$ds" == "bat" ]]; then
    
        by_replicate="FALSE"
        coding="TRUE"
    
        elif [[ "$ds" == "shark" ]]; then
    
        by_replicate="TRUE"
        coding="FALSE"
    
        else
    
        by_replicate="TRUE"
        coding="TRUE"
    
        fi
        
        /usr/bin/time -f \
        "vtamR, $i, filter_default, ${ds}, %E, %U, %S" \
        -a -o runtime/time_vtamR.csv \
        Rscript scripts/vtamR_filter.R \
            ${ds} filter_default $vsearch_path $cutadapt_path $pcr_error_var_prop $sample_cutoff $read_count_cutoff $variant_cutoff $by_replicate $coding $outdir\
        > /dev/null

    
    #################################################
    # STEP suggest_parameters
    #################################################

    
        if [[ "$ds" == "bat" ]]; then
    
        pcr_error_var_prop="0.2"
        sample_cutoff="0.002"
        by_replicate="FALSE"
    
        elif [[ "$ds" == "shark" ]]; then
    
        pcr_error_var_prop="0.2"
        sample_cutoff="0.001"
        by_replicate="TRUE"
    
        else
    
        pcr_error_var_prop="0.1"
        sample_cutoff="0.003"
        by_replicate="TRUE"
    
        fi

        /usr/bin/time -f \
        "vtamR, $i, suggest_parameters, ${ds}, %E, %U, %S" \
        -a -o runtime/time_vtamR.csv \
        Rscript scripts/vtamR_suggest_parameters.R \
            ${ds} $vsearch_path $pcr_error_var_prop $sample_cutoff $by_replicate $outdir\
        > /dev/null

    #################################################
    # STEP taxassign - RDP for shark, LTG for the others
    #################################################
    
        /usr/bin/time -f \
        "vtamR, $i, taxassign, ${ds}, %E, %U, %S" \
        -a -o runtime/time_vtamR.csv \
        Rscript scripts/vtamR_taxassign.R \
            ${ds} $blast_path $taxonomy $blast_db $outdir\
        > /dev/null

    done

    #################################################
    # STEP pool mfzr and zfzr markers for fish
    #################################################
    
    mkdir -p $outdir/fish_pooled_markers/vtamR
    
    /usr/bin/time -f \
    "vtamR, $i, pool_markers, fish, %E, %U, %S" \
    -a -o runtime/time_vtamR.csv \
    Rscript scripts/vtamR_pool_markers.R \
        $vsearch_path $outdir \
    > /dev/null
    
done

