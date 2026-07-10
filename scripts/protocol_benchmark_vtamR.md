# Benchmark vtamR
Emese Meglécz 
2026/05/18

## Input

Working directory: **~/vtamR_benchmark/**

**bat**

Download fastq files from https://datadryad.org/stash/dataset/doi:10.5061/dryad.kv02g, untar move fastq files to **data/fastq/bat**

The fastq files have already been demultiplexed, but primers should be trimmed.
vtamR  can handle this case, but not VTAM. 
Use vtamR to prepare the files for filtering and place them to **bat/vtamR/demultiplex_trim** and reformat the sampleinfo.csv to sortedinfo.tsv for VTAM filtering

~~~
bash ./scripts/preprocess_bat.sh
~~~

**fish**

Download fastq files from https://datadryad.org/stash/dataset/doi:10.5061/dryad.f40v5

Untar and move fastq files to **data/fastq/fish_mfzr** and **data/fastq/fish_zfzr**
Can be used as is for both VTAM and vtamR

**shark**

Download fastq files from https://osf.io/3txv7/overview, untar move fastq files to 
**data/fastq/shark/200309** and **data/fastq/shark/200513**
The files of the two sequencing runs should be pooled into one. This can be done by vtamR. The output can be used for both VTAM and vtamR.

~~~
bash ./scripts/vtamR_concatenate_files.sh
~~~

**metadata**

Fastqinfo and mock_composition, sample_types and parameter files for VTAM are in data/metadata: 

- vtam_fastqinfo_DATASET.tsv
- vtam_mock_composition_DATASET.tsv
- vtamR_fastqinfo_DATASET.csv
- vtamR_mock_composition_DATASET.csv
- sample_types_DATASET.tsv
- params_STEP_DATASET.yml

## VTAM

Run VTAM for 
- merge (skip for bat, use vtamR)
- sortreads (skip for bat, use vtamR)
- filter default 
  - sorteddir, sortedinfo from vtamR for bat
  - lfn_variant_replicate for fish and shark
  - skip FilterIndel, FilterCodonStop for shark 
- make_known_occurrences
- OptimizePCRerror
  - sorteddir, sortedinfo from vtamR for bat
- OptimizeLFNsampleReplicate
  - sorteddir, sortedinfo from vtamR for bat
- OptimizeLFNreadCountAndLFNvariant
  - sorteddir, sortedinfo from vtamR for bat
  - lfn_variant_replicate for fish and shark
  - use optimized value for lfn_sample_replicate_cutoff and pcr_error_var_prop
- taxassign
  - skip for shark
- pool fish_mfzr and fish_zfzr results


~~~
bash ./scripts/vtam.sh
~~~



## vtamR

**asv_list**
Extract complete list of asv (variant) and asv_id (id) from VTAM database (db.sqlite) and use it in vtamR to get **coherent asv_id across pipelines**.
~~~
mkdir output_info
sqlite3 -header -csv db.sqlite "SELECT * FROM Variant;" > output_info/asv_all_datasets.csv
~~~
**vtamR functions**
Run the following  vtamR functions using the order of the filters and parametres that match as much as possible the VTAM pipeline:

- merge_fastq_pairs 
- demultiplex_and_trim (for bat use trim_primers)
- filter matching default maramaters of VTAM 
  - dereplicate (VTAM: included in filter )
  - filter_asv_global (VTAM: filter global_read_count_cutoff )
  - filter_occurrence_sample (VTAM: LFN_sample_replicate)
  - filter_occurrence_variant (VTAM: LFN_variant_replicate for fish and shark LFN_variant for bat) 
  - filter_occurrence_read_count (VTAM: LFN_read_count)
  - pool_filters (VTAM: Pool LFN)
  - filter_min_replicate  (VTAM: FilterMinReplicateNumber)
  - filter_pcr_error (VTAM: FilterPCRerror)
  - filter_chimera  (VTAM: FilterChimera)
  - filter_min_replicate - 2  (VTAM: FilterMinReplicateNumber)
  - filter_replicate (VTAM: FilterRenkonen)
  - filter_min_replicate - 3 (VTAM: FilterMinReplicateNumber)
  - filter_indel (VTAM: FilterIndel; skip for shark)
  - filter_stop_codon (VTAM: FilterCodonStop; skip for shark)
  - pool_replicates (VTAM: ReadCountAverageOverReplicates)
  - write_asv_table (VTAM: MakeAsvTable)
  - classify_control_occurrences (VTAM: NA)
- suggest parametres
  - suggest_sample_cutoff (on unfiltered data; VTAM optimize --until OptimizeLFNsampleReplicate)
  - suggest_pcr_error_cutoff (on unfiltered data; VTAM optimize --until OptimizePCRerror)
  - prefilter with optimized falues
    - filter_occurrence_sample  (on unfiltered data)
    - filter_pcr_error  (on unfiltered data)
    - pool_filters
    - filter_min_replicate
  - classify_control_occurrences (on prefiltered data; VTAM: make_known_occurrences)
  - suggest_variant_readcount_cutoffs (on prefiltered data; VTAM optimize --until OptimizeLFNreadCountAndLFNvariant)
- taxassign
  - assign_taxonomy_rdp for shark
  - assign_taxonomy_ltg for the other datasets (VTAM: taxassign)

```
bash ./scripts/vtamR.sh
```



VTAM pool takes the sum of the read count of the markers

## Compare VTAM and vtamR output

### merge, demultiplex

~~~
Rscript scripts/compare_read_counts.R
~~~

Output:  **output_info/read_count_merge_demultiplex.csv**

VTAM and vtamR merge/merge_and_trim and sortreads/demultiplex_and_trim produce files with the same number of sequences. Only one exception with only one read difference:

| filename        | read_count_vtam | read_count_vtamr | dataset   | step             |
| --------------- | --------------- | ---------------- | --------- | ---------------- |
| 14Cro13-1.fasta | 41538           | 41539            | fish_mfzr | demultiplex_trim |

### filter

**Make tables for each VTAM filter from db.sqlite**

~~~
mkdir -p output_info/db_tables

sqlite3 -header -csv db.sqlite "SELECT * FROM Sample;" > output_info/db_tables/samples.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM Run;" > output_info/db_tables/runs.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM VariantReadCount;" > output_info/db_tables/VariantReadCount.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterLFN;" > output_info/db_tables/FilterLFN.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterMinReplicateNumber;" > output_info/db_tables/FilterMinReplicateNumber.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterPCRerror;" > output_info/db_tables/FilterPCRerror.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterChimera;" > output_info/db_tables/FilterChimera.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterMinReplicateNumber2;" > output_info/db_tables/FilterMinReplicateNumber2.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterRenkonen;" > output_info/db_tables/FilterRenkonen.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterMinReplicateNumber3;" > output_info/db_tables/FilterMinReplicateNumber3.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterIndel;" > output_info/db_tables/FilterIndel.csv

sqlite3 -header -csv db.sqlite "SELECT * FROM FilterCodonStop;" > output_info/db_tables/FilterCodonStop.csv
~~~

**compare_filter.R**

Compares the final output files produced by VTAM and vtamR and identifies occurrences that are present in one output but absent from the other.

For these occurrences, the script generates files reporting their presence or absence at each filtering step of the pipeline in which they occur:

- occurrences_missing_in_vtam_DATASET.csv
- occurrences_missing_in_vtamR_DATASET.csv

It also generates a file containing the read counts from both pipelines for each occurrence detected in the final filtered outputs:

- read_count_VTAM_vtamR_DATASET.csv

~~~
Rscript scripts/compare_filter.R
~~~

**Conclusions**

There are only a small number of occurrences that differ between the two pipelines. All discrepancies originate from the *filter_chimera* or *filter_pcr_error* steps.

For PCR error and chimera filtering, vtamR removes a variant when it is classified as a PCR error in more than 80% of the samples in which it occurs. In contrast, VTAM processes samples independently and makes filtering decisions on a sample-by-sample basis, without enforcing a consistent decision for a given ASV across all samples.

As a result, when analyses are performed independently for each sample, obtaining exactly identical results between the two pipelines is not possible. However, we believe that the behavior implemented in vtamR is more biologically justified, since a given ASV is unlikely to represent an artefact in one sample while being a true biological variant in another.

### Taxassign

**compare_taxassig.R**

The poprortion of different assignment between vtam and vtamR id low (2-3 %). Almost all assignements that differ betwwen the pielines are at <=80% identity, which are doubious anyway.

| dataset   | total_asv | different | prop               |
| --------- | --------- | --------- | ------------------ |
| fish_mfzr | 153       | 5         | 0.0326797385620915 |
| fish_zfzr | 105       | 2         | 0.019047619047619  |
| bat       | 1676      | 41        | 0.0244630071599045 |

### pool markers

**compare_pool_markers.R**

VTAM pool function and vtamR pool_markers funtion give identical results.

### Runtime

**compare_runtime.R**

For almost all steps and dataset vtamR is faster, both regarding elapsed time and total CPU time.

Ovarall the whole pipline takes 4484.65 seconds for vtamR and 6510.25 for VTAM, thus vtamR increses the speed by 30%.