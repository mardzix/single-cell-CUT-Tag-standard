PIPELINE_FILE="/data/proj/GCB_MB/standard_scCT_preprocess/workflow/snakefile_preprocess"
CONFIG_FILE="/data/proj/GCB_MB/standard_scCT_preprocess/config/config.yaml"

snakemake --snakefile $PIPELINE_FILE  --cores 32 -p --configfile $CONFIG_FILE $@


