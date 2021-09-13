PIPELINE_PREFIX=$PWD/$(dirname "$0")/

PIPELINE_FILE=$PIPELINE_PREFIX"/../workflow/snakefile_preprocess"
CONFIG_FILE=$PIPELINE_PREFIX"/../config/config.yaml"

snakemake --snakefile "$PIPELINE_FILE" --cores 32 -p --configfile "$CONFIG_FILE" "$@"




