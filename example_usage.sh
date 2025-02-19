#!/usr/bin/env bash

# This script is an example of how to run Talos, and is not intended to be run as-is
# It can be run from any directory within the Talos image, or any environment where the
# Talos package is installed
#
# If this is run end-to-end in one

set -e

# set the path to use for an output directory
OUTPUT_DIR=${1:-$(date +%F)}
mkdir -p "$OUTPUT_DIR"

# pass the populated Config TOML file, and export as an environment variable
CONFIG_FILE=${2:-config.toml}
export TALOS_CONFIG="$CONFIG_FILE"

# pass the Pedigree file to the script
PED_FILE=${3:-example_ped.txt}

# pass the MatrixTable of variants to the script
VARIANT_MT=${4:-variants.mt}

# pass both ClinVar tables to the script
CLINVAR_DECISIONS=${5:-clinvar_decisions.ht}
CLINVAR_PM5=${6:-clinvar_pm5.ht}

# pass the HPO OBO file from http://purl.obolibrary.org/obo/hp.obo
HPO_OBO=${7:-hp.obo}

# [optional] pass the MatrixTable of SVs to the script
SV_MT=${8:-sv.mt}

# identify the PanelApp data to use
MATCHED_PANELS="${OUTPUT_DIR}/matched_panels.json"
GeneratePanelData \
  -i "$PED_FILE" \
  --hpo "${HPO_OBO}" \
  --out "$MATCHED_PANELS"

# query PanelApp for panels
PANELAPP_RESULTS="${OUTPUT_DIR}/panelapp_results.json"
QueryPanelapp \
  --panels "$MATCHED_PANELS" \
  --out "$PANELAPP_RESULTS"

# find the symbol~ENSG map for all genes in the panelapp results
GENE_LOOKUP="${OUTPUT_DIR}/gene_lookup.json"
FindGeneSymbolMap \
  --panelapp "$PANELAPP_RESULTS" \
  --out_path "$GENE_LOOKUP"

# run Hail filtering on the small variant MatrixTable
# this step is the most resource-intensive, so I'd recommend running it on a VM with more resources
# aim for a machine with at least 8-cores, 16GB RAM
# the config entry 'RunHailFiltering.cores.small_variants' is consulted when setting up the PySpark cluster
# with a default of 8, but you can override this if you have more cores available.
SMALL_VARIANT_VCF="${OUTPUT_DIR}/small_variants.vcf.bgz"
RunHailFiltering \
  --mt "$VARIANT_MT" \
  --panelapp "$PANELAPP_RESULTS" \
  --pedigree "$PED_FILE" \
  --vcf_out "$SMALL_VARIANT_VCF" \
  --clinvar "$CLINVAR_DECISIONS" \
  --pm5 "$CLINVAR_PM5" \
  --checkpoint small_var_checkpoint.mt

# If you have SVs, run Hail filtering on the SV MatrixTable
SV_VARIANT_VCF="${OUTPUT_DIR}/sv_variants.vcf.bgz"
if [ -n "$SV_MT" ]; then
  RunHailFilteringSV \
    --mt "$SV_MT" \
    --panelapp "$PANELAPP_RESULTS" \
    --pedigree "$PED_FILE" \
    --vcf_out "$SV_VARIANT_VCF"
fi

# run the MOI validation
MOI_RESULTS="${OUTPUT_DIR}/moi_results.json"
# If the SV MatrixTable was provided, run the SV version of the validation
if [ -n "$SV_VARIANT_VCF" ]; then
    ValidateMOI \
      --labelled_vcf "$SMALL_VARIANT_VCF" \
      --labelled_sv "$SV_VARIANT_VCF" \
      --out_json "$MOI_RESULTS" \
      --panelapp "$PANELAPP_RESULTS" \
      --pedigree "$PED_FILE" \
      --participant_panels "$MATCHED_PANELS"
else
    ValidateMOI \
      --labelled_vcf "$SMALL_VARIANT_VCF" \
      --out_json "$MOI_RESULTS" \
      --panelapp "$PANELAPP_RESULTS" \
      --pedigree "$PED_FILE" \
      --participant_panels "$MATCHED_PANELS"
fi

# run the Phenotype-based flagging
# pass the HPO OBO file from http://purl.obolibrary.org/obo/hp.obo
# downloaded from https://data.monarchinitiative.org/monarch-kg/latest/phenio.db.gz and unzipped
PHENIO_DB=PLACEHOLDER
# download Genes to Phenotypes from https://hpo.jax.org/app/data/annotations
GENE_TO_PHENOTYPE=PLACEHOLDER
PHENO_ANNOTATED_RESULTS="${OUTPUT_DIR}/phenotype_annotated_results.json"
PHENO_FILTERED_RESULTS="${OUTPUT_DIR}/phenotype_filtered_results.json"
HPOFlagging \
  --results "$MOI_RESULTS" \
  --gene_map "$GENE_LOOKUP" \
  --gen2phen "$GENE_TO_PHENOTYPE" \
  --phenio "$PHENIO_DB" \
  --out "$PHENO_ANNOTATED_RESULTS" \
  --phenout "$PHENO_FILTERED_RESULTS"

# generate the HTML report
HTML_REPORT="${OUTPUT_DIR}/talos_results.html"
LATEST_HTML_REPORT="${OUTPUT_DIR}/latest_talos_results.html"
CreateTalosHTML \
  --results "$PHENO_ANNOTATED_RESULTS" \
  --panelapp "$PANELAPP_RESULTS" \
  --output "$HTML_REPORT" \
  --latest "$LATEST_HTML_REPORT"

# generate the Seqr file
SEQR_LABELS="${OUTPUT_DIR}/seqr_labels.json"
PHENOTYPE_SPECIFIC_SEQR_LABELS="${OUTPUT_DIR}/phenotype_match_seqr_labels.json"
GenerateSeqrFile "$MOI_RESULTS" "$PHENOTYPE_SPECIFIC_SEQR_LABELS"
