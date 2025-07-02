#!/usr/bin/env bash
# Module: run_pipeline1.sh
# Author: Nancy Anderson
# Description: Entry point script for pipeline1

################################################################################
# USAGE:
# With default values:
#   bash run_pipeline1.sh [options]
#
# Example with mice genome:
#   bash run_pipeline1.sh -t 4 -r mm10
#
# OPTIONS:
#   -t THREADS      Number of threads to use (default: 4)
#   -r REFERENCE    Reference genome prefix (default: hg38)
#                   Allowed values: hg38, mm10
#   -a ADAPTER      Adapter file name to use for trimming (default: tn5_truseq)
#   -p PLATFORM     Sequencing platform (default: ILLUMINA)
#                   Valid platforms: ILLUMINA, IONTORRENT, PACBIO, ONT, SOLID, NANOPORE
#   -m MAPPING      Path to metadata mapping file (default: metadata/mapping.tsv)
#   -f FORMAT       Report format for alignment module (default: human)
#                   Allowed formats: human, csv, tsv
#
# EXAMPLE:
#   bash run_pipeline1.sh -t 8 -r mm10 -a tn5_nextera -p IONTORRENT -m data/sample_metadata.tsv -f csv
#
################################################################################


set -uo pipefail

# === SCRIPT ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}" 
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Color Definitions ===
if [[ -t 1 ]]; then
    readonly RED='\033[0;31m'
    readonly GREEN='\033[0;32m'
    readonly YELLOW='\033[0;33m'
    readonly BLUE='\033[0;34m'
    readonly NC='\033[0m'
else
    readonly RED='' GREEN='' YELLOW='' BLUE='' NC=''
fi

#=== log Color function ===

log() {
    local level="$1"
    local message="$2"
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    
    case "$level" in
        ERROR) echo -e "${RED}[${timestamp}] ERROR: ${message}${NC}" >&2 ;;
        WARN)  echo -e "${YELLOW}[${timestamp}] WARNING: ${message}${NC}" >&2 ;;
        INFO)  echo -e "${BLUE}[${timestamp}] INFO: ${message}${NC}" ;;
        *)     echo "[${timestamp}] ${message}" ;;
    esac
}

# Default values
THREADS=4
REFERENCE="hg38"
ADAPTER="tn5_truseq"
REPORT_FORMAT="${REPORT_FORMAT:-human}"
PLATFORM="ILLUMINA"
VALID_PLATFORMS=("ILLUMINA" "IONTORRENT" "PACBIO" "ONT" "SOLID" "NANOPORE")
MAPPING="metadata/mapping.tsv"  # default
# Help message
usage() {
  echo "Usage: $0 [-t THREADS] [-r REFERENCE] [-a ADAPTER] [-p PLATFORM] [-m MAPPING]"
  echo ""
  echo "Options:"
  echo "  -t THREADS    Number of threads to use (default: 4)"
  echo "  -r REFERENCE  Reference genome prefix (default: hg38). Must be 'hg38' or 'mm10'"
  echo "  -a ADAPTER    Adapter file name (default: tn5_truseq)"
  echo "  -p PLATFORM   Sequencing platform (default: ILLUMINA)"
  echo "               Valid platforms: ${VALID_PLATFORMS[*]}"
  echo "  -m MAPPING    Path to metadata mapping file (default: metadata/mapping.tsv)"
  echo ""
  echo "Example:"
  echo "  bash run_pipeline1.sh -t 8 -r mm10 -a tn5_nextera -p IONTORRENT -m data/sample_metadata.tsv"
  exit 1
}


# Help override before parsing
if [[ "${1:-}" == "--help" || "${1:-}" == "-h" ]]; then
  usage
fi

# === Valid report formats ===
VALID_FORMATS=("human" "csv" "tsv")

# === Parse command-line arguments ===
while getopts ":t:r:a:p:m:f:" opt; do
  case $opt in
    t) THREADS="$OPTARG" ;;
    r) REFERENCE="$OPTARG" ;;
    a) ADAPTER="$OPTARG" ;;
    p) PLATFORM="${OPTARG^^}" ;;  # Force uppercase
    m) MAPPING="$OPTARG" ;;
    f) REPORT_FORMAT="$OPTARG" ;;
    *) usage ;;
  esac
done


#=== Validate Platform ===
if [[ ! " ${VALID_PLATFORMS[*]} " =~ " $PLATFORM " ]]; then
  echo "Invalid platform: $PLATFORM. Must be one of: ${VALID_PLATFORMS[*]}"
  exit 1
fi

# === Validate REFERENCE ===
if [[ "$REFERENCE" != "hg38" && "$REFERENCE" != "mm10" ]]; then
  echo "Invalid reference: $REFERENCE. Must be 'hg38' or 'mm10'."
  exit 1
fi

# === Validate report format if set ===
if [[ -n "${REPORT_FORMAT:-}" ]]; then
  if [[ ! " ${VALID_FORMATS[*]} " =~ " $REPORT_FORMAT " ]]; then
    echo "❌ Invalid report format: $REPORT_FORMAT. Allowed: ${VALID_FORMATS[*]}"
    exit 1
  fi
else
  REPORT_FORMAT="human"  # default if not passed
fi

# Export environment variables used by modules
export PIPELINE_RUNNING=true
export PIPELINE_LOG_MODE=hybrid  # Change to 'centralized' if needed

log INFO "Running pipeline with THREADS=$THREADS and REFERENCE=$REFERENCE"
log INFO "Starting pipeline..."
start_time=$(date +%s)
# === Validate metadata mapping file before running pipeline ===
log INFO "🔎 Validating metadata file: $MAPPING"
python3 modules/utils/validate_mapping.py --mapping "$MAPPING" --schema templates/mapping_schema.yaml

if [[ $? -ne 0 ]]; then
  log ERROR "❌ Metadata validation failed. Exiting."
  exit 1
fi


# Step-by-step execution
log INFO " 🗃️ Creating pipeline1 directories and folders ..."
# === MODULE 01: Directories ===
modules/pipeline1/01_init_directories.sh

log INFO " 🧬 Checking genome references, blacklist, gtf and TSS bed files ..."
# === MODULE 02: Checking Reference ===
modules/pipeline1/02_reference_check.sh "Reference/${REFERENCE}" "$REFERENCE"

log INFO " 🧺 Fetching samples from NCBI list or checking samples in folders ..."
# === MODULE 03: Samples fetching ===
modules/pipeline1/03_input_fetch.sh

log INFO " 🔬 Performing raw fastqc analysis.."
# === MODULE 04: Checking Reference ===
modules/pipeline1/04_fastqc_parallel.sh --threads "$THREADS"

log INFO " ✂️ Trimming samples fastqc..."
#=== Validate adapter ====
ADAPTER_FILE="adapters/${ADAPTER}.adapters"
if [[ ! -f "$ADAPTER_FILE" ]]; then
  log ERROR "❌ Adapter file not found: $ADAPTER_FILE"
  log INFO  "📂 Available adapter profiles:"
  ls adapters/*.adapters | sed 's|adapters/||;s/.adapters$//'
  exit 1
fi
# === MODULE 05: Checking Reference ===
modules/pipeline1/05_Cutadapt_trimming_phix_parallel.sh --threads "$THREADS" "$ADAPTER"

log INFO " 🔍 Performing trimmed fastqc analysis.. ..."
# === MODULE 06: Checking Reference ===
modules/pipeline1/06_fastqc_trimmed_parallel.sh --threads "$THREADS"

log INFO " 🪤 Detecting exogenous spike 🐝 if present ..."
# === MODULE 07: Checking Reference ===
modules/pipeline1/07_spike_detect.sh --ref-prefix "$REFERENCE"

log INFO " 🔗 Bwa-mem paired alignment ..."
# === MODULE 08: Checking Reference ===
modules/pipeline1/08_alignment_bwa_spike.sh --format "$REPORT_FORMAT" -t "$THREADS" Reference "$REFERENCE"

log INFO " 📇  Adding platform read groups ..."
# === MODULE 09: Checking Reference ===
modules/pipeline1/09_readgroups_add.sh "$PLATFORM"

log INFO " ⚗️ Filtering and cleaning 🛁 🧼 BAM files ..."
# === MODULE 10: Checking Reference ===
modules/pipeline1/10_bam_cleaning.sh --ref Reference --prefix "$REFERENCE"

log INFO "renaming samples based in metadata for further analyis and grouping ..."
# === MODULE 11: Checking Reference ===

modules/pipeline1/11_Renaming_bam.sh -m "$MAPPING"

log INFO "🎉 Pipeline completed."
end_time=$(date +%s)
elapsed=$((end_time - start_time))
log INFO "⏱️ Pipeline completed in ${elapsed}s."
