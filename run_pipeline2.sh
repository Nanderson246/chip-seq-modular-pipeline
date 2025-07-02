#!/usr/bin/env bash
# Script: run_pipeline2.sh
# Author: Nancy Anderson • 2025-06-22
# Description: Entry point for Pipeline 2 — QC → merge → peak calling (MACS3 or HOMER) → IDR → annotation
################################################################################
# USAGE:
#   bash run_pipeline2.sh --caller macs3 --genome mm10
#
# OPTIONS:
#   --caller CALLER    Required: peak caller to use [macs3 or homer]
#
#   -t THREADS         Number of threads to use (default: 4)
#   -m MAPPING         Path to metadata file (default: metadata/mapping.tsv)
#   -f FORCE_PBC       Force recalculation of PBC complexity check (default: false)
#   -c CORR_METHOD     Correlation method for plots and DeepTools (default: pearson)
#
# EXAMPLES:
#   bash run_pipeline2.sh --caller macs3 -t 4  --genome mm10
#   bash run_pipeline2.sh --caller homer -m my_mapping.tsv -f true -c spearman
#
# DESCRIPTION:
#   This script runs the second stage of the ChIP-seq pipeline:
#   1. Performs replicate QC
#   2. Merges replicates and creates pseudoreplicates
#   3. Branches into either MACS3 or HOMER for peak calling
#   4. Runs IDR consistency analysis
#   5. Annotates reproducible peaks
#
#   Requires a validated metadata mapping file.
#   Output directories are auto-generated inside the analysis/ and logs/ folders.
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

# === Default Values ===
THREADS=4
MAPPING="metadata/mapping.tsv"
FORCE_PBC=false
CORRELATION_METHOD="pearson"
CALLER=""
GENOME_NAME="hg38"

usage() {
 echo "Usage: $0 --caller [macs3|homer] [-t THREADS] [-m MAPPING] [-f FORCE_PBC] [-c CORR_METHOD] [--genome GENOME]"
  exit 1
}

# === Parse CLI Arguments ===
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t) THREADS="$2"; shift 2 ;;
    -m) MAPPING="$2"; shift 2 ;;
    -f) FORCE_PBC="$2"; shift 2 ;;
    -c) CORRELATION_METHOD="$2"; shift 2 ;;
    --genome) GENOME_NAME="$2"; shift 2 ;;
    --caller) CALLER="$2"; shift 2 ;;
    -h|--help) usage ;;
    *) echo "❌ Unknown argument: $1"; usage ;;
  esac
done

# === Validate Required Caller ===
if [[ "$CALLER" != "macs3" && "$CALLER" != "homer" ]]; then
  echo "❌ --caller must be 'macs3' or 'homer'"
  usage
fi

# === Start ===
log INFO "🧬 Running Pipeline 2"
log INFO "📁 Metadata: $MAPPING"
log INFO "⚙️ Threads: $THREADS"
log INFO "🧪 Force PBC: $FORCE_PBC"
log INFO "📊 Correlation Method: $CORRELATION_METHOD"
log INFO "🔔 Peak Caller: $CALLER"
log INFO "🧬 Genome: $GENOME_NAME"

# === Metadata Validation ===
python3 modules/utils/validate_mapping.py --mapping "$MAPPING" --schema templates/mapping_schema.yaml
if [[ $? -ne 0 ]]; then
  log ERROR  "❌ Metadata validation failed"
  exit 1
fi

#=== Genome validation ===

if [[ -z "$GENOME_NAME" ]]; then
  log ERROR "Genome name (--genome) is missing or empty."
  exit 1
fi

log INFO "🧪 Running replicate QC module..."
# === MODULE 01: Replicate QC ===
modules/pipeline2/01_replicate_qc.sh \
  --meta "$MAPPING" \
  --threads "$THREADS" \
  --correlation-method "$CORRELATION_METHOD" \
  $( [[ "$FORCE_PBC" == "true" ]] && echo "--force-pbc" )

log INFO "🔗 Merging pooled and pseudoreplicates..."
# === MODULE 02: Merge + Pseudoreps ===
modules/pipeline2/02_merge_pseudoreplicates.sh \
  --meta metadata/mapping_filtered.tsv \
  --correlation-method "$CORRELATION_METHOD"

# ===BRANCH: MACS3 vs HOMER ===
if [[ "$CALLER" == "macs3" ]]; then
  log INFO "🚨 MACS3 PATH SELECTED"

log INFO "🧫 Peak calling on filtered replicates..."
  # 03_1: Peak calling on filtered replicates
  modules/pipeline2/03_1_MACS3_peak_calling.sh --mapping metadata/mapping_filtered.tsv

log INFO "🧫 Peak calling on pooled and pseudoreplicates..."
  # 03_3: Peak calling on pooled/pseudoreps
  modules/pipeline2/03_3_MACS3_peak_calling_pooled_pseudoreps.sh

elif [[ "$CALLER" == "homer" ]]; then
  log INFO "🚨 HOMER PATH SELECTED"

 log INFO "🧫 Peak calling on filtered replicates..."
  # 03_2: Peak calling on filtered replicates
  modules/pipeline2/03_2_homer_peak_calling_fold_fdr_relaxed.sh --mapping metadata/mapping_filtered.tsv

log INFO "🧫 Peak calling on pooled and pseudoreplicates..."
  # 03_4: Peak calling on pooled/pseudoreps
  modules/pipeline2/03_4_homer_peak_calling_pooled_pseudoreps.sh
fi

log INFO "🔍 Running IDR Analysis..."
# === MODULE 04: IDR Analysis (for both branches) ===
modules/pipeline2/04_run_idr.sh --source "$CALLER"

log INFO "📌 Performing peak annotation and enrichment..."
# === MODULE 05: Peak Annotation ===
# or pull from metadata or args in future if needed
modules/pipeline2/05_peak_annotation_parallel.sh --idr "$CALLER" --genome "$GENOME_NAME"


# === Done ===
log INFO "✅ Pipeline 2 completed for caller: $CALLER"

