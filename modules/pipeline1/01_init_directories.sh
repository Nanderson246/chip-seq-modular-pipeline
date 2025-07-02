#!/usr/bin/env bash


# Module: 01_init_directories.sh
# Author: Nancy Anderson
# Description: Create all necessary directories for the pipeline
################################################################################
# SOFTWARE REQUIREMENTS
#
# This script requires the following tools installed and accessible in PATH:
#
# Required Tools:
#   • bash      – Standard Unix shell interpreter (version 4+ recommended)
#   • mkdir     – For creating directory structure
#   • date      – For generating timestamps
#   • tee       – For logging stdout/stderr to file (optional, for hybrid logging)
#
# Optional:
#   • color-capable terminal – For colored logging output
#
# Notes:

################################################################################
# USAGE:
#   Normal run: bash modules/pipeline1/01_init_directories.sh 
# Standalone usage:
# PIPELINE_RUNNING=true PIPELINE_LOG_MODE=hybrid bash modules/pipeline1/01_init_directories.sh

################################################################################

set -e  # Exit on error
set -x  # Echo commands for debugging/logging

# === Safe defaults for hybrid execution ===
: "${PIPELINE_RUNNING:=}"
: "${PIPELINE_LOG_MODE:=local}"

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

# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}"  
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Logging Configuration ===
LOG_DIR="logs"
readonly MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"

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


# === Ensure logs directory exists ===
mkdir -p "$LOG_DIR"
log "INFO" "📎 Ensured logs directory exists: $LOG_DIR"


# === Informational Banner ===
log "INFO" "==========================================="
log "INFO" "🧬 Module: ${SCRIPT_NAME} v${VERSION}"
log "INFO" "📌 Purpose: Create directory structure"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "=============================================="


# === Logging mode ===
if [[ -z "$PIPELINE_RUNNING" || "$PIPELINE_LOG_MODE" == "hybrid" ]]; then
    exec > >(tee -a "$MODULE_LOG") 2>&1
    
    log "INFO" "Logging to: $MODULE_LOG"
else
    log "INFO" "Centralized logging enabled, skipping per-module log."
fi

log "INFO" "[MODULE: 01_init_directories] Started at $(date +%F_%H-%M-%S)"

# === Define directory structure ===
DIRS=(
  samples logs tmp
  results/{QC_fastqc,Trimmed,BAM,Filtered/{Deduplicated,Cleaned,Metrics}}
  analysis/{Renamed_BAMs,Replicate_QC}
)

# ===Create all directories ===
for d in "${DIRS[@]}"; do
  mkdir -p "$d"
done


echo ""
log "INFO" "=============================================="
log "INFO" "✅ Directories created"
log "INFO" "[MODULE: ${SCRIPT_BASE_NAME}] Completed at $(date)"
log "INFO" "=============================================="
