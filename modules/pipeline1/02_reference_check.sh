#!/usr/bin/env bash

#Module: 02_reference_check.sh
# Author: Nancy Anderson
# Description: Validates reference genome files with auto-recovery and dry-run support
################################################################################
# SOFTWARE REQUIREMENTS:
#   • bash                 >= 4.x
#   • samtools             >= 1.10       (used for quickcheck on .fa files)
#   • coreutils            (basename, date, mkdir, etc.)
#   • awk, grep, find      (standard POSIX tools, assumed available)
#   • Optional: reference_setup.sh      (local helper script for auto-recovery)
################################################################################
################################################################################
# Usage: 
#   Normal run: bash modules/pipeline1/02_reference_check.sh /path/to/reference [genome_prefix]
#   Dry run:    DRY_RUN=true ./02_reference_check.sh /path/to/reference [genome_prefix]
# standalone usage 
# ./02_reference_check.sh /path/to/reference [genome_prefix]
################################################################################

set -uo pipefail

# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}"  
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Logging Configuration ===
LOG_DIR="logs/pipeline1"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"


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

# ------------------------------------------------------------------------------
# 0. Resolve script location and package root
# ------------------------------------------------------------------------------

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_ROOT="${SCRIPT_DIR%/modules/pipeline1}"
UTIL_DIR="${PACKAGE_ROOT}/modules/utils"
DEFAULT_REF_ROOT="${PACKAGE_ROOT}/Reference"

mkdir -p "$LOG_DIR"
log "INFO" "📎 Ensured logs directory exists: $LOG_DIR"

# === Tool availability check ===
command -v samtools >/dev/null || {
    log "ERROR" "❌ Required tool 'samtools' not found in PATH"
    exit 1
}

[ "${DRY_RUN:-false}" = true ] && set -x  # Only enable command echo in dry-run

SECONDS=0  # Timer starts

# ------------------------------------------------------------------------------
# 1. Argument parsing and default reference path
# ------------------------------------------------------------------------------

REF_PREFIX="${2:-hg38}"

# -----------------------------------------------------------
#   Map genome → expected GENCODE GTF filename (unzipped)
# -----------------------------------------------------------
case "$REF_PREFIX" in
    hg38) GTF_REQUIRED="gencode.v45.annotation.gtf" ;;
    mm10) GTF_REQUIRED="gencode.vM31.annotation.gtf" ;;
    *)    log "ERROR" "Unsupported genome prefix: $REF_PREFIX"; exit 1 ;;
esac

if [[ $# -ge 1 ]]; then
    REF_FOLDER="$1"
else
    REF_FOLDER="${DEFAULT_REF_ROOT}/${REF_PREFIX}"
    log "INFO" "📁 No reference path provided — defaulting to ${REF_FOLDER}"
fi

log "INFO" "📍 Checking reference folder: $REF_FOLDER"
log "INFO" "🔓 Using reference prefix: $REF_PREFIX"

# === Reference constants ===
declare -a REQUIRED_FILES=(
    "${REF_PREFIX}.fa"
    "${REF_PREFIX}.fa.fai"
    "${REF_PREFIX}.dict"
    "${REF_PREFIX}.fa.amb"
    "${REF_PREFIX}.fa.ann"
    "${REF_PREFIX}.fa.bwt"
    "${REF_PREFIX}.fa.pac"
    "${REF_PREFIX}.fa.sa"
    "${REF_PREFIX}-blacklist.bed"
    "$GTF_REQUIRED"          # genome-specific GENCODE GTF
    "tss_annotations.bed"
)

{
log "INFO" ""
log "INFO" "==============================================="
log "INFO" "🧬 Module: ${SCRIPT_NAME} v${VERSION}"
log "INFO" "🕒 Started: $(date '+%F %T')"
log "INFO" "📁 Reference folder: $REF_FOLDER"
log "INFO" "🔤 Reference prefix: $REF_PREFIX"
log "INFO" "💡 Dry-run mode: ${DRY_RUN:-false}"
log "INFO" "=============================================="


# === Dry-run Simulation ===
if [ "${DRY_RUN:-false}" = true ]; then
    log "INFO" "🔍 [DRY-RUN] Would check these files:"
    for file in "${REQUIRED_FILES[@]}"; do
        log "INFO" " - ${REF_FOLDER}/${file}"
    done
    
    log "INFO" "🔧 [DRY-RUN] Would attempt recovery if missing:"
    log "INFO" " - reference_setup.sh ${REF_FOLDER} ${REF_PREFIX}"
    
    log "INFO" "✅ [DRY-RUN] Simulation completed (no changes made)"
    exit 0
fi

# === Core Validation Function ===
validate_reference() {
    local missing=0
    
    log "INFO" "🔍 Checking reference files:"
    for file in "${REQUIRED_FILES[@]}"; do
        if [[ -f "${REF_FOLDER}/${file}" ]]; then
            log "INFO" "✅  ${file}"
        else
            log "ERROR" "❌  MISSING: ${file}"
            ((missing++))
        fi
    done
    
 
    if ((missing > 0)); then
        log "ERROR" "⚠️  Missing $missing required file(s)"
        return 1
    fi
    
    log "INFO" "✅  All reference files validated"
    return 0
}

# === Main Execution ===
if ! validate_reference; then
    log "INFO" "🛠️  Attempting auto-recovery..."
    if [[ ! -f "${REF_FOLDER}/${REF_PREFIX}.fa" || ! -f "${REF_FOLDER}/${REF_PREFIX}-blacklist.bed" ]]; then
        log "ERROR" "❌  Critical files missing - running reference setup"
          mkdir -p "$REF_FOLDER"  # ✅ Only create if truly needed
        bash "${BASH_SOURCE%/*}/reference_setup.sh" "$REF_FOLDER" "$REF_PREFIX" || {
            log "ERROR" "❌  FAILED: Auto-recovery unsuccessful"
            exit 1
        }
        validate_reference || exit 1
    else
        log "WARN" "⚠️  Only auxiliary files missing - please generate indices manually"
        exit 1
    fi
fi


# === Final Checks ===
log "INFO" "🧪 Verifying genome FASTA integrity manually..."

if [[ ! -s "${REF_FOLDER}/${REF_PREFIX}.fa" ]]; then
    log "ERROR" "❌ Genome FASTA is empty or missing"
    exit 1
fi

if ! grep -q "^>chr" "${REF_FOLDER}/${REF_PREFIX}.fa"; then
    log "WARN" "⚠️ Genome headers may be non-standard. Expected >chr*"
fi

if [[ ! -s "${REF_FOLDER}/${REF_PREFIX}.fa.fai" ]]; then
    log "ERROR" "❌ FASTA index (.fai) file missing or empty"
    exit 1
fi

log "INFO" "✅ Genome FASTA and index look okay"



# ✅ NEW: Print example structure if passed   <-- keep as comment or delete
log "INFO" ""
log "INFO" "📁 Expected folder structure:"
log "INFO" "└── ${REF_FOLDER}/"
for file in "${REQUIRED_FILES[@]}"; do log "INFO" "    ├── $file"; done


log "INFO" ""
log "INFO" "=============================================="
log "INFO" "✅ MODULE COMPLETED: All checks passed"
log "INFO" "🕒 Duration: ${SECONDS}s"
log "INFO" "🕒 [MODULE: ${SCRIPT_BASE_NAME}] Completed at $(date)"
log "INFO" "=============================================="
} | if [ "${DRY_RUN:-false}" != true ]; then 
    tee -a "$MODULE_LOG" 
else
    cat  # Show dry-run output directly
fi


