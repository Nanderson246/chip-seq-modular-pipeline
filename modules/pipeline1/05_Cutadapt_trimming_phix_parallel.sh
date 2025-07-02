#!/usr/bin/env bash
# Module: 05_Cutadapt_trimming_phix.sh (Parallel Optimized)
# Author: Nancy Anderson
# Description: Trim paired-end FASTQ files using parallel processing
#cutadapt --version 4.4

################################################################################
# SOFTWARE REQUIREMENTS:
#   • bash                >= 4.x
#   • cutadapt            >= 4.4
#   • GNU parallel        >= 20221122
#   • coreutils           (for: date, basename, nproc, etc.)
#   • findutils           (for: find)
#   • awk, ps, df         (standard POSIX tools for performance logging)
#   • gzip                (used for validating .gz files)
#   • flock               (from util-linux, for safe parallel log writes)
################################################################################

################################################################################
# USAGE:
# Inside the packages
#   Normal run: bash modules/pipeline1/05_Cutadapt_trimming_phix_paralel.sh [options] [adapter_profile]
#Standalone use:
#bash 05_Cutadapt_trimming_phix.sh --input-dir /data/myfastqs --output-dir /mnt/trimmed tn5_truseq
# OPTIONS:
#   -d, --dry-run    Simulate execution without making changes
#   -h, --help       Show this help message
#   --input-dir <path>    Override input directory (default: samples)
#   --output-dir <path>   Override output directory (default: results/Trimmed)
#
# PARAMETERS:
#   adapter_profile  Name of adapter profile (default: tn5_truseq)
#
# EXAMPLES:
#   bash modules/pipeline1/05_Cutadapt_trimming_phix.sh  tn5_truseq
#   bash modules/pipeline1/05_Cutadapt_trimming_phix.sh --help
################################################################################
set -uo pipefail

# === SCRIPT ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}" 
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Path Resolution ===
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# === Color Definitions ===
if [[ -t 1 ]]; then
     RED='\033[0;31m'
     GREEN='\033[0;32m'
     YELLOW='\033[0;33m'
     BLUE='\033[0;34m'
     NC='\033[0m'
else
     RED='' GREEN='' YELLOW='' BLUE='' NC=''
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




# === CLI Parser ===
DRY_RUN=false
CUSTOM_INPUT=""
CUSTOM_OUTPUT=""
ADAPTER_PROFILE="tn5_truseq"
HELP_MSG=$'\nUsage: 05_Cutadapt_trimming.sh [options] [adapter_profile]\nOptions:\n  -d, --dry-run  Simulate execution\n  -h, --help     Show help'

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "$HELP_MSG"
            exit 0
            ;;
        --input-dir)                         
            CUSTOM_INPUT="$2"
            shift 2
            ;;
        --output-dir)                        
            CUSTOM_OUTPUT="$2"
            shift 2
            ;;
        --threads)
            MAX_JOBS_OVERRIDE="$2"
            shift 2
            ;;    
        *)
            ADAPTER_PROFILE="$1"
            shift
            ;;
            
    esac
done

# === Initialization ===
INPUT_DIR="${PROJECT_ROOT}/samples"
OUTPUT_DIR="${PROJECT_ROOT}/results/Trimmed"
[[ -n "$CUSTOM_INPUT" ]] && INPUT_DIR="$CUSTOM_INPUT"        
[[ -n "$CUSTOM_OUTPUT" ]] && OUTPUT_DIR="$CUSTOM_OUTPUT" 

ADAPTERS_DIR="${PROJECT_ROOT}/adapters"
LOG_DIR="${PROJECT_ROOT}/logs/pipeline1"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"
META_TSV="${PROJECT_ROOT}/metadata/mapping.tsv"
PHIX_ADAPTER="${ADAPTERS_DIR}/phix.adapters"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR"

# === Parallel Processing Setup ===
PARALLEL_AVAILABLE=$(command -v parallel &>/dev/null && echo true || echo false)
if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    SYSTEM="SLURM"
    CORES="${SLURM_CPUS_ON_NODE:-$(nproc)}"
elif [[ -n "${PBS_JOBID:-}" ]]; then
    SYSTEM="PBS"
    CORES="${PBS_NUM_PPN:-$(nproc)}"
elif [[ -n "${LSB_JOBID:-}" ]]; then
    SYSTEM="LSF"
    CORES="${LSB_DJOB_NUMPROC:-$(nproc)}"
else
    SYSTEM="Local"
    CORES="$(nproc)"
fi


DEFAULT_JOBS=$(( CORES * 50 / 100 ))
DEFAULT_JOBS=$(( DEFAULT_JOBS < 1 ? 1 : DEFAULT_JOBS ))
MAX_JOBS="${MAX_JOBS_OVERRIDE:-$DEFAULT_JOBS}"

readonly SYSTEM CORES MAX_JOBS

log "INFO" "🧠 Scheduler detected:$PARALLEL_AVAILABLE, $SYSTEM, using $MAX_JOBS parallel jobs"


# === Resource Monitoring Functions ===
record_metrics() {
    local message="$1"
    local timestamp=$(date +%s)
    local cpu_usage=$(ps -p $$ -o %cpu | tail -n 1 | awk '{print $1}')
    local mem_usage=$(ps -p $$ -o %mem | tail -n 1 | awk '{print $1}')
    local disk_usage=$(df -h "${OUTPUT_DIR}" | tail -n 1)
    local trimmed_count=$(find "${OUTPUT_DIR}" -name '*_trimmed.fastq.gz' 2>/dev/null | wc -l)
    
    # Use locking for parallel-safe writes
    (
        flock -x 200
        echo "[${timestamp}] ${message}" >> "$PERF_LOG"
        echo "  CPU: ${cpu_usage}% | Memory: ${mem_usage}% | Trimmed Files: ${trimmed_count}" >> "$PERF_LOG"
        echo "  Disk: ${disk_usage}" >> "$PERF_LOG"
    ) 200>"$PERF_LOG.lock"
}

start_timer() {
    SECONDS=0
}

get_elapsed_time() {
    local duration=$SECONDS
    echo "$((duration / 60))m $((duration % 60))s"
}

# === Processing Function for Parallel ===
process_sample() {
    local read1="$1"
    local base=$(basename "$read1" _1.fastq.gz)
    local read2="${INPUT_DIR}/${base}_2.fastq.gz"
    local trim1="${OUTPUT_DIR}/${base}_1_trimmed.fastq.gz"
    local trim2="${OUTPUT_DIR}/${base}_2_trimmed.fastq.gz"
    local done_marker="${OUTPUT_DIR}/${base}.trimmed.done"
    local clean_r1="${OUTPUT_DIR}/${base}_1_clean.fastq.gz"
    local clean_r2="${OUTPUT_DIR}/${base}_2_clean.fastq.gz"

    log "INFO" "🔍 Processing sample: $base"
    
    # Skip if already processed
    if [[ -f "$done_marker" ]]; then
        if [[ -s "$trim1" && -s "$trim2" ]] && gzip -t "$trim1" &>/dev/null && gzip -t "$trim2" &>/dev/null; then
            log "INFO" "✅ Already trimmed: $base - skipping"
            echo "skipped"
            return 0
        else
            log "WARN" "⚠️ Output files missing or corrupted for $base - reprocessing"
            rm -f "$trim1" "$trim2" "$done_marker"
        fi
    fi

    # Check for R2 file
    if [[ ! -f "$read2" ]]; then
        log "ERROR" "❌ Missing R2 file for ${base}, skipping"
        echo "error"
        return 1
    fi

    # Determine if we need to remove phiX
    local REMOVE_PHIX=true
    if [[ -f "$META_TSV" ]]; then
        local spike_type=$(awk -F'\t' -v sample="$base" '
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if ($i == "Sample") sample_col = i;
                    if ($i == "Spike_Type") spike_col = i;
                }
                next
            }
            $sample_col == sample && spike_col > 0 {print $spike_col; exit}' "$META_TSV")
        [[ "$spike_type" == "phiX174" ]] && REMOVE_PHIX=false
    fi

    # Remove phiX if needed
    if $REMOVE_PHIX; then
        log "INFO" "🧽 Removing phiX from $base"
        if ! cutadapt -j 1 -a file:"$PHIX_ADAPTER" -A file:"$PHIX_ADAPTER" \
            --discard-trimmed -o "$clean_r1" -p "$clean_r2" "$read1" "$read2"; then
            log "ERROR" "❌ phiX removal failed for $base"
            echo "error"
            return 1
        fi
    else
        cp "$read1" "$clean_r1"
        cp "$read2" "$clean_r2"
    fi

    # Main trimming with adapter removal
    log "INFO" "✂️ Trimming adapters from $base"
    local CMD=(cutadapt -j 1 -a "$ADAPT_R1" -A "$ADAPT_R2" -m 20 -o "$trim1" -p "$trim2")
    [ -n "$GTRIM_R1" ] && CMD+=(-g "$GTRIM_R1")
    [ -n "$GTRIM_R2" ] && CMD+=(-G "$GTRIM_R2")
    CMD+=("$clean_r1" "$clean_r2")

    if ! "${CMD[@]}"; then
        log "ERROR" "❌ Trimming failed for $base"
        rm -f "$trim1" "$trim2" "$done_marker"
        echo "error"
        return 1
    fi

    # Cleanup and mark completion
    touch "$done_marker"
    rm -f "$clean_r1" "$clean_r2"
    log "INFO" "✅ Successfully trimmed: $base"
    echo "processed"
    return 0
}

# Export function for parallel
export -f process_sample log record_metrics
export ADAPT_R1 ADAPT_R2 GTRIM_R1 GTRIM_R2 OUTPUT_DIR INPUT_DIR PERF_LOG META_TSV PHIX_ADAPTER
export RED GREEN YELLOW BLUE NC  # Export color variables

# === Main Execution ===
main() {
    start_timer
    touch "$PERF_LOG.lock"  # Create lock file
    
    log "INFO" ""
    log "INFO" "=============================================="
    log "INFO" "🧬 Module: ${SCRIPT_BASE_NAME} (Parallel v${VERSION})"
    log "INFO" "📌 Adapter profile: $ADAPTER_PROFILE"
    log "INFO" "📂 Input: $INPUT_DIR"                       
    log "INFO" "📂 Output: $OUTPUT_DIR"                     
    log "INFO" "💡 Dry-run mode: $DRY_RUN"
    log "INFO" "💡 Parallel mode: $PARALLEL_AVAILABLE (Max jobs: $MAX_JOBS)"
    log "INFO" "🕒 Start time: $(date '+%F %T')"
    log "INFO" "=============================================="

    # === Verify Required Tool ===
    command -v cutadapt >/dev/null || { log "ERROR" "❌ cutadapt not found in PATH"; exit 1; }
    
    # === Log Tool Versions ===
    log "INFO" "🛠️  cutadapt path: $(command -v cutadapt)"
    log "INFO" "🔧 Tool version: ✂️ Cutadapt: $(cutadapt --version 2>&1 | head -n 1)"

    # === Dry run ===
    if [ "$DRY_RUN" = true ]; then
        echo ""
        log "INFO" "🛑 DRY-RUN MODE: No files will be modified"
        record_metrics "DRY-RUN START"
        # ... (dry run simulation remains the same)
        exit 0
    fi

    record_metrics "PROCESS START"

    ADAPTER_FILE="${ADAPTERS_DIR}/${ADAPTER_PROFILE}.adapters"
    if [[ ! -f "$ADAPTER_FILE" ]]; then
        log "ERROR" "❌ Error: Adapter profile not found: ${ADAPTER_FILE}"
        log "INFO" "🎨 Available profiles:"
        ls "${ADAPTERS_DIR}/"*.adapters | sed "s|${ADAPTERS_DIR}/||;s/.adapters$//"
        record_metrics "PROCESS FAILED - INVALID ADAPTER PROFILE"
        exit 1
    fi

    source "$ADAPTER_FILE"
    export ADAPT_R1 ADAPT_R2 GTRIM_R1 GTRIM_R2

    log "INFO" "🔍 Scanning input directory: ${INPUT_DIR}"
    mapfile -t read1_files < <(find "${INPUT_DIR}" -maxdepth 1 -name "*_1.fastq.gz" -print | sort)
    total_pairs=${#read1_files[@]}
    log "INFO" "📊 Total read pairs detected: $total_pairs"
    record_metrics "INPUT DETECTED: $total_pairs pairs"

    # === Parallel Execution ===
    processed=0
    skipped=0
    errors=0
    
    if $PARALLEL_AVAILABLE; then
        log "INFO" "🚀 Starting parallel processing with $MAX_JOBS jobs"
        
        # Create temp file for results
        result_file=$(mktemp)
        
        # Run processing in parallel
        printf "%s\n" "${read1_files[@]}" | parallel -j $MAX_JOBS --bar --joblog "$LOG_DIR/cutadapt_joblog" \
            "process_sample {} >> $result_file"
        
        # Collect results
        while read -r result; do
            case $result in
                processed) ((processed++)) ;;
                skipped) ((skipped++)) ;;
                error) ((errors++)) ;;
            esac
        done < "$result_file"
        rm "$result_file"
    else
        log "INFO" "⏳ Starting sequential processing"
        for read1 in "${read1_files[@]}"; do
            result=$(process_sample "$read1")
            case $result in
                processed) ((processed++)) ;;
                skipped) ((skipped++)) ;;
                error) ((errors++)) ;;
            esac
        done
    fi

    echo ""
    log "INFO" "=============================================="
    log "INFO" "✅ TRIMMING COMPLETED"
    log "INFO" "📊 Newly processed: $processed pairs"
    log "INFO" "📊 Skipped existing: $skipped pairs"
    log "INFO" "📊 Errors encountered: $errors pairs"
    log "INFO" "⏱️  Elapsed time: $(get_elapsed_time)"
    log "INFO" "📂 Output location: ${OUTPUT_DIR}/"
    log "INFO" "📊 Performance metrics saved to: ${PERF_LOG}"
    log "INFO" "🕒 End time: $(date '+%F %T')"
    log "INFO" "=============================================="

    record_metrics "PROCESS COMPLETE: Processed $processed pairs"
}

main | tee -a "$MODULE_LOG"
