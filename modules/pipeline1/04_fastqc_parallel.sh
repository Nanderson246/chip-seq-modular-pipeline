#!/usr/bin/env bash
# Module: 04_fastqc.sh (Parallel Optimized)
# Author: Nancy Anderson
# Description: Quality control with FastQC and MultiQC using parallel processing
################################################################################
# SOFTWARE REQUIREMENTS:
#   • bash                >= 4.x
#   • fastqc              >= 0.11.9
#   • multiqc             >= 1.12
#   • GNU parallel        >= 20221122
#   • coreutils           (for commands: date, basename, nproc, etc.)
#   • unzip               (used for validating FastQC zip outputs)
#   • findutils           (for 'find')
#   • awk, ps, df         (standard POSIX tools for performance logging)
################################################################################
################################################################################
# Usage:
#   Normal run: bash modules/pipeline1/04_fastqc.sh
#   Dry run:    DRY_RUN=true bash modules/pipeline1/04_fastqc.sh
#   Standalone: ./04_fastqc.sh --skip-multiqc
#   Single thread use 
#   ./bash modules/pipeline1/04_fastqc.sh --threads 1 --skip-multiqc
################################################################################
set -uo pipefail

: "${DRY_RUN:=false}"

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

#=== Log Function ===
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

# === SCRIPT ===
readonly VERSION="3.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}" 
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
SKIP_MULTIQC=false

# === Path Resolution ===
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"


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

readonly SYSTEM
readonly CORES
readonly MAX_JOBS

log INFO "🧠 Scheduler detected: $SYSTEM, using $MAX_JOBS parallel jobs"



# === Argument Parsing ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-multiqc)
            SKIP_MULTIQC=true
            shift
            ;;
        --threads)
            MAX_JOBS_OVERRIDE="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage:"
            echo "  ./04_fastqc.sh                # Normal run"
            echo "  DRY_RUN=true ./04_fastqc.sh   # Dry-run simulation"
            echo "  ./04_fastqc.sh --skip-multiqc # Run without MultiQC"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            exit 1
            ;;
    esac
done

# === Constants ===
INPUT_DIR="${PROJECT_ROOT}/samples"
OUTPUT_DIR="${PROJECT_ROOT}/results/QC_fastqc"
MULTIQC_REPORT="raw_multiqc_report"
LOG_DIR="${PROJECT_ROOT}/logs/pipeline1"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"

# === Resource Monitoring ===
record_metrics() {
    local message="$1"
    local timestamp=$(date +%s)
    local cpu_usage=$(ps -p $$ -o %cpu | tail -n 1 | awk '{print $1}')
    local mem_usage=$(ps -p $$ -o %mem | tail -n 1 | awk '{print $1}')
    local disk_usage=$(df -h "${OUTPUT_DIR}" | tail -n 1)
    local file_count=$(find "${OUTPUT_DIR}" -name '*_fastqc.*' 2>/dev/null | wc -l)
    
    echo "[${timestamp}] ${message}" >> "$PERF_LOG"
    echo "  CPU: ${cpu_usage}% | Memory: ${mem_usage}% | QC Reports: ${file_count}" >> "$PERF_LOG"
    echo "  Disk: ${disk_usage}" >> "$PERF_LOG"
}

start_timer() {
    SECONDS=0
}

get_elapsed_time() {
    local duration=$SECONDS
    echo "$((duration / 60))m $((duration % 60))s"
}

# === Initialize Environment ===
mkdir -p "$LOG_DIR" "$OUTPUT_DIR"
start_timer

{
log "INFO" ""
log "INFO" "=============================================="
log "INFO" "🧬 Module: ${SCRIPT_BASE_NAME} (Parallel v${VERSION})"
log "INFO" "📌 Purpose: Quality Control of raw data"
log "INFO" "💡 Dry-run mode: ${DRY_RUN:-false}"
log "INFO" "💡 Parallel mode: ${PARALLEL_AVAILABLE} (Max jobs: ${MAX_JOBS})"
log "INFO" "📁 Input dir: ${INPUT_DIR}/"
log "INFO" "📁 Output dir: ${OUTPUT_DIR}/"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "=============================================="

# === Verify Required Tools ===
for tool in fastqc; do
    if ! command -v "$tool" &>/dev/null; then
        log "ERROR" "❌ Required tool '$tool' not found in PATH"
        exit 1
    fi
done

log "INFO" "🔧 Tool versions:"
log "INFO" "• fastqc:  $(fastqc --version 2>&1 | head -n 1)"
[ "$SKIP_MULTIQC" = false ] && log "INFO" "• multiqc: $(multiqc --version 2>&1 | head -n 1)"

# === Dry-run Simulation ===
if [ "${DRY_RUN:-false}" = true ]; then
    echo ""
    log "INFO" "🔍 [DRY-RUN] Simulation:"
    record_metrics "DRY-RUN START"
    
    # Check input files
    if compgen -G "${INPUT_DIR}/*.fastq.gz" > /dev/null; then
        file_count=$(find "${INPUT_DIR}" -name "*.fastq.gz" | wc -l)
        log "INFO" "📋 Found ${file_count} FASTQ files that would be processed"
        
        echo ""
        log "INFO" "⚡ Would execute:"
        log "INFO" "fastqc -o ${OUTPUT_DIR} ${INPUT_DIR}/*.fastq.gz"
        [ "$SKIP_MULTIQC" = false ] && log "INFO" "multiqc -o ${OUTPUT_DIR} -n ${MULTIQC_REPORT} ${OUTPUT_DIR}/"
        
        record_metrics "DRY-RUN: Would process ${file_count} files"
    else
        log "ERROR" "❌ No FASTQ files found in ${INPUT_DIR}/"
        record_metrics "DRY-RUN: No inputs found"
    fi
    
    echo ""
    log "INFO" "✅ [DRY-RUN] Simulation completed (no changes made)"
    exit 0
fi

# === Main Execution ===
log "INFO" "🔍 Checking input files..."
record_metrics "PROCESS START"

# Find FASTQ files
mapfile -t fastq_files < <(find "$INPUT_DIR" -maxdepth 1 -name "*.fastq.gz" -type f 2>/dev/null | sort -V)

if [ ${#fastq_files[@]} -eq 0 ]; then
    log "ERROR" "❌ Error: No FASTQ files found in ${INPUT_DIR}/"
    record_metrics "PROCESS FAILED - NO INPUTS"
    exit 1
fi

log "INFO" "✅ Found ${#fastq_files[@]} FASTQ file(s)"
record_metrics "INPUT DETECTED: ${#fastq_files[@]} files"

# === Parallel Processing Function ===
process_fastq() {
    local fq="$1"
    local fq_base=$(basename "$fq")
    local fq_prefix="${fq_base%.fastq.gz}"
    local fq_zip="${OUTPUT_DIR}/${fq_prefix}_fastqc.zip"
    local fq_html="${OUTPUT_DIR}/${fq_prefix}_fastqc.html"
    local done_marker="${OUTPUT_DIR}/${fq_prefix}.done"

    # Check for existing valid output
    if [[ -f "$fq_zip" && -f "$done_marker" ]]; then
        if [[ -s "$fq_zip" ]] && unzip -tq "$fq_zip" &>/dev/null; then
            log "INFO" "✓ Already processed: $fq_base"
            return 0
        else
            log "WARN" "⚠️ $fq_zip is corrupted — deleting for reprocessing"
            rm -f "$fq_zip" "$fq_html" "$done_marker"
        fi
    fi

    log "INFO" "⚡ Processing: $fq_base"
    if fastqc -o "$OUTPUT_DIR" "$fq"; then
        touch "$done_marker"
        log "INFO" "✅ Completed: $fq_base"
        return 0
    else
        log "WARN" "⚠️ FastQC failed for $fq_base"
        rm -f "$fq_zip" "$fq_html" "$done_marker"
        return 1
    fi
}

# Export function and variables for parallel
export -f process_fastq log
export OUTPUT_DIR RED GREEN YELLOW BLUE NC

# === Parallel Execution ===
processed_count=0
if $PARALLEL_AVAILABLE; then
    log "INFO" "🚀 Starting parallel processing (${MAX_JOBS} jobs)"
    
    # Create file list for parallel
    tmp_file_list=$(mktemp)
    printf "%s\n" "${fastq_files[@]}" > "$tmp_file_list"
    
    # Run parallel processing
    parallel -j $MAX_JOBS --bar --joblog "$OUTPUT_DIR/fastqc_joblog" \
        "process_fastq {}" < "$tmp_file_list"
    
    # Get processed count from job log
    processed_count=$(awk 'NR>1 && $7==0 {count++} END {print count+0}' "$OUTPUT_DIR/fastqc_joblog")
    rm "$tmp_file_list"
else
    log "INFO" "⏳ Starting sequential processing"
    for fq in "${fastq_files[@]}"; do
        if process_fastq "$fq"; then
            ((processed_count++))
        fi
    done
fi

# === Run MultiQC (unless skipped) ===
if [ "$SKIP_MULTIQC" = true ]; then
    log "INFO" "⏭️ Skipping MultiQC as requested via --skip-multiqc"
    record_metrics "⏭️ MULTIQC SKIPPED - Flag set"
else
    echo ""
    log "INFO" "🔗 Aggregating results with MultiQC..."
    record_metrics "MULTIQC START"

    if [ "$processed_count" -gt 0 ] || [ ! -f "${OUTPUT_DIR}/${MULTIQC_REPORT}.html" ]; then
        [ -f "${OUTPUT_DIR}/${MULTIQC_REPORT}.html" ] && rm -f "${OUTPUT_DIR}/${MULTIQC_REPORT}"*

        if multiqc -o "$OUTPUT_DIR" -n "$MULTIQC_REPORT" "$OUTPUT_DIR/"; then
            log "INFO" "✅ MultiQC report generated"
            record_metrics "MULTIQC COMPLETE"
        else
            log "ERROR" "❌ Error: MultiQC failed"
            record_metrics "MULTIQC FAILED"
            exit 1
        fi
    else
        log "INFO" "✅ MultiQC report already exists and no new files processed"
        record_metrics "⏭️ MULTIQC SKIPPED - No new data"
    fi
fi

# === Completion ===
log "INFO" ""
log "INFO" "=============================================="
log "INFO" "✅ MODULE COMPLETED"
log "INFO" "📊 Processed ${processed_count} new FASTQ file(s)"
log "INFO" "📂 Output FASTQC files: $(find "$OUTPUT_DIR" -name '*_fastqc.zip' | wc -l)"
log "INFO" "⏱️  Elapsed time: $(get_elapsed_time)"
log "INFO" "📁 QC reports directory: ${OUTPUT_DIR}/"
[ "$SKIP_MULTIQC" = false ] && log "INFO" "📄 MultiQC report: ${OUTPUT_DIR}/${MULTIQC_REPORT}.html"
log "INFO" "📊 Performance metrics saved to: ${PERF_LOG}"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "=============================================="

record_metrics "PROCESS COMPLETE"

} | tee -a "$MODULE_LOG"
