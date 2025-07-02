#!/usr/bin/env bash
# Module: 06_fastqc_trimmed.sh (Parallel Optimized)
# Description: Run FastQC and MultiQC on trimmed FASTQ files

################################################################################
# SOFTWARE REQUIREMENTS
#
# Required Tools:
#   • bash           – Unix shell interpreter (v4+ recommended)
#   • fastqc         – Quality control for FASTQ files
#   • multiqc        – Aggregate reports from FastQC (can be skipped)
#   • parallel       – GNU Parallel for multicore processing (optional but recommended)
#   • unzip          – For validating FastQC ZIP outputs
#   • tee            – For hybrid logging (if enabled)
#   • ps, top, free  – System resource monitoring (from procps package)
#   • df, find       – Disk space and file discovery
#
# Installation Example (Debian/Ubuntu):
#   sudo apt install fastqc multiqc parallel unzip coreutils procps
#
# Notes:
#   - The script detects the job scheduler (SLURM, PBS, LSF) and adjusts thread usage.
#   - MultiQC execution can be skipped with --skip-multiqc.
#   - Logs performance metrics and supports dry-run for simulation mode.
################################################################################

################################################################################
# Usage:
#   Normal run: bash modules/pipeline1/06_fastqc_trimmed_parallel.sh
#   Dry run:    DRY_RUN=true bash modules/pipeline1/06_fastqc_trimmed_parallel.sh
# standalone usage: 
# ./06_fastqc_trimmed_parallel.sh --skip-multiqc
# DRY_RUN=true ./06_fastqc_trimmed_parallel --threads 1 --skip-multiqc [false/true]
# default  false
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

# === Path Resolution ===
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

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

readonly SYSTEM
readonly CORES
readonly MAX_JOBS

log INFO "🧠 Scheduler detected: $SYSTEM, using $MAX_JOBS parallel jobs"


# === Argument Parsing ===
SKIP_MULTIQC=false
MAX_JOBS_OVERRIDE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --skip-multiqc)
            SKIP_MULTIQC=true
            shift
            ;;
        --threads)
            if [[ -n "${2:-}" && "$2" =~ ^[0-9]+$ ]]; then
                MAX_JOBS_OVERRIDE="$2"
                shift 2
            else
                log "ERROR" "❌ Invalid or missing value for --threads"
                exit 1
            fi
            ;;
        --help|-h)
            echo "Usage:"
            echo "  $SCRIPT_NAME [--threads N] [--skip-multiqc]"
            echo "  DRY_RUN=true $SCRIPT_NAME  # Dry-run simulation"
            exit 0
            ;;
        *)
            log "ERROR" "❌ Unknown option: $1"
            exit 1
            ;;
    esac
done


# === Constants ===
INPUT_DIR="${PROJECT_ROOT}/results/Trimmed"
OUTPUT_DIR="${PROJECT_ROOT}/results/QC_trimmed_fastqc"
MULTIQC_REPORT="trimmed_multiqc_report"
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
echo ""
log "INFO" "=============================================="
log "INFO" "🧬 MODULE: ${SCRIPT_BASE_NAME} (Parallel v${VERSION})"
log "INFO" "📌 Purpose: Quality Control of trimmed data"
log "INFO" "💡 Dry-run mode: ${DRY_RUN:-false}"
log "INFO" "💡 Parallel mode: ${PARALLEL_AVAILABLE} (Max jobs: ${MAX_JOBS})"
log "INFO" "📁 Input dir: ${INPUT_DIR}/"
log "INFO" "📁 Output dir: ${OUTPUT_DIR}/"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "=============================================="

# === Verify Required Tools ===
for tool in fastqc; do
    if ! command -v "$tool" &>/dev/null; then
        log "ERROR" "❌ Required tool '$tool' not found"
        exit 1
    fi
done

log "INFO" "🔧 Tool versions:"
log "INFO" "• fastqc: $(fastqc --version | head -n 1)"
[ "$SKIP_MULTIQC" = false ] && log "INFO" "• multiqc: $(multiqc --version | head -n 1)"

# === Dry-run Simulation ===
if [ "${DRY_RUN:-false}" = true ]; then
    echo ""
    log "INFO" "🔍 [DRY-RUN] Simulation:"
    record_metrics "DRY-RUN START"
    
    # Dry-run input check (match real execution logic)
    if compgen -G "${INPUT_DIR}/*_trimmed.fastq.gz" > /dev/null; then
        file_count=$(find "${INPUT_DIR}" -name "*_trimmed.fastq.gz" | wc -l)
        log "INFO" "📋 Found ${file_count} trimmed FASTQ files that would be processed"

        
        echo ""
        log "INFO" "⚡ Would execute:"
        log "INFO" "fastqc -o ${OUTPUT_DIR} ${INPUT_DIR}/*_trimmed.fastq.gz"
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
echo ""
log "INFO" "🔍 Checking input files..."
record_metrics "PROCESS START"

mapfile -t fastq_files < <(find "$INPUT_DIR" -name "*_trimmed.fastq.gz" -type f | sort -V)

if [ ${#fastq_files[@]} -eq 0 ]; then
    log "ERROR" "❌ No trimmed FASTQ files found in ${INPUT_DIR}/"
    exit 1
fi

log "INFO" "✅ Found ${#fastq_files[@]} trimmed FASTQ file(s)"
record_metrics "INPUT DETECTED: ${#fastq_files[@]} files"

# === Parallel Processing Function ===
process_fastq() {
    local fq="$1"
    local fq_base=$(basename "$fq")
    local fq_prefix="${fq_base%.fastq.gz}"
    local fq_zip="${OUTPUT_DIR}/${fq_prefix}_fastqc.zip"
    local fq_html="${OUTPUT_DIR}/${fq_prefix}_fastqc.html"
    local done_marker="${OUTPUT_DIR}/${fq_prefix}.done"
    
    # Validate existing output
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

# Export function for parallel
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
            record_metrics "MULTIQC COMPLETE"
        else
            log "ERROR" "❌ MultiQC failed"
            record_metrics "MULTIQC FAILED"
            exit 1
        fi
    else
        log "INFO" "✅ MultiQC report already exists and no new files processed"
        record_metrics "⏭️ MULTIQC SKIPPED - No new data"
    fi
fi

# === Completion ===
echo ""
log "INFO" "=============================================="
log "INFO" "✅ MODULE COMPLETED"
log "INFO" "📦 FASTQC on trimmed files completed"
log "INFO" "📊 Processed: ${processed_count} files"
log "INFO" "🧾 Total FASTQC reports: $(find "$OUTPUT_DIR" -name '*_fastqc.zip' | wc -l)"
[ "$SKIP_MULTIQC" = false ] && log "INFO" "📄 MultiQC: ${OUTPUT_DIR}/${MULTIQC_REPORT}.html"
log "INFO" "⏱️ Duration: $(get_elapsed_time)"
log "INFO" "📁 QC reports directory: ${OUTPUT_DIR}/"
log "INFO" "📈 Performance metrics saved to: ${PERF_LOG}"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "=============================================="

record_metrics "PROCESS COMPLETE"
} | tee -a "$MODULE_LOG"
