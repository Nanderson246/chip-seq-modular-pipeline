#!/usr/bin/env bash
# Module: 03_input_fetch.sh
# Author: Nancy Anderson
# Description: Fetch samples from SRA or detect pre-existing FASTQs with comprehensive monitoring

################################################################################
# SOFTWARE REQUIREMENTS
#
# Required Tools:
#   • bash            – Unix shell interpreter (v4+ recommended)
#   • prefetch        – From NCBI SRA Toolkit (for downloading SRA files)
#   • fasterq-dump    – From NCBI SRA Toolkit (for FASTQ conversion)
#   • gzip            – For compressing FASTQ files
#   • tee             – For logging output (used if PIPELINE_LOG_MODE=hybrid)
#   • ps, top, free   – For CPU/memory monitoring (from procps package)
#   • uptime, df      – For system load and disk usage reporting
#
# Optional Tools:
#   • fastq-dump      – Legacy FASTQ extractor (not used but logged)
#
# Installation (Debian/Ubuntu example):
#   sudo apt install sra-toolkit coreutils procps
#
# Notes:
#   - fasterq-dump is preferred over fastq-dump for performance
#   - Script assumes access to the internet for SRA downloads
#   - This script supports both test (dry-run) and actual data fetch modes
################################################################################

################################################################################
# Usage:
#   Normal run: ./03_input_fetch.sh
#   Dry run:    DRY_RUN=true ./03_input_fetch.sh
# Standalone:
# ./03_input_fetch.sh --sra-list path/to/custom_list.txt
################################################################################

set -euo pipefail

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

print_header() {
    echo "========================================"
    echo "$1"
    echo "========================================"
} 


[ "${DRY_RUN:-false}" = true ] && set -x  # Command echo only in dry-run

# === Optional Argument Parsing ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        --sra-list)
            SRA_LIST="$2"
            shift 2
            ;;
        --help|-h)
            echo "Usage: $0 [--sra-list path/to/SRR_Acc_List.txt]"
            exit 0
            ;;
        *)
            echo "❌ Unknown option: $1"
            exit 1
            ;;
    esac
done

# === Constants ===
SAMPLES_DIR="samples"
SRA_LIST="${SRA_LIST:-metadata/SRR_Acc_List.txt}"
LOG_DIR="logs"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"

# === Initialize Performance Log ===
mkdir -p "$LOG_DIR" "$SAMPLES_DIR"
if [[ ! -f "$PERF_LOG" ]]; then
    log "INFO" "Timestamp | Stage | CPU (%) | MEM (%) | LOAD | DISK | SAMPLES" > "$PERF_LOG"
    echo "---------------------------------------------------------------" >> "$PERF_LOG"
fi

# === Verify Required Commands ===
command -v prefetch >/dev/null || { echo "❌ Error: prefetch not found in PATH"; exit 1; }
command -v fastq-dump >/dev/null || { echo "❌ Error: fastq-dump not found in PATH"; exit 1; }
#command -v fasterq-dump >/dev/null || { echo "❌ Error: fasterq-dump not found in PATH"; exit 1; }

# === Log Tool Versions ===
log "INFO" "🔧 Tool versions:"
log "INFO" "• prefetch:     $(prefetch --version 2>&1 | head -n 1)"
#log "INFO" "• fastq-dump:   $(fastq-dump --version 2>&1 | head -n 1)"
log "INFO" "• fasterq-dump: $(fasterq-dump --version 2>&1 | head -n 1)"


# === Resource Monitoring Functions ===
record_metrics() {
    local message="$1"
    local timestamp=$(date +"%Y-%m-%d %H:%M:%S")
    
    local cpu_usage=$(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1}')
    local mem_usage=$(free -m | awk '/Mem:/ {printf "%.1f", $3/$2*100}')
    local load_avg=$(uptime | awk -F'load average: ' '{print $2}' | cut -d, -f1)
    local disk_usage=$(df -h "${SAMPLES_DIR}" | tail -1 | awk '{print $5}')
    local sample_count=$(ls -1 "${SAMPLES_DIR}"/*.fastq.gz 2>/dev/null | wc -l)
    
    printf "%-19s | %-20s | %6s | %6s | %4s | %4s | %5d\n" \
        "$timestamp" "$message" "$cpu_usage" "$mem_usage" "$load_avg" "$disk_usage" "$sample_count" >> "$PERF_LOG"
}

monitor_process() {
    local pid=$1
    local tool_name=$2
    local sample_id=$3
    
    while kill -0 "$pid" 2>/dev/null; do
        local process_metrics=$(ps -p "$pid" -o %cpu=,%mem= --no-headers 2>/dev/null || echo "0.0 0.0")
        printf "[%s] %s metrics for %s: CPU %s MEM %s\n" \
            "$(date +"%H:%M:%S")" "$tool_name" "$sample_id" $process_metrics >> "$PERF_LOG"
        sleep 10
    done
}

run_with_metrics() {
    local tool_name="$1"
    local sample_id="$2"
    shift 2
    local args=("$@")

    record_metrics "START $tool_name:$sample_id"

    "$tool_name" "${args[@]}" &
    local pid=$!

    monitor_process "$pid" "$tool_name" "$sample_id" &
    local monitor_pid=$!

    wait "$pid"
    local exit_code=$?

    kill "$monitor_pid" 2>/dev/null || true

    record_metrics "END $tool_name:$sample_id"

    return $exit_code
}

start_timer() {
    SECONDS=0
}

get_elapsed_time() {
    local duration=$SECONDS
    echo "$((duration / 60))m $((duration % 60))s"
}

# === Initialize Environment ===
mkdir -p "$SAMPLES_DIR"
start_timer

{
print_header "🧬 Module: ${SCRIPT_NAME} v${VERSION}"
log "INFO" "📌 Purpose: Fetch samples from NCBI or detect existing FASTQs"
log "INFO" "💡 Dry-run mode: ${DRY_RUN:-false}"
print_header "🕒 Start time: $(date '+%F %T')"


# === Dry-run Simulation ===
if [ "${DRY_RUN:-false}" = true ]; then
    log "INFO" "🔍 [DRY-RUN] Simulation:"
    record_metrics "DRY-RUN START"
    
    if compgen -G "${SAMPLES_DIR}/*_1.fastq.gz" > /dev/null; then
        log "INFO" "✅ Would use existing FASTQ files in ${SAMPLES_DIR}/"
        ls -1 "${SAMPLES_DIR}/"*_[12].fastq.gz | head -2
        echo "..."
    elif [ -f "$SRA_LIST" ]; then
        log "INFO" "📋 Would download SRA accessions from: ${SRA_LIST}"
        log "INFO" "Sample IDs that would be processed:"
        head -n 3 "$SRA_LIST"
        [ $(wc -l < "$SRA_LIST") -gt 3 ] && echo "..."
        log "INFO" "→ Output dir: ${SAMPLES_DIR}/"
    else
        log "ERROR" "❌ Would fail: No FASTQs found and no ${SRA_LIST}"
    fi
    
    record_metrics "DRY-RUN COMPLETE"
    echo ""
    log "INFO" "✅ [DRY-RUN] Simulation completed (no changes made)"
    exit 0
fi

# === Main Execution ===
echo ""
log "INFO" "🔍 Checking input sources..."
record_metrics "PROCESS START"

if [ -f "$SRA_LIST" ]; then
    log "INFO" "📋 SRA accession list found. Processing..."
    record_metrics "SRA DOWNLOAD START"

    total_samples=$(wc -l < "$SRA_LIST")
    processed=0

    while IFS= read -r sample_id || [ -n "$sample_id" ]; do
        processed=$((processed + 1))
        out_prefix="${SAMPLES_DIR}/${sample_id}"

        if [[ -f "${out_prefix}_1.fastq.gz" && -f "${out_prefix}_2.fastq.gz" ]]; then
            log "INFO" "✅ Already processed: ${sample_id} (${processed}/${total_samples})"
            continue
        fi

        if [[ -f "${SAMPLES_DIR}/${sample_id}.sra" && ! -f "${out_prefix}_1.fastq.gz" ]]; then
            log "INFO" "🔁 Resuming FASTQ conversion for ${sample_id}..."
            if ! run_with_metrics "fasterq-dump" "$sample_id" --split-files --threads 4 --outdir "$SAMPLES_DIR" "${SAMPLES_DIR}/${sample_id}"; then
                log "ERROR" "❌ FASTQ conversion failed for ${sample_id}"
                record_metrics "CONVERSION FAILED: ${sample_id}"
                continue
            fi
            
            if ! gzip -f "${SAMPLES_DIR}/${sample_id}_1.fastq" "${SAMPLES_DIR}/${sample_id}_2.fastq"; then
                log "ERROR" "❌ Gzip compression failed for ${sample_id}"
                record_metrics "GZIP FAILED: ${sample_id}"
                continue
            fi
            rm -f "${SAMPLES_DIR}/${sample_id}.sra"
            record_metrics "SAMPLE COMPLETE: ${sample_id}"
            continue
        fi
        
        echo ""
        log "INFO" "⬇️ Downloading ${sample_id} (${processed}/${total_samples})..."
        
        # 🧹 Clean up potential stale lock file
        lock_file="${SAMPLES_DIR}/${sample_id}/${sample_id}.sra.lock"
        if [[ -f "$lock_file" ]]; then
            log "WARN" "🔓 Removing stale lock file: $lock_file"
             rm -f "$lock_file"
        fi

        if ! run_with_metrics "prefetch" "$sample_id" -O "$SAMPLES_DIR" "$sample_id"; then
            log "WARN" "⚠️ Download timed out, retrying..."
            if ! run_with_metrics "prefetch" "$sample_id" -O "$SAMPLES_DIR" "$sample_id"; then
                log "ERROR" "❌ Failed to download ${sample_id} after retry"
                record_metrics "DOWNLOAD FAILED: ${sample_id}"
                continue
            fi
        fi

        log "INFO" "🔧 Converting to FASTQ..."
        if ! run_with_metrics "fasterq-dump" "$sample_id" --split-files --threads 4 --outdir "$SAMPLES_DIR" "${SAMPLES_DIR}/${sample_id}"; then
            log "ERROR" "❌ FASTQ conversion failed for ${sample_id}"
            record_metrics "CONVERSION FAILED: ${sample_id}"
            continue
        fi
        
        if ! gzip -f "${SAMPLES_DIR}/${sample_id}_1.fastq" "${SAMPLES_DIR}/${sample_id}_2.fastq"; then
            log "ERROR" "❌ Gzip compression failed for ${sample_id}"
            record_metrics "GZIP FAILED: ${sample_id}"
            continue
        fi

        rm -f "${SAMPLES_DIR}/${sample_id}.sra"
        record_metrics "SAMPLE COMPLETE: ${sample_id}"

    done < "$SRA_LIST"

    record_metrics "SRA DOWNLOAD COMPLETE"

else
    echo ""
    log "ERROR" "❌ ERROR: No input sources found"
    log "INFO" "🧰 Possible solutions:"
    log "INFO" "1. Place existing FASTQ files in ${SAMPLES_DIR}/ (named as *_1.fastq.gz and *_2.fastq.gz)"
    log "INFO" "2. Provide an ${SRA_LIST} file with SRA accessions"
    record_metrics "PROCESS FAILED - NO INPUTS"
    exit 1
fi

#=== Remove NCBI binary folder to save space ===
for dir in SRR*/; do
    sample=${dir%/}
    if [[ -f "${sample}_1.fastq.gz" && -f "${sample}_2.fastq.gz" ]]; then
        log "INFO" "🧹 Removing original SRA folder: $sample"
        rm -r "$sample"
    else
        log "INFO" "⏭️ Skipping $sample — FASTQ pair missing"
    fi
done

# === Completion ===
echo ""
print_header "✅ MODULE COMPLETED"
log "INFO" "🕒 [MODULE: ${SCRIPT_BASE_NAME}] Completed at $(date)"
log "INFO" "⏱️  Elapsed time: $(get_elapsed_time)"
log "INFO" "📂 Output location: ${SAMPLES_DIR}/"
print_header "📊 Performance metrics saved to: ${PERF_LOG}"

record_metrics "PROCESS COMPLETE"
} | if [ "${DRY_RUN:-false}" != true ]; then 
    tee -a "$MODULE_LOG" 
else
    cat
fi

