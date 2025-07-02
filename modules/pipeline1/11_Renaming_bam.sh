#!/usr/bin/env bash
# Module: 11_Renaming_bam.sh
# Author: Nancy Anderson
# Description: Rename cleaned BAM files according to metadata mapping with resource monitoring

################################################################################
# SOFTWARE REQUIREMENTS
#
# This script requires the following tools to be installed and accessible in PATH:
#
# Required Tools:
#   • bash        - Unix shell interpreter (v4+ recommended)
#   • awk         - For metadata field parsing
#   • ps          - For CPU/memory usage tracking
#   • top         - System-wide CPU stats
#   • free        - Memory usage stats (from `procps`)
#   • df          - Disk usage reporting
#   • iostat      - Disk I/O metrics (optional; from `sysstat` package)
#
# Optional:
#   • tee         - For dual logging (stdout and file)
#
# Installation (Debian/Ubuntu-based systems):
#   sudo apt install coreutils procps sysstat gawk
#
# Notes:
#   - 'iostat' is optional; script will warn if unavailable but still function.
#   - Qualimap and MultiQC are not required here, but are expected to have
#     already been run in previous pipeline steps if QC directories are used.
#
################################################################################

################################################################################
# USAGE:
#   Normal run: bash modules/pipeline1/11_Renaming_bam.sh -m metadata/mapping.tsv
#   Dry run:    bash 11_Renaming_bam.sh -m metadata/mapping.tsv --dry-run
#USAGE STANDALONE:
#bash 09_Renaming_bam.sh \
#  -m my_mapping.tsv \
#  --input-dir ./bam_inputs \
#  --output-dir ./renamed_bams \
#  --log-dir ./custom_logs \
#  --metrics-dir ./qc_metrics
# OPTIONS:
#   -m, --mapping FILE    Metadata mapping file (required)
#   -d, --dry-run         Simulate execution without changes
#   -h, --help            Show this help message
#
# EXAMPLES:
#   bash 11_Renaming_bam.sh -m metadata/mapping.tsv
#   bash 11_Renaming_bam.sh --mapping metadata.csv --dry-run
################################################################################
set -uo pipefail

# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
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
# === Resource Tracking Variables ===
declare -A PROCESS_METRICS=()
declare -A SYSTEM_METRICS=()
declare -a TEMP_FILES=()

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
start_monitoring() {
    local stage="$1"
    local pid="$2"
    
    PROCESS_METRICS["${stage}_start"]=$(date +%s.%N)
    PROCESS_METRICS["${stage}_pid"]=$pid
    
    # Capture initial system metrics
    capture_system_metrics "${stage}_start"
}

stop_monitoring() {
    local stage="$1"
    local end_time=$(date +%s.%N)
    local start_time=${PROCESS_METRICS["${stage}_start"]}
    local pid=${PROCESS_METRICS["${stage}_pid"]}
    
    # Calculate process metrics
    PROCESS_METRICS["${stage}_time"]=$(awk -v s="$start_time" -v e="$end_time" 'BEGIN {print e - s}')
    PROCESS_METRICS["${stage}_cpu"]=$(ps -p $pid -o %cpu --no-headers | awk '{print $1}')
    PROCESS_METRICS["${stage}_mem"]=$(ps -p $pid -o rss --no-headers | awk '{printf "%.2f", $1/1024}') # MB
    
    # Capture final system metrics
    capture_system_metrics "${stage}_end"
}

capture_system_metrics() {
    local prefix="$1"
    local timestamp=$(date +%s.%N)
    
    # CPU usage (system-wide)
    local cpu=$(top -bn1 | grep "Cpu(s)" | sed "s/.*, *\([0-9.]*\)%* id.*/\1/" | awk '{print 100 - $1}')
    
    # Memory usage
    local mem_total=$(free -m | awk '/Mem:/ {print $2}')
    local mem_used=$(free -m | awk '/Mem:/ {print $3}')
    local mem_pct=$(awk -v used="$mem_used" -v total="$mem_total" 'BEGIN {printf "%.1f", (used/total)*100}')
    
    if command -v iostat >/dev/null; then
        local disk_read=$(iostat -d | awk '/^[a-z]+/ {print $3}' | head -1)
        local disk_write=$(iostat -d | awk '/^[a-z]+/ {print $4}' | head -1)
    else
        log "WARN" "iostat not available — skipping disk IO metrics"
        local disk_read=0
        local disk_write=0
    fi
    
    SYSTEM_METRICS["${prefix}_timestamp"]=$timestamp
    SYSTEM_METRICS["${prefix}_cpu"]=$cpu
    SYSTEM_METRICS["${prefix}_mem_used"]=$mem_used
    SYSTEM_METRICS["${prefix}_mem_pct"]=$mem_pct
    SYSTEM_METRICS["${prefix}_disk_read"]=$disk_read
    SYSTEM_METRICS["${prefix}_disk_write"]=$disk_write
}



# === CLI Parser ===
DRY_RUN=false
MAPPING_FILE=""
HELP_MSG=$'\nUsage: Renaming_bam.sh -m <mapping_file> [options]
Options:
  -m, --mapping FILE     Metadata mapping file
  -d, --dry-run          Simulate execution
  --input-dir DIR        Path to cleaned BAMs (default: results/Filtered/Cleaned)
  --output-dir DIR       Path to renamed BAM output (default: analysis/Renamed_Cleaned)
  --log-dir DIR          Log directory (default: logs)
  --metrics-dir DIR      QC metrics directory (default: results/Filtered/Metrics)
  -h, --help             Show help'

while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mapping)
            MAPPING_FILE="$2"
            shift 2
            ;;       
        --input-dir)
            INPUT_DIR="$2"
            shift 2
            ;;
        --output-dir)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        --log-dir)
            LOG_DIR="$2"
            shift 2
            ;;   
        --metrics-dir)
            METRICS_DIR="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "$HELP_MSG"
            exit 0
            ;;
        -*)
            echo "❌ Error: Unknown option $1"
            echo "$HELP_MSG"
            exit 1
            ;;
        *)
            echo "❌ Error: Positional arguments not allowed"
            echo "$HELP_MSG"
            exit 1
            ;;
    esac
done

# === Initialization ===
INPUT_DIR="${INPUT_DIR:-results/Filtered/Cleaned}"
OUTPUT_DIR="${OUTPUT_DIR:-analysis/Renamed_Cleaned}"
LOG_DIR="${LOG_DIR:-logs/pipeline1}"
METRICS_DIR="${METRICS_DIR:-results/Filtered/Metrics}"
QC_DIR="${OUTPUT_DIR}/QC"
MODULE_LOG="${LOG_DIR}/09_Renaming_bam_$(date +%F_%H-%M-%S).log"
PERF_LOG="${LOG_DIR}/09_Renaming_bam_performance_metrics.log"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR" "$QC_DIR"

# === Resource Monitoring Functions ===
record_metrics() {
    local message="$1"
    local timestamp=$(date +%s)
    local cpu_usage=$(ps -p $$ -o %cpu | tail -n 1 | awk '{print $1}')
    local mem_usage=$(ps -p $$ -o %mem | tail -n 1 | awk '{print $1}')
    local disk_usage=$(df -h "${OUTPUT_DIR}" | tail -n 1)
    local renamed_count=$(ls -1 "${OUTPUT_DIR}"/*.bam 2>/dev/null | wc -l)
    
    echo "[${timestamp}] ${message}" >> "$PERF_LOG"
    echo "  CPU: ${cpu_usage}% | Memory: ${mem_usage}% | Renamed BAMs: ${renamed_count}" >> "$PERF_LOG"
    echo "  Disk: ${disk_usage}" >> "$PERF_LOG"
}

start_timer() {
    SECONDS=0
}

get_elapsed_time() {
    local duration=$SECONDS
    echo "$((duration / 60))m $((duration % 60))s"
}

# === Dry-run Setup ===
if [ "$DRY_RUN" = true ]; then
    log "INFO" -e "\n🛑 DRY-RUN MODE: No files will be modified"
    set -x  # Enable command echo
fi

start_timer

{
echo "=============================================="
echo "🧬 Module: 11_Renaming_bam.sh"
echo "📌 Purpose: Rename BAMs according to metadata"
echo "💡 Dry-run mode: $DRY_RUN"
echo "📁 Input mapping: $MAPPING_FILE"
echo "🕒 Start time: $(date '+%F %T')"
echo "=============================================="

record_metrics "PROCESS START"

# === Validate mapping file ===
if [[ -z "$MAPPING_FILE" ]]; then
    log "ERROR" "❌ Error: Mapping file not provided"
    record_metrics "PROCESS FAILED - NO MAPPING FILE"
    exit 1
elif [[ ! -f "$MAPPING_FILE" ]]; then
    log "ERROR" "❌ Error: Mapping file not found: $MAPPING_FILE"
    record_metrics "PROCESS FAILED - MAPPING FILE NOT FOUND"
    exit 1
fi

# === Detect delimiter ===
case "$MAPPING_FILE" in
    *.csv) DELIM=',' ;;
    *.tsv|*.txt) DELIM=$'\t' ;;
    *) 
        log "ERROR" "❌ Error: Unsupported file extension. Use .tsv, .csv, or .txt"
        record_metrics "PROCESS FAILED - INVALID FILE EXTENSION"
        exit 1
        ;;
esac

# === Check for duplicate Sample_IDs ===
if [[ $(cut -d"$DELIM" -f1 "$MAPPING_FILE" | sort | uniq -d | wc -l) -gt 0 ]]; then
    log "ERROR" "❌ Duplicate Sample_IDs found in mapping file"
    exit 1
fi



# === Dry-run Simulation ===
if [ "$DRY_RUN" = true ]; then
    log "INFO" "🔍 SIMULATION RESULTS:"
    record_metrics "DRY-RUN START"
    
    total_samples=$(($(wc -l < "$MAPPING_FILE") - 1))
    log "INFO" "Found ${total_samples} samples in mapping file"
    
    log "INFO" "\nFirst 5 sample IDs that would be processed:"
    head -n 5 "$MAPPING_FILE" | cut -f1 | sed '1d; s/^/  /'
    [ "$total_samples" -gt 5 ] && echo "  ... (and $((total_samples - 5)) more)"
    
    log "INFO" "⚡ Example renaming that would occur:"
    log "INFO" "  From: ${INPUT_DIR}/SAMPLE_clean.bam"
    log "INFO" "  To:   ${OUTPUT_DIR}/SAMPLE_CONDITION_REP_IP_TARGET.bam"
    
    log "INFO" "📊 Would generate outputs in:"
    log "INFO" "  Renamed BAMs: ${OUTPUT_DIR}/"
    log "INFO""  QC Reports: ${QC_DIR}/"
    
    record_metrics "DRY-RUN COMPLETE: Would process ${total_samples} samples"
    log "INFO" "\n✅ [DRY-RUN] Simulation completed (no changes made)"
    exit 0
fi

# === Main Processing ===
processed=0
skipped=0
errors=0
total_samples=$(($(wc -l < "$MAPPING_FILE") - 1))

log "INFO" "🔧 Processing ${total_samples} samples from mapping file..."
record_metrics "PROCESSING START: ${total_samples} samples"

HEADER=$(head -n 1 "$MAPPING_FILE")
IFS="$DELIM" read -ra COL_NAMES <<< "$HEADER"

while IFS="$DELIM" read -ra FIELDS; do
    declare -A meta
#    for i in "${!COL_NAMES[@]}"; do
#        meta["${COL_NAMES[$i]}"]="${FIELDS[$i]}"
#    done

     for i in "${!COL_NAMES[@]}"; do 
         meta["${COL_NAMES[$i]}"]="${FIELDS[$i]:-}"
    done

    Sample_ID="${meta[Sample_ID]}"
    
    # Validate required fields
    if [[ -z "$Sample_ID" || -z "${meta[Sample_Type]}" || -z "${meta[Condition]}" || -z "${meta[Target]}" ]]; then
        log "WARN" "⚠️  Missing required field(s) in line: ${FIELDS[*]}"
        record_metrics "PROCESS WARNING: Missing fields in ${Sample_ID:-UNKNOWN}"
        ((errors++))
        continue
    fi
    
    # Validate Sample_Type values
#if [[ ! "${meta[Sample_Type]}" =~ ^(IP|ChIP|Input|IgG|Mock) ]]; then
if [[ ! "${meta[Sample_Type]}" =~ ^(Input|IgG|Mock|ChIP|IP(_rep[0-9]+)?)$ ]]; then
    log "WARN" "⚠️ Invalid Sample_Type '${meta[Sample_Type]}' for ${Sample_ID}"
fi
        
    # Input files
    input_bam="${INPUT_DIR}/${Sample_ID}_clean.bam"
    input_bai="${input_bam}.bai"
    
    if [[ ! -f "$input_bam" ]]; then
        log "WARN" "⚠️  BAM not found: $input_bam"
        record_metrics "PROCESS WARNING: Missing BAM for ${Sample_ID}"
        ((errors++))
        continue
    fi
    

    # Build new name
parts=("$Sample_ID")
    [[ -n "${meta[Instrument]}" ]]  && parts+=("$(echo "${meta[Instrument]}" | tr -s ' ' '_' | tr -d '/')")
    [[ -n "${meta[Condition]}" ]]   && parts+=("$(echo "${meta[Condition]}" | tr -s ' ' '_' | tr -d '/')")
    [[ -n "${meta[Replicate]}" ]]   && parts+=("$(echo "${meta[Replicate]}" | tr -s ' ' '_' | tr -d '/')")
    [[ -n "${meta[Sample_Type]}" ]] && parts+=("$(echo "${meta[Sample_Type]}" | tr -s ' ' '_' | tr -d '/')")
    [[ -n "${meta[Target]}" ]]      && parts+=("$(echo "${meta[Target]}" | tr -s ' ' '_' | tr -d '/')")


    new_base=$(IFS="_"; echo "${parts[*]}")
    output_bam="${OUTPUT_DIR}/${new_base}.bam"
    output_bai="${output_bam}.bai"

    # Skip if already exists
    if [[ -f "$output_bam" ]]; then
        log "INFO" "✅Already renamed: $output_bam"
        record_metrics "SAMPLE SKIPPED: ${Sample_ID} (already exists)"
        ((skipped++))
        continue
    fi

    # Perform copy with metrics
    log "INFO" "🗃️ Renaming: ${Sample_ID} → ${new_base}.bam"
    record_metrics "RENAMING START: ${Sample_ID} → ${new_base}"
    
    if ! cp "$input_bam" "$output_bam"; then
        log "ERROR" "❌ Error: Failed to copy BAM for ${Sample_ID}"
        record_metrics "PROCESS FAILED: BAM copy for ${Sample_ID}"
        ((errors++))
        continue
    fi
    
    # Copy BAI index if it exists
    if [[ -f "$input_bai" ]]; then
        cp "$input_bai" "$output_bai"
    else
        log "WARN" "⚠️ Missing index for ${input_bam}"
    fi
    
    # Copy QC files
    old_flagstat="${METRICS_DIR}/${Sample_ID}_clean_flagstat.txt"
    new_flagstat="${QC_DIR}/${new_base}_flagstat.txt"
    [[ -f "$old_flagstat" ]] && cp "$old_flagstat" "$new_flagstat"

    old_qualimap_dir="${METRICS_DIR}/${Sample_ID}_qualimap_clean"
    new_qualimap_dir="${QC_DIR}/${new_base}_qualimap"
    [[ -d "$old_qualimap_dir" ]] && cp -r "$old_qualimap_dir" "$new_qualimap_dir"

    ((processed++))
    record_metrics "SAMPLE COMPLETE: ${Sample_ID} → ${new_base}"
done < <(tail -n +2 "$MAPPING_FILE" | grep -v '^[[:space:]]*$')


# === Completion ===
log "INFO" "\n=============================================="
log "INFO" "🎉🟢 RENAMING COMPLETED"
log "INFO" "🔬 Total samples: $total_samples"
log "INFO" "📊 Processed: $processed files"
log "INFO" "⏭️ Skipped: $skipped already renamed files"
log "INFO" "⚠️  Errors: $errors"
log "INFO" "⏱️  Elapsed time: $(get_elapsed_time)"
log "INFO" "📂 Output directory: $OUTPUT_DIR/"
log "INFO" "📝 Performance metrics saved to: ${PERF_LOG}"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "=============================================="

record_metrics "PROCESS COMPLETE: ${processed} samples renamed"
} | tee -a "$MODULE_LOG"
