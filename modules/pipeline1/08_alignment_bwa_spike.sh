#!/usr/bin/env bash
# Module: 08_alignment_bwa_spike.sh
# Author: Nancy Anderson
# Description: BWA alignment with detailed resource monitoring and multiple output formats
################################################################################
# SOFTWARE REQUIREMENTS:
#   - BWA (≥ 0.7.17)
#   - SAMtools (≥ 1.9)
#   - GNU Awk (≥ 4.0)
#   - coreutils (for readlink, dirname, date, etc.)
#   - procps-ng (for `free`, `top`)
#   - Optional: sysstat (for `iostat` if detailed disk I/O is desired)
#   - Bash (≥ 4.x) with associative arrays and advanced features
#
# Ensure all tools are available in $PATH.
################################################################################

################################################################################
# USAGE:
#   bash modules/pipeline1/ 08_alignment_bwa_spike [options] <reference_folder> [genome_prefix] 
#bash modules/pipeline1/08_alignment_bwa_spike.sh (when default hg38 and all foldere inside package)
#bash modules/pipeline1/08_alignment_bwa_spike.sh Reference mm10
#Standalone usage:
#bash 06_alignment_bwa.sh --in <folder_with-trimmed_fastq.gz> /
#     --out <otput_folder_bam> <Reference_folder>  <genome_prefix>
# NEW FEATURES:
#   - Precise CPU/Memory/Disk monitoring
#   - Nanosecond timing
#   - CSV/TSV performance reports
#   - System-wide resource logging
################################################################################

set -uo pipefail
shopt -s nullglob

# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}"

# === Path Resolution ===
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"

# === Configuration ===
OUTPUT_FORMATS=("human" "csv" "tsv")  # Supported report formats
REPORT_FORMAT="human"                 # Default format
LOG_DIR="logs/pipeline1"
LOG_DIR="${PROJECT_ROOT}/logs/pipeline1"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"
SYSTEM_METRICS_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_system_metrics_${TIMESTAMP}.csv"


# === CLI Parsing ===
parse_args() {
    THREADS=$(nproc --all 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)
    DRY_RUN=false
    INPUT_DIR="results/Trimmed"
    BAM_DIR="results/BAM"

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --format)
                if [[ " ${OUTPUT_FORMATS[@]} " =~ " $2 " ]]; then
                    REPORT_FORMAT="$2"
                    shift 2
                else
                    echo "Invalid format: $2. Supported: ${OUTPUT_FORMATS[*]}"
                    exit 1
                fi
                ;;
            -t|--threads)
                if [[ "$2" =~ ^[0-9]+$ ]] && [[ "$2" -gt 0 ]]; then
                    THREADS="$2"
                    shift 2
                else
                    echo "Invalid thread count: $2"
                    exit 1
                fi
                ;;
            --in)
                INPUT_DIR="$2"
                shift 2
                ;;
            --out)
                BAM_DIR="$2"
                shift 2
                ;;
            --meta-tsv)
                META_TSV="$2"
                shift 2
                ;;
            -d|--dry-run)
                DRY_RUN=true
                shift
                ;;
            -v|--version)
                echo "${SCRIPT_NAME} v${VERSION}"
                exit 0
                ;;
            -h|--help)
                echo "${SCRIPT_NAME} v${VERSION}"
                echo "Usage: $SCRIPT_NAME [options] <reference_folder> [genome_prefix]"
                echo "Options: ..."
                exit 0
                ;;
            *)
                if [[ -z "${REF_FOLDER:-}" ]]; then
                    REF_FOLDER="$1"
                else
                    REF_PREFIX="$1"
                fi
                shift
                ;;
        esac
    done

    [[ -z "${REF_FOLDER:-}" ]] && {
        echo "Reference folder not specified."
        exit 1
    }

    REF_PREFIX="${REF_PREFIX:-hg38}"
    readonly REF_GENOME="${REF_FOLDER}/${REF_PREFIX}/${REF_PREFIX}.fa"
    readonly FILTERED_DIR="results/Filtered"
    readonly DEDUP_DIR="${FILTERED_DIR}/Deduplicated"
}

# === Color Definitions ===
if [[ -t 1 ]]; then
    readonly RED='\033[0;31m'
    readonly GREEN='\033[0;32m'
    readonly YELLOW='\033[0;33m'
    readonly BLUE='\033[0;34m'
    readonly NC='\033[0m' # No Color
else
    readonly RED='' GREEN='' YELLOW='' BLUE='' NC=''
fi

# === Resource Tracking Variables ===
declare -A PROCESS_METRICS=()
declare -A SYSTEM_METRICS=()
declare -a TEMP_FILES=()

# === Functions ===
cleanup() {
    for file in "${TEMP_FILES[@]}"; do
        [[ -f "$file" ]] && rm -f "$file"
    done
    find "$BAM_DIR" -name "*.tmp" -mtime +1 -delete
    find "$BAM_DIR" -name "*.sort.tmp*" -mtime +1 -delete
}

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


check_requirements() {
    local tools=("bwa" "samtools" "awk" "top" "free")
    for tool in "${tools[@]}"; do
        if ! command -v "$tool" &>/dev/null; then
            log "ERROR" "Required tool not found: $tool"
            exit 1
        fi
    done
    log "INFO" "🔧 Tool versions:"
    log "INFO" "• bwa:          $(bwa 2>&1 |grep -i version | head -n 1)"
    log "INFO" "• samtools:     $(samtools --version | head -n 1)"
    log "INFO" "• awk:          $(awk --version 2>&1 | head -n 1 | cut -d ' ' -f1,2)"
    log "INFO" "• free:         $(free --version 2>&1 | head -n 1 || echo 'built-in')"
}


check_disk_space() {
    local required=$(( $(du -s "$INPUT_DIR" | awk '{print $1}') * 3 ))
    local available=$(df -P "$BAM_DIR" | awk 'NR==2 {print $4}')
    [[ $available -gt $required ]] || {
        log "ERROR" "Insufficient disk space (need $((required/1024))MB, have $((available/1024))MB)"
        exit 1
    }
}

validate_inputs() {
    [[ -d "$INPUT_DIR" ]] || {
        log "ERROR" "Input directory not found: $INPUT_DIR"
        exit 1
    }
    [[ $(ls "${INPUT_DIR}/"*_1_trimmed.fastq.gz 2>/dev/null | wc -l) -gt 0 ]] || {
        log "ERROR" "No input FASTQs found in $INPUT_DIR"
        exit 1
    }
}

validate_sample_pair() {
    local base="$1"
    local TRIM1="$2"
    local TRIM2="$3"
    [[ -f "$TRIM2" ]] || {
        log "ERROR" "Missing paired-end file for $base"
        return 1
    }
    return 0
}

start_monitoring() {
    local stage="$1"
    local pid="$2"
    
    PROCESS_METRICS["${stage}_start"]=$(date +%s.%N)
    PROCESS_METRICS["${stage}_pid"]=$pid
    
    # Capture initial system metrics
    if [[ "$REPORT_FORMAT" != "human" ]]; then
        capture_system_metrics "${stage}_start"
    fi
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
    if [[ "$REPORT_FORMAT" != "human" ]]; then
        capture_system_metrics "${stage}_end"
    fi
}

capture_system_metrics() {
    local prefix="$1"
    local timestamp=$(date +%s.%N)

    # CPU usage
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


generate_report() {
    case "$REPORT_FORMAT" in
        human)
            generate_human_report | tee -a "$PERF_LOG"
            ;;
        csv)
            generate_csv_report >> "$PERF_LOG"
            ;;
        tsv)
            generate_tsv_report >> "$PERF_LOG"
            ;;
    esac
}

generate_human_report() {
    echo -e "\n${BLUE}=== PERFORMANCE REPORT ===${NC}"
    printf "%-12s %10s %8s %10s %8s %8s\n" "Stage" "Time(s)" "CPU(%)" "Mem(MB)" "DiskR(MB/s)" "DiskW(MB/s)"
    echo "--------------------------------------------------"
    
    for stage in alignment sorting indexing; do
        printf "%-12s %10.2f %8.1f %10.1f %8.1f %8.1f\n" \
            "$stage" \
            "${PROCESS_METRICS["${stage}_time"]:-0}" \
            "${PROCESS_METRICS["${stage}_cpu"]:-0}" \
            "${PROCESS_METRICS["${stage}_mem"]:-0}" \
            "${SYSTEM_METRICS["${stage}_start_disk_read"]:-0}" \
            "${SYSTEM_METRICS["${stage}_start_disk_write"]:-0}"
    done
    
    local total_time=$(awk -v a="${PROCESS_METRICS["alignment_time"]:-0}" \
                          -v s="${PROCESS_METRICS["sorting_time"]:-0}" \
                          -v i="${PROCESS_METRICS["indexing_time"]:-0}" \
                          'BEGIN {print a + s + i}')
    
    echo "--------------------------------------------------"
    printf "%-12s %10.2f\n" "TOTAL" "$total_time"
    echo -e "${BLUE}==================================================${NC}"
}

generate_csv_report() {
    # CSV Header
    echo "stage,time_sec,cpu_percent,mem_mb,disk_read_mb,disk_write_mb,system_cpu,system_mem"
    
    # Data Rows
    for stage in alignment sorting indexing; do
        echo "${stage},${PROCESS_METRICS["${stage}_time"]:-0},${PROCESS_METRICS["${stage}_cpu"]:-0},"\
             "${PROCESS_METRICS["${stage}_mem"]:-0},${SYSTEM_METRICS["${stage}_start_disk_read"]:-0},"\
             "${SYSTEM_METRICS["${stage}_start_disk_write"]:-0},${SYSTEM_METRICS["${stage}_start_cpu"]:-0},"\
             "${SYSTEM_METRICS["${stage}_start_mem_pct"]:-0}"
    done
}

generate_tsv_report() {
    # TSV Header
    echo -e "stage\ttime_sec\tcpu_percent\tmem_mb\tdisk_read_mb\tdisk_write_mb\tsystem_cpu\tsystem_mem"
    
    # Data Rows
    for stage in alignment sorting indexing; do
        echo -e "${stage}\t${PROCESS_METRICS["${stage}_time"]:-0}\t${PROCESS_METRICS["${stage}_cpu"]:-0}\t"\
                "${PROCESS_METRICS["${stage}_mem"]:-0}\t${SYSTEM_METRICS["${stage}_start_disk_read"]:-0}\t"\
                "${SYSTEM_METRICS["${stage}_start_disk_write"]:-0}\t${SYSTEM_METRICS["${stage}_start_cpu"]:-0}\t"\
                "${SYSTEM_METRICS["${stage}_start_mem_pct"]:-0}"
    done
}

write_system_metrics_csv() {
    # Initialize system metrics log
    echo "timestamp,stage,event,cpu_percent,mem_used_mb,mem_percent,disk_read_mb,disk_write_mb" > "$SYSTEM_METRICS_LOG"
    
    # Append data for each stage
    for stage in alignment sorting indexing; do
        echo "${SYSTEM_METRICS["${stage}_start_timestamp"]},$stage,start,"\
             "${SYSTEM_METRICS["${stage}_start_cpu"]},"\
             "${SYSTEM_METRICS["${stage}_start_mem_used"]},"\
             "${SYSTEM_METRICS["${stage}_start_mem_pct"]},"\
             "${SYSTEM_METRICS["${stage}_start_disk_read"]},"\
             "${SYSTEM_METRICS["${stage}_start_disk_write"]}" >> "$SYSTEM_METRICS_LOG"
             
        echo "${SYSTEM_METRICS["${stage}_end_timestamp"]},$stage,end,"\
             "${SYSTEM_METRICS["${stage}_end_cpu"]},"\
             "${SYSTEM_METRICS["${stage}_end_mem_used"]},"\
             "${SYSTEM_METRICS["${stage}_end_mem_pct"]},"\
             "${SYSTEM_METRICS["${stage}_end_disk_read"]},"\
             "${SYSTEM_METRICS["${stage}_end_disk_write"]}" >> "$SYSTEM_METRICS_LOG"
    done
}


# === Main Alignment Function ===
validate_reference() {
    [[ -f "$REF_GENOME" ]] || {
        log "ERROR" "Reference genome not found at: $REF_GENOME"
        exit 1
    }
}

validate_bam() {
    local bam="$1"
    
    # Check BAM validity
    samtools quickcheck -v "$bam" || {
        log "ERROR" "Invalid BAM file: $bam"
        return 1
    }

    # Check for BAM index
    if [[ ! -f "${bam}.bai" ]]; then
        log "ERROR" "Missing BAM index file: ${bam}.bai"
        return 1
    fi
}


# === Alignment Function with spike-aware logic ===
# === Safe & Sequential Alignment Block ===
align_sample() {
    local base="$1"
    local TRIM1="$2"
    local TRIM2="$3"

    log "INFO" "Processing sample $base..."

    if [[ -f "${BAM_DIR}/${base}.sorted.bam.bai" ]]; then
        log "INFO" "Skipping sample $base (already indexed)"
        return 0
    fi

    local SPIKE_TYPE="none"
    if [[ -f "$META_TSV" ]]; then
        SPIKE_TYPE=$(awk -F'\t' -v sample="$base" '
            NR==1 {
                for (i=1; i<=NF; i++) {
                    if ($i == "Sample_ID") sample_col = i;
                    if ($i == "Spike_Type") spike_col = i;
                }
                next
            }
            $sample_col == sample && spike_col > 0 { print $spike_col; exit }' "$META_TSV")
        [[ -z "$SPIKE_TYPE" ]] && SPIKE_TYPE="none"
    else
        log "WARN" "⚠️ META_TSV file not found: $META_TSV — skipping spike-in reference logic"
    fi

local REF_PATH="$REF_GENOME"

if [[ "$SPIKE_TYPE" != "none" ]]; then
    local COMBINED_REF="$REF_FOLDER/spikein_${SPIKE_TYPE}_${REF_PREFIX}.fa"
    local RENAMED_SPIKE="${REF_FOLDER}/tmp/renamed_${SPIKE_TYPE}.fa"
    
    mkdir -p "${REF_FOLDER}/tmp"

    if [[ ! -f "$COMBINED_REF" ]]; then
        log "INFO" "Combined reference not found — generating: $COMBINED_REF"

        log "INFO" "Renaming spike-in contigs to avoid conflicts..."
        awk '{if ($0 ~ /^>/) print ">spike_" substr($0,2); else print $0}' \
            "SpikeinReference/${SPIKE_TYPE}/${SPIKE_TYPE}.fa" > "$RENAMED_SPIKE" || {
            log "ERROR" "❌ Failed to rename spike-in contigs"
            return 1
        }

        cat "${REF_FOLDER}/${REF_PREFIX}.fa" "$RENAMED_SPIKE" > "$COMBINED_REF" || {
            log "ERROR" "❌ Failed to create combined reference FASTA"
            return 1
        }
    fi

    if [[ ! -f "${COMBINED_REF}.bwt" ]]; then
        log "INFO" "BWA index not found for $COMBINED_REF — generating index..."
        bwa index "$COMBINED_REF" || {
            log "ERROR" "❌ Failed to index combined reference: $COMBINED_REF"
            return 1
        }
    fi

    REF_PATH="$COMBINED_REF"
    log "INFO" "✅ Using combined reference for $base: $COMBINED_REF"
else
    if [[ ! -f "${REF_PATH}.bwt" ]]; then
        log "WARN" "⚠️ BWA index not found for $REF_PATH — generating index..."
        bwa index "$REF_PATH" || {
            log "ERROR" "❌ Failed to index reference: $REF_PATH"
            return 1
        }
    fi
fi

# Check if SAM already exists

    log "INFO" "🧩 Starting alignment..."
    start_monitoring "alignment" $$
    bwa mem -M -t "$THREADS" "$REF_PATH" "$TRIM1" "$TRIM2" > "${BAM_DIR}/${base}.sam" || {
        log "ERROR" "❌ Alignment failed for sample $base"
        return 1
    }
    stop_monitoring "alignment"


    log "INFO" "📦 Starting sorting..."
    start_monitoring "sorting" $$
    samtools view -bS "${BAM_DIR}/${base}.sam" | \
        samtools sort -@ "$((THREADS/2))" -o "${BAM_DIR}/${base}.sorted.bam" -T "${BAM_DIR}/${base}.sort.tmp" || {
        log "ERROR" "❌ Sorting failed for sample $base"
        return 1
    }
    stop_monitoring "sorting"

    log "INFO" "🏷️ Starting indexing..."
    start_monitoring "indexing" $$
    samtools index "${BAM_DIR}/${base}.sorted.bam" || {
        log "ERROR" "❌ Indexing failed for sample $base"
        return 1
    }
    stop_monitoring "indexing"

    validate_bam "${BAM_DIR}/${base}.sorted.bam" || {
        log "ERROR" "❌ Generated BAM failed validation"
        return 1
    }

    rm "${BAM_DIR}/${base}.sam"
}

# === Main Script ===
main() {
    check_requirements
    validate_inputs
    check_disk_space
    validate_reference
    # Clean old logs (now safely inside main execution block)
    
    # Ensure necessary directories exist
    if [[ ! -d "$INPUT_DIR" ]]; then
        log "ERROR" "❌ Provided input directory does not exist: $INPUT_DIR"
    exit 1
    fi


    mkdir -p  "$LOG_DIR" "$BAM_DIR" "$FILTERED_DIR" "$DEDUP_DIR"
    
    # Clean old logs
    find "$LOG_DIR" -name "*performance_*.log" -mtime +30 -delete
    find "$LOG_DIR" -name "*system_metrics_*.csv" -mtime +30 -delete   
    
echo "=============================================="
log "INFO" "Checking required tools..."
log "INFO" "Checking disk space availability..."
log "INFO" "🧬 ${SCRIPT_NAME} v${VERSION}"
log "INFO" "📌 Reference: $REF_GENOME"
log "INFO" "💡 Threads: $THREADS"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "🕒 Start time: $(date '+%F %T')"
echo "=============================================="

if "$DRY_RUN"; then
    log "INFO" "🧯 SIMULATION RESULTS:"
    log "INFO" "🌀 Input FASTQs that would be processed:"
    ls "${INPUT_DIR}/"*_1_trimmed.fastq.gz | head -3
    [[ $(ls "${INPUT_DIR}/"*_1_trimmed.fastq.gz | wc -l) -gt 3 ]] && echo "..."
    
    log "INFO" "📝 Example commands that would run:"
    echo "bwa mem -M -t $THREADS $REF_GENOME input_R1.fq input_R2.fq > output.sam"
    echo "samtools view -bS output.sam | samtools sort -@ $((THREADS/2)) -o output.sorted.bam"
    echo "samtools index output.sorted.bam"
    
    log "INFO" "🧯 DRY-RUN completed successfully"
    exit 0
fi

# Process samples
processed=0
skipped=0
failed=0

log "INFO" "Scanning input: ${INPUT_DIR}/ *_1_trimmed.fastq.gz"
ls -1 "${INPUT_DIR}/"*_1_trimmed.fastq.gz
for TRIM1 in "${INPUT_DIR}/"*_1_trimmed.fastq.gz; do
    base=$(basename "$TRIM1" _1_trimmed.fastq.gz)
    TRIM2="${INPUT_DIR}/${base}_2_trimmed.fastq.gz"

    validate_sample_pair "$base" "$TRIM1" "$TRIM2" || continue

    if "$DRY_RUN"; then
        log "INFO" "🧯Would process sample: $base"
        continue
    fi

    log "INFO" "👻 Invoking align_sample for: $base"
    if align_sample "$base" "$TRIM1" "$TRIM2"; then
        ((processed++))
    else
        ((skipped++))
    fi
done

    # === Completion Message ===
    echo "=============================================="
    log "INFO" "✅ ALIGNMENT COMPLETED"
    log "INFO" "📦 Newly processed samples: $processed"
    log "INFO" "⏭️ Skipped samples: $skipped"
    log "INFO" "📂 Output directory: $BAM_DIR"
    log "INFO" "🕒 End time: $(date '+%F %T')"
    echo "=============================================="

    
   # Generate reports (FIXED)
generate_report

# Always write system metrics if collected
[[ -n "${SYSTEM_METRICS[*]}" ]] && write_system_metrics_csv

# Only log file paths if files exist
[[ -f "$PERF_LOG" ]] && log "INFO" "Performance report: $PERF_LOG"
[[ -f "$SYSTEM_METRICS_LOG" ]] && log "INFO" "System metrics: $SYSTEM_METRICS_LOG"

#check_disk_space 
   } 
# Generate reports
# generate_report
# [[ "$REPORT_FORMAT" != "human" ]] && write_system_metrics_csv
#log "INFO" "Performance report saved to $PERF_LOG"
#    [[ "$REPORT_FORMAT" != "human" ]] && log "INFO" "System metrics saved to $SYSTEM_METRICS_LOG"


# === Entry Point ===

parse_args "$@"

# === Set default if --meta-tsv was not passed ===
META_TSV=${META_TSV:-metadata/mapping.tsv}

# Ensure necessary directories exist
if [[ ! -d "$INPUT_DIR" ]]; then
    log "ERROR" "Provided input directory does not exist: $INPUT_DIR"
    exit 1
fi

mkdir -p logs "$BAM_DIR" "$FILTERED_DIR" "$DEDUP_DIR"

# Register cleanup on script exit
trap cleanup EXIT

# Launch main logic
main | tee -a "$MODULE_LOG"


