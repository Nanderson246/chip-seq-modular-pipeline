#!/bin/bash
# Module: #!/bin/bash
# Module: 09_readgroups_add.sh
# Author: Nancy Anderson
# Description: Add readgroups and deduplicate BAM files with resource monitoring
################################################################################
# SOFTWARE REQUIREMENTS
#
# This script depends on the following tools being available in your system or
# in the 'tools/' directory (for wrapped versions):
#
# Required Tools:
#   • bash       - Unix shell interpreter (version ≥ 4)
#   • java       - Java Runtime Environment (version ≥ 8)
#   • samtools   - SAM/BAM file handling (version ≥ 1.9)
#   • picard     - Toolkit for BAM manipulation (picard.jar ≥ 2.23.0)
#   • qualimap   - BAM QC and summary metrics (version ≥ 2.2.1)
#
# Optional (for system monitoring and reporting):
#   • awk        - Pattern scanning & processing language
#   • top        - CPU monitoring (usually pre-installed)
#   • free       - Memory usage tool (part of procps-ng)
#   • iostat     - Disk I/O statistics (from 'sysstat' package)
#
# Installation Tips (Debian/Ubuntu):
#   sudo apt install samtools default-jre procps sysstat gawk
#
# Wrapper Support:
#   Place 'picard.jar' in tools/
#   Place Qualimap binary in tools/qualimap/qualimap and ensure it is executable
#
# Example:
#   chmod +x tools/qualimap/qualimap
################################################################################

################################################################################
# USAGE:
#   Pipeline
#      [options] [PLATFORM]
# bash modules/pipeline1/09_readgroups_add.sh [PLATFORM]
#default: bash modules/pipeline1/09_readgroups_add.sh ILLUMINA
# 
#Standalone usage:
#   bash 09_readgroups_add.sh --in <bams_folder> --out <output_folder> [PLATFORM]
# NEW FEATURES:
#   - CPU/Memory/Disk monitoring
#   - Nanosecond timing
#   - CSV/TSV performance reports
#   - System-wide resource logging
################################################################################
set -uo pipefail

: "${TOOLS_DIR:=tools}"

# === Load wrappers if not already loaded ===

echo "[INFO] Checking for local tool wrappers..."

# === Helper to resolve project root by climbing directory tree ===
get_root_dir() {
    local dir
    dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
    while [[ "$dir" != "/" ]]; do
        if [[ -d "$dir/tools" ]]; then
            echo "$dir"
            return
        fi
        dir="$(dirname "$dir")"
    done
    echo "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
}

# Only define if not already set (support for pipeline-wide sourcing)
if [[ -z "$TOOLS_DIR" ]]; then
    ROOT_DIR="$(get_root_dir)"
    TOOLS_DIR="$ROOT_DIR/tools"
    echo "[DEBUG] 🔧 TOOLS_DIR resolved to: $TOOLS_DIR"
fi

# === Picard wrapper ===
if ! command -v picard >/dev/null 2>&1; then
    picard() {
        jar="$TOOLS_DIR/picard.jar"
        echo "[DEBUG] Looking for picard.jar at: $jar"
        if [[ ! -f "$jar" ]]; then
            echo "[ERROR]❌  picard.jar not found at $jar"
            exit 1
        fi
        java -jar "$jar" "$@"
    }
    export -f picard
fi

# === Qualimap wrapper ===
if ! command -v qualimap >/dev/null 2>&1; then
    qualimap() {
        bin="$TOOLS_DIR/qualimap/qualimap"
        echo "[DEBUG] Looking for qualimap at: $bin"
        if [[ ! -f "$bin" ]]; then
            echo "[ERROR]❌  qualimap not found at $bin"
            exit 1
        fi
        if [[ ! -x "$bin" ]]; then
            echo "[ERROR] ❌ qualimap is not executable at $bin"
            exit 1
        fi
        "$bin" "$@"
    }
    export -f qualimap
fi

    
echo "[INFO] 🔧Tool wrappers loaded successfully."
echo "[INFO] ✅Using system-wide tool wrappers."


# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}"  
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Logging Configuration ===
LOG_DIR="logs/pipeline1"
readonly PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"
readonly SYSTEM_METRICS_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_metrics_${TIMESTAMP}.csv"
readonly MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
mkdir -p "$LOG_DIR"

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

check_and_log_tool_versions() {
    local tools=("samtools" "java")

    log "INFO" "🔧 Verifying and logging required tool versions..."

    for tool in "${tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            log "ERROR" "❌ Required tool not found in PATH: $tool"
            exit 1
        fi
    done

    log "INFO" "• samtools:     $(samtools --version 2>&1 | head -n 1)"
    log "INFO" "• java:         $(java -version 2>&1 | head -n 1)"
    log "INFO" "• picard:       $(java -jar ${TOOLS_DIR}/picard.jar MarkDuplicates --version 2>&1 | grep -Ei 'version|MarkDuplicates' | head -n 1 || echo 'unknown')"
    log "INFO" "• qualimap:     $(${TOOLS_DIR}/qualimap/qualimap --version 2>&1 | head -n 1 || echo 'unknown')"
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
        log "WARN" "⚠️iostat not available — skipping disk IO metrics"
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

generate_performance_report() {
    echo -e "\n${BLUE}=== PERFORMANCE REPORT ===${NC}"
    printf "%-15s %10s %8s %10s\n" "Stage" "Time(s)" "CPU(%)" "Mem(MB)"
    echo "--------------------------------------------"
    
    local total_time=0
    for stage in readgroup deduplication indexing_rg indexing_dedup qualimap; do
        local stage_time=${PROCESS_METRICS["${stage}_time"]:-0}
        printf "%-15s %10.2f %8.1f %10.1f\n" \
            "$stage" \
            "$stage_time" \
            "${PROCESS_METRICS["${stage}_cpu"]:-0}" \
            "${PROCESS_METRICS["${stage}_mem"]:-0}"
        total_time=$(awk -v t="$total_time" -v s="$stage_time" 'BEGIN {print t + s}')
    done
    
    echo "--------------------------------------------"
    printf "%-15s %10.2f\n" "TOTAL" "$total_time"
    echo -e "${BLUE}============================================${NC}"
}

write_system_metrics_csv() {
    echo "timestamp,stage,event,cpu_percent,mem_used_mb,mem_percent,disk_read_mb,disk_write_mb" > "$SYSTEM_METRICS_LOG"
    for stage in readgroup deduplication indexing_rg indexing_dedup qualimap; do
        for point in start end; do
            echo "${SYSTEM_METRICS["${stage}_${point}_timestamp"]},$stage,$point,"\
                 "${SYSTEM_METRICS["${stage}_${point}_cpu"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_mem_used"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_mem_pct"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_disk_read"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_disk_write"]}" >> "$SYSTEM_METRICS_LOG"
        done
    done
}

# === Main Script ===

# === Default directories ===
INPUT_DIR="results/BAM"
OUTPUT_DIR="results/Filtered/Deduplicated"


# Initialize directories
VALID_PLATFORMS=("ILLUMINA" "IONTORRENT" "PACBIO" "ONT" "SOLID" "NANOPORE")
BAM_DIR="$INPUT_DIR"
DEDUP_DIR="$OUTPUT_DIR"
FILTERED_DIR=$(dirname "$OUTPUT_DIR")
METRICS_DIR="${FILTERED_DIR}/Metrics"

mkdir -p  "$BAM_DIR" "$FILTERED_DIR" "$DEDUP_DIR" "$METRICS_DIR"

# === CLI Parser ===
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            INPUT_DIR="$2"
            shift 2
            ;;
        --out)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "${SCRIPT_NAME} v${VERSION}"
            echo "Usage: $SCRIPT_NAME [--in DIR] [--out DIR] [--dry-run] [PLATFORM]"
            echo "Platforms: ${VALID_PLATFORMS[*]}"
            exit 0
            ;;
        -v|--version)
            echo "${SCRIPT_NAME} v${VERSION}"
            exit 0
            ;;
        -*)
            log "ERROR" "Unknown option: $1"
            exit 1
            ;;
        *)
            PLATFORM="${1^^}"
            shift
            ;;
    esac
done

# === Platform Validation ===
if [[ -z "$PLATFORM" ]]; then
    log "ERROR" "❌ PLATFORM not specified. Valid options: ${VALID_PLATFORMS[*]}"
    exit 1
fi

if [[ ! " ${VALID_PLATFORMS[*]} " =~ " ${PLATFORM} " ]]; then
    log "ERROR" "❌ Invalid PLATFORM: $PLATFORM. Valid options: ${VALID_PLATFORMS[*]}"
    exit 1
fi

{
echo "=============================================="
log "INFO" "🧬 ${SCRIPT_NAME} v${VERSION}"
log "INFO" "🔧 Platform: $PLATFORM"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "📥 Input directory: $INPUT_DIR"
log "INFO" "🕒 Start time: $(date '+%F %T')"
echo "=============================================="

# === Tools verification and version ===

check_and_log_tool_versions

# === Dry-run Simulation ===
if "$DRY_RUN"; then
    log "INFO" "\n🔍 SIMULATION RESULTS:"
    log "INFO" "BAM files that would be processed:"
    ls "${BAM_DIR}/"*.sorted.bam | head -3
    [[ $(ls "${BAM_DIR}/"*.sorted.bam | wc -l) -gt 3 ]] && echo "..."
    
    log "INFO" "Example commands that would run:"
    echo "picard AddOrReplaceReadGroups I=input.bam O=output.rg.bam \\"
    echo "  RGID=sample1 RGLB=lib1 RGPL=$PLATFORM RGPU=unit1 RGSM=sample1"
    echo "picard MarkDuplicates I=input.rg.bam O=output.dedup.bam \\"
    echo "  M=metrics.txt REMOVE_DUPLICATES=true"
    
    log "INFO" "DRY-RUN completed successfully"
    exit 0
fi

# === Main Processing ===

processed=0
skipped=0

for SORTBAM in "${BAM_DIR}/"*.sorted.bam; do
    [ -e "$SORTBAM" ] || continue
    base=$(basename "$SORTBAM" .sorted.bam)
    RG_BAM="${BAM_DIR}/${base}.rg.bam"
    DEDUP_BAM="${DEDUP_DIR}/${base}.dedup.bam"
    
    # Skip if already processed
    if [[ -f "$DEDUP_BAM" ]] && [[ -f "${DEDUP_BAM}.bai" ]]; then
        log "INFO" "Skipping completed sample: $base"
        ((skipped++))
        continue
    fi
    
    log "INFO" "⚙️ Processing sample $base..."
    
    # Read Group Handling
    if ! samtools view -H "$SORTBAM" | grep -q '@RG'; then
        log "INFO" "📇 Adding read groups..."
        picard AddOrReplaceReadGroups \
            I="$SORTBAM" \
            O="$RG_BAM" \
            RGID="$base" \
            RGLB=lib1 \
            RGPL="$PLATFORM" \
            RGPU=unit1 \
            RGSM="$base" &
        rg_pid=$!
        start_monitoring "readgroup" $rg_pid
        wait $rg_pid || {
            log "ERROR" "❌ Read group addition failed for $base"
            exit 1
        }
        stop_monitoring "readgroup"
        
        samtools index "$RG_BAM" &
        index_pid=$!
        start_monitoring "indexing_rg" $index_pid
        wait $index_pid || {
            log "ERROR" "❌ Indexing failed for $base"
            exit 1
        }
        stop_monitoring "indexing_rg"
        
        INPUT_BAM="$RG_BAM"
    else
        log "INFO" "Read groups already present"
        INPUT_BAM="$SORTBAM"
    fi
    
    # Deduplication
    log "INFO" "🔁 Marking duplicates..."
    picard MarkDuplicates \
        I="$INPUT_BAM" \
        O="${BAM_DIR}/${base}.dedup.bam" \
        M="${METRICS_DIR}/${base}_metrics.txt" \
        REMOVE_DUPLICATES=true &
    dedup_pid=$!
    start_monitoring "deduplication" $dedup_pid
    wait $dedup_pid || {
        log "ERROR" "❌ Deduplication failed for $base"
        exit 1
    }
    stop_monitoring "deduplication"
    
    # Post-processing
    log "INFO" "📈 Generating QC metrics..."
    samtools index "${BAM_DIR}/${base}.dedup.bam" &
    index_pid=$!
    start_monitoring "indexing_dedup" $index_pid
    wait $index_pid || {
        log "ERROR" "❌ Indexing failed for $base"
        exit 1
    }
    stop_monitoring "indexing_dedup"
    
    QUALIMAP_OUTDIR="${METRICS_DIR}/${base}_qualimap_dedup_$(date +%s)"
    
    qualimap bamqc \
        -bam "${BAM_DIR}/${base}.dedup.bam" \
        -outdir "$QUALIMAP_OUTDIR" \
        -outformat PDF:HTML &
    qc_pid=$!
    start_monitoring "qualimap" $qc_pid
    wait $qc_pid || {
        log "ERROR" "❌ Qualimap failed for $base"
        exit 1
    }
    stop_monitoring "qualimap"
    
    # Final organization
    log "INFO" "🗂️ Moving to final location..."
    mv "${BAM_DIR}/${base}.dedup.bam" "$DEDUP_DIR/"
    mv "${BAM_DIR}/${base}.dedup.bam.bai" "$DEDUP_DIR/"
    
    ((processed++))
done

# Generate performance report
generate_performance_report | tee -a "$PERF_LOG"

# Write detailed system metrics if desired for analysis
[[ -n "${SYSTEM_METRICS[*]}" ]] && write_system_metrics_csv


# === Completion ===
echo "=============================================="
log "INFO" "✅ PROCESSING COMPLETED"
log "INFO" "🧪 Processed samples: $processed"
log "INFO" "⏭️ Skipped samples: $skipped"
log "INFO" "📂 Output directory: $DEDUP_DIR"
log "INFO" "📊 Performance log: $PERF_LOG"
log "INFO" "🕒 End time: $(date '+%F %T')"
echo "=============================================="
} | tee -a "$MODULE_LOG"
#Standalone usage:
#   bash 09_readgroups_add.sh --in <bams_folder> --out <output_folder> [PLATFORM]
# NEW FEATURES:
#   - CPU/Memory/Disk monitoring
#   - Nanosecond timing
#   - CSV/TSV performance reports
#   - System-wide resource logging
################################################################################
set -uo pipefail

: "${TOOLS_DIR:=tools}"

# === Load wrappers if not already loaded ===

echo "[INFO] Checking for local tool wrappers..."

# === Helper to resolve project root by climbing directory tree ===
get_root_dir() {
    local dir
    dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
    while [[ "$dir" != "/" ]]; do
        if [[ -d "$dir/tools" ]]; then
            echo "$dir"
            return
        fi
        dir="$(dirname "$dir")"
    done
    echo "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
}

# Only define if not already set (support for pipeline-wide sourcing)
if [[ -z "$TOOLS_DIR" ]]; then
    ROOT_DIR="$(get_root_dir)"
    TOOLS_DIR="$ROOT_DIR/tools"
    echo "[DEBUG] 🔧 TOOLS_DIR resolved to: $TOOLS_DIR"
fi

# === Picard wrapper ===
if ! command -v picard >/dev/null 2>&1; then
    picard() {
        jar="$TOOLS_DIR/picard.jar"
        echo "[DEBUG] Looking for picard.jar at: $jar"
        if [[ ! -f "$jar" ]]; then
            echo "[ERROR]❌  picard.jar not found at $jar"
            exit 1
        fi
        java -jar "$jar" "$@"
    }
    export -f picard
fi

# === Qualimap wrapper ===
if ! command -v qualimap >/dev/null 2>&1; then
    qualimap() {
        bin="$TOOLS_DIR/qualimap/qualimap"
        echo "[DEBUG] Looking for qualimap at: $bin"
        if [[ ! -f "$bin" ]]; then
            echo "[ERROR]❌  qualimap not found at $bin"
            exit 1
        fi
        if [[ ! -x "$bin" ]]; then
            echo "[ERROR] ❌ qualimap is not executable at $bin"
            exit 1
        fi
        "$bin" "$@"
    }
    export -f qualimap
fi

    
echo "[INFO] 🔧Tool wrappers loaded successfully."
echo "[INFO] ✅Using system-wide tool wrappers."


# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}"  
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Logging Configuration ===
LOG_DIR="logs"
readonly PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"
readonly SYSTEM_METRICS_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_metrics_${TIMESTAMP}.csv"
readonly MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
mkdir -p "$LOG_DIR"

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

check_and_log_tool_versions() {
    local tools=("samtools" "java")

    log "INFO" "🔧 Verifying and logging required tool versions..."

    for tool in "${tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            log "ERROR" "❌ Required tool not found in PATH: $tool"
            exit 1
        fi
    done

    log "INFO" "• samtools:     $(samtools --version 2>&1 | head -n 1)"
    log "INFO" "• java:         $(java -version 2>&1 | head -n 1)"
    log "INFO" "• picard:       $(java -jar ${TOOLS_DIR}/picard.jar MarkDuplicates --version 2>&1 | grep -Ei 'version|MarkDuplicates' | head -n 1 || echo 'unknown')"
    log "INFO" "• qualimap:     $(${TOOLS_DIR}/qualimap/qualimap --version 2>&1 | head -n 1 || echo 'unknown')"
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
        log "WARN" "⚠️iostat not available — skipping disk IO metrics"
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

generate_performance_report() {
    echo -e "\n${BLUE}=== PERFORMANCE REPORT ===${NC}"
    printf "%-15s %10s %8s %10s\n" "Stage" "Time(s)" "CPU(%)" "Mem(MB)"
    echo "--------------------------------------------"
    
    local total_time=0
    for stage in readgroup deduplication indexing_rg indexing_dedup qualimap; do
        local stage_time=${PROCESS_METRICS["${stage}_time"]:-0}
        printf "%-15s %10.2f %8.1f %10.1f\n" \
            "$stage" \
            "$stage_time" \
            "${PROCESS_METRICS["${stage}_cpu"]:-0}" \
            "${PROCESS_METRICS["${stage}_mem"]:-0}"
        total_time=$(awk -v t="$total_time" -v s="$stage_time" 'BEGIN {print t + s}')
    done
    
    echo "--------------------------------------------"
    printf "%-15s %10.2f\n" "TOTAL" "$total_time"
    echo -e "${BLUE}============================================${NC}"
}

write_system_metrics_csv() {
    echo "timestamp,stage,event,cpu_percent,mem_used_mb,mem_percent,disk_read_mb,disk_write_mb" > "$SYSTEM_METRICS_LOG"
    for stage in readgroup deduplication indexing_rg indexing_dedup qualimap; do
        for point in start end; do
            echo "${SYSTEM_METRICS["${stage}_${point}_timestamp"]},$stage,$point,"\
                 "${SYSTEM_METRICS["${stage}_${point}_cpu"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_mem_used"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_mem_pct"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_disk_read"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_disk_write"]}" >> "$SYSTEM_METRICS_LOG"
        done
    done
}

# === Main Script ===

# === Default directories ===
INPUT_DIR="results/BAM"
OUTPUT_DIR="results/Filtered/Deduplicated"


# Initialize directories
VALID_PLATFORMS=("ILLUMINA" "IONTORRENT" "PACBIO" "ONT" "SOLID" "NANOPORE")
BAM_DIR="$INPUT_DIR"
DEDUP_DIR="$OUTPUT_DIR"
FILTERED_DIR=$(dirname "$OUTPUT_DIR")
METRICS_DIR="${FILTERED_DIR}/Metrics"

mkdir -p  "$BAM_DIR" "$FILTERED_DIR" "$DEDUP_DIR" "$METRICS_DIR"

# === CLI Parser ===
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            INPUT_DIR="$2"
            shift 2
            ;;
        --out)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "${SCRIPT_NAME} v${VERSION}"
            echo "Usage: $SCRIPT_NAME [--in DIR] [--out DIR] [--dry-run] [PLATFORM]"
            echo "Platforms: ${VALID_PLATFORMS[*]}"
            exit 0
            ;;
        -v|--version)
            echo "${SCRIPT_NAME} v${VERSION}"
            exit 0
            ;;
        -*)
            log "ERROR" "Unknown option: $1"
            exit 1
            ;;
        *)
            PLATFORM="${1^^}"
            shift
            ;;
    esac
done

# === Platform Validation ===
if [[ -z "$PLATFORM" ]]; then
    log "ERROR" "❌ PLATFORM not specified. Valid options: ${VALID_PLATFORMS[*]}"
    exit 1
fi

if [[ ! " ${VALID_PLATFORMS[*]} " =~ " ${PLATFORM} " ]]; then
    log "ERROR" "❌ Invalid PLATFORM: $PLATFORM. Valid options: ${VALID_PLATFORMS[*]}"
    exit 1
fi

{
echo "=============================================="
log "INFO" "🧬 ${SCRIPT_NAME} v${VERSION}"
log "INFO" "🔧 Platform: $PLATFORM"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "📥 Input directory: $INPUT_DIR"
log "INFO" "🕒 Start time: $(date '+%F %T')"
echo "=============================================="

# === Tools verification and version ===

check_and_log_tool_versions

# === Dry-run Simulation ===
if "$DRY_RUN"; then
    log "INFO" "\n🔍 SIMULATION RESULTS:"
    log "INFO" "BAM files that would be processed:"
    ls "${BAM_DIR}/"*.sorted.bam | head -3
    [[ $(ls "${BAM_DIR}/"*.sorted.bam | wc -l) -gt 3 ]] && echo "..."
    
    log "INFO" "Example commands that would run:"
    echo "picard AddOrReplaceReadGroups I=input.bam O=output.rg.bam \\"
    echo "  RGID=sample1 RGLB=lib1 RGPL=$PLATFORM RGPU=unit1 RGSM=sample1"
    echo "picard MarkDuplicates I=input.rg.bam O=output.dedup.bam \\"
    echo "  M=metrics.txt REMOVE_DUPLICATES=true"
    
    log "INFO" "DRY-RUN completed successfully"
    exit 0
fi

# === Main Processing ===

processed=0
skipped=0

for SORTBAM in "${BAM_DIR}/"*.sorted.bam; do
    [ -e "$SORTBAM" ] || continue
    base=$(basename "$SORTBAM" .sorted.bam)
    RG_BAM="${BAM_DIR}/${base}.rg.bam"
    DEDUP_BAM="${DEDUP_DIR}/${base}.dedup.bam"
    
    # Skip if already processed
    if [[ -f "$DEDUP_BAM" ]] && [[ -f "${DEDUP_BAM}.bai" ]]; then
        log "INFO" "Skipping completed sample: $base"
        ((skipped++))
        continue
    fi
    
    log "INFO" "⚙️ Processing sample $base..."
    
    # Read Group Handling
    if ! samtools view -H "$SORTBAM" | grep -q '@RG'; then
        log "INFO" "📇 Adding read groups..."
        picard AddOrReplaceReadGroups \
            I="$SORTBAM" \
            O="$RG_BAM" \
            RGID="$base" \
            RGLB=lib1 \
            RGPL="$PLATFORM" \
            RGPU=unit1 \
            RGSM="$base" &
        rg_pid=$!
        start_monitoring "readgroup" $rg_pid
        wait $rg_pid || {
            log "ERROR" "❌ Read group addition failed for $base"
            exit 1
        }
        stop_monitoring "readgroup"
        
        samtools index "$RG_BAM" &
        index_pid=$!
        start_monitoring "indexing_rg" $index_pid
        wait $index_pid || {
            log "ERROR" "❌ Indexing failed for $base"
            exit 1
        }
        stop_monitoring "indexing_rg"
        
        INPUT_BAM="$RG_BAM"
    else
        log "INFO" "Read groups already present"
        INPUT_BAM="$SORTBAM"
    fi
    
    # Deduplication
    log "INFO" "🔁 Marking duplicates..."
    picard MarkDuplicates \
        I="$INPUT_BAM" \
        O="${BAM_DIR}/${base}.dedup.bam" \
        M="${METRICS_DIR}/${base}_metrics.txt" \
        REMOVE_DUPLICATES=true &
    dedup_pid=$!
    start_monitoring "deduplication" $dedup_pid
    wait $dedup_pid || {
        log "ERROR" "❌ Deduplication failed for $base"
        exit 1
    }
    stop_monitoring "deduplication"
    
    # Post-processing
    log "INFO" "📈 Generating QC metrics..."
    samtools index "${BAM_DIR}/${base}.dedup.bam" &
    index_pid=$!
    start_monitoring "indexing_dedup" $index_pid
    wait $index_pid || {
        log "ERROR" "❌ Indexing failed for $base"
        exit 1
    }
    stop_monitoring "indexing_dedup"
    
    QUALIMAP_OUTDIR="${METRICS_DIR}/${base}_qualimap_dedup_$(date +%s)"
    
    qualimap bamqc \
        -bam "${BAM_DIR}/${base}.dedup.bam" \
        -outdir "$QUALIMAP_OUTDIR" \
        -outformat PDF:HTML &
    qc_pid=$!
    start_monitoring "qualimap" $qc_pid
    wait $qc_pid || {
        log "ERROR" "❌ Qualimap failed for $base"
        exit 1
    }
    stop_monitoring "qualimap"
    
    # Final organization
    log "INFO" "🗂️ Moving to final location..."
    mv "${BAM_DIR}/${base}.dedup.bam" "$DEDUP_DIR/"
    mv "${BAM_DIR}/${base}.dedup.bam.bai" "$DEDUP_DIR/"
    
    ((processed++))
done

# Generate performance report
generate_performance_report | tee -a "$PERF_LOG"

# Write detailed system metrics if desired for analysis
[[ -n "${SYSTEM_METRICS[*]}" ]] && write_system_metrics_csv


# === Completion ===
echo "=============================================="
log "INFO" "✅ PROCESSING COMPLETED"
log "INFO" "🧪 Processed samples: $processed"
log "INFO" "⏭️ Skipped samples: $skipped"
log "INFO" "📂 Output directory: $DEDUP_DIR"
log "INFO" "📊 Performance log: $PERF_LOG"
log "INFO" "🕒 End time: $(date '+%F %T')"
echo "=============================================="
} | tee -a "$MODULE_LOG"
# Author: Nancy Anderson
# Description: Add readgroups and deduplicate BAM files with resource monitoring
################################################################################
# USAGE:
#   Pipeline
#      [options] [PLATFORM]
#
#Standalone usage:
#   bash 09_readgroups_add.sh --in <bams_folder> --out <output_folder> [PLATFORM]
# NEW FEATURES:
#   - CPU/Memory/Disk monitoring
#   - Nanosecond timing
#   - CSV/TSV performance reports
#   - System-wide resource logging
################################################################################
set -uo pipefail

: "${TOOLS_DIR:=tools}"

# === Load wrappers if not already loaded ===

echo "[INFO] Checking for local tool wrappers..."

# === Helper to resolve project root by climbing directory tree ===
get_root_dir() {
    local dir
    dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
    while [[ "$dir" != "/" ]]; do
        if [[ -d "$dir/tools" ]]; then
            echo "$dir"
            return
        fi
        dir="$(dirname "$dir")"
    done
    echo "$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
}

# Only define if not already set (support for pipeline-wide sourcing)
if [[ -z "$TOOLS_DIR" ]]; then
    ROOT_DIR="$(get_root_dir)"
    TOOLS_DIR="$ROOT_DIR/tools"
    echo "[DEBUG] 🔧 TOOLS_DIR resolved to: $TOOLS_DIR"
fi

# === Picard wrapper ===
if ! command -v picard >/dev/null 2>&1; then
    picard() {
        jar="$TOOLS_DIR/picard.jar"
        echo "[DEBUG] Looking for picard.jar at: $jar"
        if [[ ! -f "$jar" ]]; then
            echo "[ERROR]❌  picard.jar not found at $jar"
            exit 1
        fi
        java -jar "$jar" "$@"
    }
    export -f picard
fi

# === Qualimap wrapper ===
if ! command -v qualimap >/dev/null 2>&1; then
    qualimap() {
        bin="$TOOLS_DIR/qualimap/qualimap"
        echo "[DEBUG] Looking for qualimap at: $bin"
        if [[ ! -f "$bin" ]]; then
            echo "[ERROR]❌  qualimap not found at $bin"
            exit 1
        fi
        if [[ ! -x "$bin" ]]; then
            echo "[ERROR] ❌ qualimap is not executable at $bin"
            exit 1
        fi
        "$bin" "$@"
    }
    export -f qualimap
fi

    
echo "[INFO] 🔧Tool wrappers loaded successfully."
echo "[INFO] ✅Using system-wide tool wrappers."


# === Constants ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}"  
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

# === Logging Configuration ===
LOG_DIR="logs"
readonly PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"
readonly SYSTEM_METRICS_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_metrics_${TIMESTAMP}.csv"
readonly MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
mkdir -p "$LOG_DIR"

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

check_and_log_tool_versions() {
    local tools=("samtools" "java")

    log "INFO" "🔧 Verifying and logging required tool versions..."

    for tool in "${tools[@]}"; do
        if ! command -v "$tool" >/dev/null 2>&1; then
            log "ERROR" "❌ Required tool not found in PATH: $tool"
            exit 1
        fi
    done

    log "INFO" "• samtools:     $(samtools --version 2>&1 | head -n 1)"
    log "INFO" "• java:         $(java -version 2>&1 | head -n 1)"
    log "INFO" "• picard:       $(java -jar ${TOOLS_DIR}/picard.jar MarkDuplicates --version 2>&1 | grep -Ei 'version|MarkDuplicates' | head -n 1 || echo 'unknown')"
    log "INFO" "• qualimap:     $(${TOOLS_DIR}/qualimap/qualimap --version 2>&1 | head -n 1 || echo 'unknown')"
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
        log "WARN" "⚠️iostat not available — skipping disk IO metrics"
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

generate_performance_report() {
    echo -e "\n${BLUE}=== PERFORMANCE REPORT ===${NC}"
    printf "%-15s %10s %8s %10s\n" "Stage" "Time(s)" "CPU(%)" "Mem(MB)"
    echo "--------------------------------------------"
    
    local total_time=0
    for stage in readgroup deduplication indexing_rg indexing_dedup qualimap; do
        local stage_time=${PROCESS_METRICS["${stage}_time"]:-0}
        printf "%-15s %10.2f %8.1f %10.1f\n" \
            "$stage" \
            "$stage_time" \
            "${PROCESS_METRICS["${stage}_cpu"]:-0}" \
            "${PROCESS_METRICS["${stage}_mem"]:-0}"
        total_time=$(awk -v t="$total_time" -v s="$stage_time" 'BEGIN {print t + s}')
    done
    
    echo "--------------------------------------------"
    printf "%-15s %10.2f\n" "TOTAL" "$total_time"
    echo -e "${BLUE}============================================${NC}"
}

write_system_metrics_csv() {
    echo "timestamp,stage,event,cpu_percent,mem_used_mb,mem_percent,disk_read_mb,disk_write_mb" > "$SYSTEM_METRICS_LOG"
    for stage in readgroup deduplication indexing_rg indexing_dedup qualimap; do
        for point in start end; do
            echo "${SYSTEM_METRICS["${stage}_${point}_timestamp"]},$stage,$point,"\
                 "${SYSTEM_METRICS["${stage}_${point}_cpu"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_mem_used"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_mem_pct"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_disk_read"]},"\
                 "${SYSTEM_METRICS["${stage}_${point}_disk_write"]}" >> "$SYSTEM_METRICS_LOG"
        done
    done
}

# === Main Script ===

# === Default directories ===
INPUT_DIR="results/BAM"
OUTPUT_DIR="results/Filtered/Deduplicated"


# Initialize directories
VALID_PLATFORMS=("ILLUMINA" "IONTORRENT" "PACBIO" "ONT" "SOLID" "NANOPORE")
BAM_DIR="$INPUT_DIR"
DEDUP_DIR="$OUTPUT_DIR"
FILTERED_DIR=$(dirname "$OUTPUT_DIR")
METRICS_DIR="${FILTERED_DIR}/Metrics"

mkdir -p  "$BAM_DIR" "$FILTERED_DIR" "$DEDUP_DIR" "$METRICS_DIR"

# === CLI Parser ===
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        --in)
            INPUT_DIR="$2"
            shift 2
            ;;
        --out)
            OUTPUT_DIR="$2"
            shift 2
            ;;
        -d|--dry-run)
            DRY_RUN=true
            shift
            ;;
        -h|--help)
            echo "${SCRIPT_NAME} v${VERSION}"
            echo "Usage: $SCRIPT_NAME [--in DIR] [--out DIR] [--dry-run] [PLATFORM]"
            echo "Platforms: ${VALID_PLATFORMS[*]}"
            exit 0
            ;;
        -v|--version)
            echo "${SCRIPT_NAME} v${VERSION}"
            exit 0
            ;;
        -*)
            log "ERROR" "Unknown option: $1"
            exit 1
            ;;
        *)
            PLATFORM="${1^^}"
            shift
            ;;
    esac
done

# === Platform Validation ===
if [[ -z "$PLATFORM" ]]; then
    log "ERROR" "❌ PLATFORM not specified. Valid options: ${VALID_PLATFORMS[*]}"
    exit 1
fi

if [[ ! " ${VALID_PLATFORMS[*]} " =~ " ${PLATFORM} " ]]; then
    log "ERROR" "❌ Invalid PLATFORM: $PLATFORM. Valid options: ${VALID_PLATFORMS[*]}"
    exit 1
fi

{
echo "=============================================="
log "INFO" "🧬 ${SCRIPT_NAME} v${VERSION}"
log "INFO" "🔧 Platform: $PLATFORM"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "📥 Input directory: $INPUT_DIR"
log "INFO" "🕒 Start time: $(date '+%F %T')"
echo "=============================================="

# === Tools verification and version ===

check_and_log_tool_versions

# === Dry-run Simulation ===
if "$DRY_RUN"; then
    log "INFO" "\n🔍 SIMULATION RESULTS:"
    log "INFO" "BAM files that would be processed:"
    ls "${BAM_DIR}/"*.sorted.bam | head -3
    [[ $(ls "${BAM_DIR}/"*.sorted.bam | wc -l) -gt 3 ]] && echo "..."
    
    log "INFO" "Example commands that would run:"
    echo "picard AddOrReplaceReadGroups I=input.bam O=output.rg.bam \\"
    echo "  RGID=sample1 RGLB=lib1 RGPL=$PLATFORM RGPU=unit1 RGSM=sample1"
    echo "picard MarkDuplicates I=input.rg.bam O=output.dedup.bam \\"
    echo "  M=metrics.txt REMOVE_DUPLICATES=true"
    
    log "INFO" "DRY-RUN completed successfully"
    exit 0
fi

# === Main Processing ===

processed=0
skipped=0

for SORTBAM in "${BAM_DIR}/"*.sorted.bam; do
    [ -e "$SORTBAM" ] || continue
    base=$(basename "$SORTBAM" .sorted.bam)
    RG_BAM="${BAM_DIR}/${base}.rg.bam"
    DEDUP_BAM="${DEDUP_DIR}/${base}.dedup.bam"
    
    # Skip if already processed
    if [[ -f "$DEDUP_BAM" ]] && [[ -f "${DEDUP_BAM}.bai" ]]; then
        log "INFO" "Skipping completed sample: $base"
        ((skipped++))
        continue
    fi
    
    log "INFO" "⚙️ Processing sample $base..."
    
    # Read Group Handling
    if ! samtools view -H "$SORTBAM" | grep -q '@RG'; then
        log "INFO" "📇 Adding read groups..."
        picard AddOrReplaceReadGroups \
            I="$SORTBAM" \
            O="$RG_BAM" \
            RGID="$base" \
            RGLB=lib1 \
            RGPL="$PLATFORM" \
            RGPU=unit1 \
            RGSM="$base" &
        rg_pid=$!
        start_monitoring "readgroup" $rg_pid
        wait $rg_pid || {
            log "ERROR" "❌ Read group addition failed for $base"
            exit 1
        }
        stop_monitoring "readgroup"
        
        samtools index "$RG_BAM" &
        index_pid=$!
        start_monitoring "indexing_rg" $index_pid
        wait $index_pid || {
            log "ERROR" "❌ Indexing failed for $base"
            exit 1
        }
        stop_monitoring "indexing_rg"
        
        INPUT_BAM="$RG_BAM"
    else
        log "INFO" "Read groups already present"
        INPUT_BAM="$SORTBAM"
    fi
    
    # Deduplication
    log "INFO" "🔁 Marking duplicates..."
    picard MarkDuplicates \
        I="$INPUT_BAM" \
        O="${BAM_DIR}/${base}.dedup.bam" \
        M="${METRICS_DIR}/${base}_metrics.txt" \
        REMOVE_DUPLICATES=true &
    dedup_pid=$!
    start_monitoring "deduplication" $dedup_pid
    wait $dedup_pid || {
        log "ERROR" "❌ Deduplication failed for $base"
        exit 1
    }
    stop_monitoring "deduplication"
    
    # Post-processing
    log "INFO" "📈 Generating QC metrics..."
    samtools index "${BAM_DIR}/${base}.dedup.bam" &
    index_pid=$!
    start_monitoring "indexing_dedup" $index_pid
    wait $index_pid || {
        log "ERROR" "❌ Indexing failed for $base"
        exit 1
    }
    stop_monitoring "indexing_dedup"
    
    QUALIMAP_OUTDIR="${METRICS_DIR}/${base}_qualimap_dedup_$(date +%s)"
    
    qualimap bamqc \
        -bam "${BAM_DIR}/${base}.dedup.bam" \
        -outdir "$QUALIMAP_OUTDIR" \
        -outformat PDF:HTML &
    qc_pid=$!
    start_monitoring "qualimap" $qc_pid
    wait $qc_pid || {
        log "ERROR" "❌ Qualimap failed for $base"
        exit 1
    }
    stop_monitoring "qualimap"
    
    # Final organization
    log "INFO" "🗂️ Moving to final location..."
    mv "${BAM_DIR}/${base}.dedup.bam" "$DEDUP_DIR/"
    mv "${BAM_DIR}/${base}.dedup.bam.bai" "$DEDUP_DIR/"
    
    ((processed++))
done

# Generate performance report
generate_performance_report | tee -a "$PERF_LOG"

# Write detailed system metrics if desired for analysis
[[ -n "${SYSTEM_METRICS[*]}" ]] && write_system_metrics_csv


# === Completion ===
echo "=============================================="
log "INFO" "✅ PROCESSING COMPLETED"
log "INFO" "🧪 Processed samples: $processed"
log "INFO" "⏭️ Skipped samples: $skipped"
log "INFO" "📂 Output directory: $DEDUP_DIR"
log "INFO" "📊 Performance log: $PERF_LOG"
log "INFO" "🕒 End time: $(date '+%F %T')"
echo "=============================================="
} | tee -a "$MODULE_LOG"
