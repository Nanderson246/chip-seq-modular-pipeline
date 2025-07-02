#!/usr/bin/env bash
# Module: Merge_and_Pseudoreplicates - With deepTools Integration
# Description: Metadata-driven BAM merging with pseudoreplicates and QC visualization

# === Software Requirements ===
# Required software/tools (must be in $PATH):
#   - bash (>= 4.0)
#   - samtools (>= 1.9)
#   - deepTools (>= 3.5.0)
#       - bamCoverage
#       - plotFingerprint
#       - multiBamSummary
#       - plotCorrelation
#   - GNU core utilities (cut, awk, find, etc.)
# Optional (but recommended):
#   - tee, parallel (if adapted to multithreading)
#   - ps (for CPU/memory tracking)
# Note:
#   - This script runs single-threaded by default.
#   - Set deepTools processors manually (e.g., `--numberOfProcessors 4`) if needed.

################################################################################
# USAGE:
#   Normal run: bash modules/pipeline2/02_merge_pseudoreplicates.sh
#   Dry run:    bash modules/pipeline2/02_merge_pseudoreplicates.sh --dry-run
#
# OPTIONS:
#   -d, --dry-run       Simulate execution without changes
#   -h, --help          Show this help message
# USAGE STANDALONE:
#bash 02_merge_Pseudoreplicates.sh \
#  --input-dir /my/bams \
#  --meta /my/meta.tsv \
#  --merge-dir /my/pooled \
#  --pseudo-dir /my/pseudoreps \
#  --stats-dir /my/qc \
#  --log-dir /my/logs \
#  --log-file /my/pooling_summary.tsv
# EXAMPLES:
#   # Normal execution
#   bash modules/pipeline2/02_merge_pseudoreplicates.sh
#
#   # Dry-run mode
#   bash modules/pipeline2/02_merge_pseudoreplicates.sh --dry-run
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
    local message="${2:-}"  # Prevents unbound variable error
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
INPUT_DIR=""
METADATA_TSV=""
MERGE_DIR=""
PSEUDO_DIR=""
STATS_DIR=""
LOG_FILE=""
LOG_DIR=""
HELP_MSG=$'\nUsage: 02_merge_Pseudoreplicates.sh [options]
Options:
  -d, --dry-run              Simulate execution
  --input-dir DIR            Input BAM directory (default: analysis/Renamed_Cleaned)
  --meta FILE                Metadata file (default: metadata/mapping_filtered.tsv)
  --merge-dir DIR            Output directory for pooled BAMs (default: analysis/Pooled_BAMs)
  --pseudo-dir DIR           Output for pseudoreplicates (default: analysis/Pseudoreplicates)
  --stats-dir DIR            Output for statistics and plots (default: analysis/Pool_Pseudo_QC_stats)
  --log-dir DIR              Directory for logs (default: logs/pipeline2)
  --log-file FILE            Summary log file (default: analysis/pooling_log.tsv)
  --correlation-method STR   Correlation method for plotCorrelation [default: pearson]
                             Options: pearson, spearman
  -h, --help                 Show this help message'

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dry-run) DRY_RUN=true; shift ;;
        --input-dir) INPUT_DIR="$2"; shift 2 ;;
        --meta) METADATA_TSV="$2"; shift 2 ;;
        --merge-dir) MERGE_DIR="$2"; shift 2 ;;
        --pseudo-dir) PSEUDO_DIR="$2"; shift 2 ;;
        --stats-dir) STATS_DIR="$2"; shift 2 ;;
        --log-dir) LOG_DIR="$2"; shift 2 ;;
        --log-file) LOG_FILE="$2"; shift 2 ;;
        --correlation-method) CORRELATION_METHOD="$2" ; shift 2 ;;
        -h|--help) echo "$HELP_MSG"; exit 0 ;;
        -*)
            echo "❌ Error: Unknown option $1"
            echo "$HELP_MSG"
            exit 1
            ;;
        *) echo "❌ Error: Positional arguments not allowed"; echo "$HELP_MSG"; exit 1 ;;
    esac
done


# === Configuration ===
METADATA_TSV="${METADATA_TSV:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
INPUT_DIR="${INPUT_DIR:-$PROJECT_ROOT/analysis/Renamed_Cleaned}"
MERGE_DIR="${MERGE_DIR:-$PROJECT_ROOT/analysis/Pooled_BAMs}"
PSEUDO_DIR="${PSEUDO_DIR:-$PROJECT_ROOT/analysis/Pseudoreplicates}"
STATS_DIR="${STATS_DIR:-$PROJECT_ROOT/analysis/Pool_Pseudo_QC_stats}"
CORRELATION_METHOD="${CORRELATION_METHOD:-pearson}"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/logs/pipeline2}"
LOG_FILE="${LOG_FILE:-$PROJECT_ROOT/analysis/pooling_log.tsv}"

MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"

# ===Correlation validation

log "INFO" "📊 Correlation method: $CORRELATION_METHOD"

if [[ ! "$CORRELATION_METHOD" =~ ^(pearson|spearman)$ ]]; then
    echo "❌ Error: Invalid correlation method: '$CORRELATION_METHOD'"
    echo "Allowed values: pearson, spearman"
    exit 1
fi

# === Initialize ===
mkdir -p "$MERGE_DIR" "$PSEUDO_DIR" "$STATS_DIR" "$LOG_DIR" "$(dirname "$LOG_FILE")"
echo -e "Group\tPooled_BAM\tPseudo1\tPseudo2\tInput_BAMs\tTotal_Reads" > "$LOG_FILE"

# === Resource Monitoring Functions ===
record_metrics() {
    local message="$1"
    local timestamp=$(date +%s)
    local cpu_usage=$(ps -p $$ -o %cpu | tail -n 1 | awk '{print $1}')
    local mem_usage=$(ps -p $$ -o %mem | tail -n 1 | awk '{print $1}')
    local disk_usage=$(df -h "$MERGE_DIR" | tail -n 1)
    local bam_count=$(ls -1 "$MERGE_DIR"/*.bam 2>/dev/null | wc -l)
    
    echo "[${timestamp}] ${message}" >> "$PERF_LOG"
    echo "  CPU: ${cpu_usage}% | Memory: ${mem_usage}% | Merged BAMs: ${bam_count}" >> "$PERF_LOG"
    echo "  Disk: ${disk_usage}" >> "$PERF_LOG"
}

start_timer() {
    SECONDS=0
}

get_elapsed_time() {
    local duration=$SECONDS
    echo "$((duration / 60))m $((duration % 60))s"
}

# === Debugging Functions ===
print_header() {
    echo "========================================"
    echo "$1"
    echo "========================================"
}

# === Enhanced Statistics Functions ===
generate_stats() {
    local bam=$1
    local prefix="${STATS_DIR}/$(basename "${bam%.bam}")"
    
    echo "📊 Generating statistics for $bam"
    record_metrics "STATS START: $bam"
    
    if ! samtools flagstat "$bam" > "${prefix}.flagstat"; then
        echo "⚠️ Warning: flagstat failed for $bam"
    fi
    
    if ! samtools idxstats "$bam" > "${prefix}.idxstats"; then
        echo "⚠️ Warning: idxstats failed for $bam"
    fi
    
    local read_count=$(samtools view -c "$bam" 2>/dev/null || echo "NA")
    echo "✅ $read_count reads" | tee "${prefix}.count"
    
    record_metrics "STATS COMPLETE: $bam ($read_count reads)"
}

# === deepTools Visualization ===
# === deepTools Visualization ===
generate_visualization() {
    local group="$1"
    local pooled="$2"

    local pseudo1="${PSEUDO_DIR}/${group}/${group}_pseudo1.bam"
    local pseudo2="${PSEUDO_DIR}/${group}/${group}_pseudo2.bam"

    print_header "GENERATING DEEPTOOLS VISUALIZATION"
    record_metrics "VISUALIZATION START: $group"

    # === 1. Create coverage tracks ===
    if [[ ! -s "${STATS_DIR}/${group}_coverage.bw" ]]; then
        log "INFO" "📈 Creating normalized bigWig coverage track"
        if ! bamCoverage \
            -b "$pooled" \
            -o "${STATS_DIR}/${group}_coverage.bw" \
            --binSize 10 \
            --normalizeUsing CPM \
            --numberOfProcessors 4 \
            --ignoreDuplicates 2>> "$PERF_LOG"; then
            log "WARN" "⚠️ Warning: bamCoverage failed for $group"
        fi
    else
        log "INFO" "📉 Skipping bamCoverage: ${group}_coverage.bw already exists"
    fi

    # === 2. Fingerprint plot ===
    if [[ ! -s "${STATS_DIR}/${group}_fingerprint.png" ]]; then 
        log "INFO" "🖐️ Creating fingerprint plot"
        if ! plotFingerprint \
            -b "$pseudo1" \
               "$pseudo2" \
            --labels "Pseudorep1" "Pseudorep2" \
            -plot "${STATS_DIR}/${group}_fingerprint.png" \
            --outRawCounts "${STATS_DIR}/${group}_fingerprint.tab" \
            --numberOfProcessors 4 2>> "$PERF_LOG"; then
            log "WARN" "⚠️ Warning: plotFingerprint failed for $group"
        fi
    else
        log "INFO" "🖐️ Skipping fingerprint plot: ${group}_fingerprint.png already exists"
    fi

    # === 3. Correlation heatmap ===
    if [[ ! -s "${STATS_DIR}/${group}_correlation_heatmap.png" ]]; then
        log "INFO" "🔥 Creating correlation heatmap"

        if [[ ! -s "${STATS_DIR}/${group}_correlation.npz" ]]; then
            if ! multiBamSummary bins \
                --bamfiles "$pooled" \
                           "$pseudo1" \
                           "$pseudo2" \
                --labels "Pooled" "Pseudorep1" "Pseudorep2" \
                --outFileName "${STATS_DIR}/${group}_correlation.npz" \
                --numberOfProcessors 4 2>> "$PERF_LOG"; then
                log "WARN" "⚠️ Warning: multiBamSummary failed for $group"
                return 1
            fi
        else
            log "INFO" "🧊 Correlation NPZ already exists: ${group}_correlation.npz"
        fi

        if ! plotCorrelation \
            -in "${STATS_DIR}/${group}_correlation.npz" \
            -c "${CORRELATION_METHOD}" \
            -o "${STATS_DIR}/${group}_correlation_heatmap.png" \
            --plotTitle "${group} ${CORRELATION_METHOD^} Correlation" \
            --whatToPlot heatmap \
            --colorMap viridis 2>> "$PERF_LOG"; then
            log "WARN" "⚠️ Warning: plotCorrelation failed for $group"
        fi
    else
        log "INFO" "🔥 Skipping correlation heatmap: ${group}_correlation_heatmap.png already exists"
    fi

    record_metrics "VISUALIZATION COMPLETE: $group"
}

# === Main Processing ===
process_group() {
    local group="$1"
    shift
    local bams=("$@")
    
   print_header "PROCESSING GROUP: $group"
    record_metrics "GROUP START: $group (${#bams[@]} samples)"
    
    log "INFO" "🔵 Found ${#bams[@]} IP BAMs:"
    printf ' - %s\n' "${bams[@]}"
    
    if [ ${#bams[@]} -eq 1 ]; then
        log "WARN" "⚠️  Single IP sample group - skipping pseudoreplicates"
        record_metrics "GROUP SKIPPED: $group (single sample)"
        return
    fi

   
    # === Pooled BAM Generation ===
local pooled="${MERGE_DIR}/pooled_${group}.bam"
pseudo1="${PSEUDO_DIR}/${group}/${group}_pseudo1.bam"
pseudo2="${PSEUDO_DIR}/${group}/${group}_pseudo2.bam"

# Safeguard: skip group if all files exist
if [[ -s "$pooled" && -s "$pooled.bai" && \
      -s "$pseudo1" && -s "$pseudo1.bai" && \
      -s "$pseudo2" && -s "$pseudo2.bai" ]]; then
    log "INFO" "⏩ Skipping group '$group': pooled and pseudoreplicates already exist"

    if [[ ! -f "${STATS_DIR}/${group}_coverage.bw" ]]; then
        generate_visualization "$group" "$pooled"
    fi

    record_metrics "GROUP SKIPPED: $group (already processed)"
    return
fi

log "INFO" "🧬 Merging ${#bams[@]} IP BAMs to $pooled..."

if $DRY_RUN; then
    log "INFO" "🛯 DRY RUN: Would merge to $pooled"
    record_metrics "DRY RUN: Would process group $group"
else
    if [[ -s "$pooled" && -s "$pooled.bai" ]]; then
        log "INFO" "⏩ Skipping merge/index: $pooled and .bai already exist"
    else
        log "INFO" "Merging BAMs..."
        if ! samtools merge -f "$pooled" "${bams[@]}"; then
            log "ERROR" "❌ Failed to merge BAMs for $group"
            record_metrics "GROUP FAILED: Merge error - $group"
            return 1
        fi

        log "INFO" "Indexing merged BAM..."
        if ! samtools index "$pooled"; then
            log "ERROR" "❌ Failed to index $pooled"
            record_metrics "GROUP FAILED: Index error - $group"
            return 1
        fi
    fi

    generate_stats "$pooled"
fi

# === Pseudoreplicate Generation ===
mkdir -p "${PSEUDO_DIR}/${group}"

if $DRY_RUN; then
    log "INFO" "DRY RUN: Would create pseudoreps in ${PSEUDO_DIR}/${group}"
else
    # === Count total reads in pooled BAM ===
    if [[ ! -s "$pooled" ]]; then
        log "ERROR" "❌ BAM file does not exist or is empty: $pooled"
        total_reads="NA"
    else
        log "INFO" "📏 Counting reads in pooled BAM: $pooled"
        if ! total_reads=$(samtools view -c "$pooled" 2>>"$PERF_LOG"); then
            log "ERROR" "❌ samtools view failed to count reads in $pooled"
            echo "[samtools view -c failed] $pooled" >> "$PERF_LOG"
            total_reads="NA"
        fi
    fi

    log "INFO" "CREATING PSEUDOREPLICATES"
    log "INFO" "🔀 Sampling 50% reads for two pseudoreplicates (different seeds)"
    log "INFO" "Total reads in pool: $total_reads"



    if [[ -s "$pseudo1" && -s "$pseudo1.bai" && -s "$pseudo2" && -s "$pseudo2.bai" ]]; then
        log "INFO" "⏩ Skipping pseudoreplicate creation: Both pseudoreps already exist for $group"
    else
        if [[ ! -f "$pseudo1" ]]; then
            log "INFO" "🧬 Creating pseudoreplicate 1..."
            if ! samtools view -s 42.5 -b "$pooled" > "$pseudo1"; then
                log "ERROR" "❌ Failed to create pseudoreplicate 1"
                record_metrics "GROUP FAILED: Pseudorep1 creation - $group"
                return 1
            fi
        fi

        if [[ ! -f "$pseudo2" ]]; then
            log "INFO" "🧬 Creating pseudoreplicate 2..."
            if ! samtools view -s 24.5 -b "$pooled" > "$pseudo2"; then
                log "ERROR" "❌ Failed to create pseudoreplicate 2"
                record_metrics "GROUP FAILED: Pseudorep2 creation - $group"
                return 1
            fi
        fi

        for pseudo in "$pseudo1" "$pseudo2"; do
            log "INFO" "📇 Indexing $pseudo"
            if ! samtools index "$pseudo"; then
                log "WARN" "⚠️ Warning: Failed to index $pseudo"
            fi
            generate_stats "$pseudo"
        done
    fi

    generate_visualization "$group" "$pooled"
fi

# Log everything
echo -e "$group\t$pooled\t$pseudo1\t$pseudo2\t${bams[*]}\t$total_reads" >> "$LOG_FILE"
log "INFO" "✅ Completed processing for $group"
record_metrics "GROUP COMPLETE: $group"

}

# === Main Execution ===
start_timer
{
print_header "🧬 MODULE START: ${SCRIPT_BASE_NAME}"
log "INFO" "🧬 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Metadata: $METADATA_TSV"
log "INFO" "📂 Input BAMs: $INPUT_DIR"
log "INFO" "📦 Merged output: $MERGE_DIR"
log "INFO" "🔁 Pseudoreplicates: $PSEUDO_DIR"
log "INFO" "📊 Stats output: $STATS_DIR"
log "INFO" "🗂️ Log file: $LOG_FILE"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "📦 Script version: ${VERSION} (${SCRIPT_NAME})"
print_header "INITIALIZATION COMPLETE"

record_metrics "PROCESS START"

# === Metadata Processing ===
print_header "PARSING METADATA"
declare -A groups
total_samples=0
processed_samples=0

#while IFS=$'\t' read -r -a cols; do
#    [[ "${cols[0]}" == "Sample_ID" ]] && continue

while IFS=$'\t' read -r -a cols || [[ -n "${cols[*]}" ]]; do
    [[ "${#cols[@]}" -eq 0 ]] && continue
    [[ "${cols[0]}" == "Sample_ID" ]] && continue

    # Only process IP samples (skip Input)
    if [[ ! "${cols[4]}" =~ ^(IP|ChIP|IP_rep[0-9]+)$ ]]; then
        log "WARN" "⚠️  Skipping non-IP sample: ${cols[0]} (Type: ${cols[4]})"
        continue
    fi

    # Find matching BAM file
    bam=$(find "$INPUT_DIR" -name "${cols[0]}_*.bam" -print -quit)
    if [ -z "$bam" ]; then
        log "ERROR" "❌ No BAM found for ${cols[0]} in $INPUT_DIR"
        continue
    fi

    # Group by Condition + Target + Instrument
    key="${cols[2]}"                          # Condition (e.g., WT, BLMneg)
    [ -n "${cols[6]:-}" ] && key+="_${cols[6]}"  # +Target if exists (e.g., G4)
    [ -n "${cols[5]:-}" ] && key+="_${cols[5]}"  # +Instrument if exists
    
    log "INFO" "🔑 Assigning IP sample to group: $key"
    groups["$key"]+="$bam "
    ((total_samples++))
done < "$METADATA_TSV"

# === Show Grouping Results ===
print_header "GROUPING RESULTS"
log "INFO" "🕵️‍♀️ Found ${#groups[@]} IP sample groups from $total_samples samples:"
for group in "${!groups[@]}"; do
    echo "Group '$group' contains:"
    for bam in ${groups["$group"]}; do
        echo " - $bam"
        ((processed_samples++))
    done
done

record_metrics "METADATA PROCESSED: ${#groups[@]} groups from $total_samples samples"

# === Process All Groups ===
print_header "PROCESSING GROUPS"
processed_groups=0
log "INFO" "🧬 Module: ${SCRIPT_BASE_NAME}"
skipped_groups=0
failed_groups=0

for group in "${!groups[@]}"; do
    if process_group "$group" ${groups["$group"]}; then
        ((processed_groups++))
    else
        ((failed_groups++))
    fi
done

echo -e "Group\t#Samples" > "$STATS_DIR/group_summary.tsv"
for group in "${!groups[@]}"; do
    echo -e "$group\t$(echo "${groups[$group]}" | wc -w)" >> "$STATS_DIR/group_summary.tsv"
done

# === Completion ===

print_header  "✅ SCRIPT COMPLETED"
log "INFO" ""
log "INFO" "=============================================="
log "INFO" "🎉 All done! Output locations:"
log "INFO" "🔗 Merged IP BAMs:      $MERGE_DIR"
log "INFO" "🔁 Pseudoreplicates:     $PSEUDO_DIR"
log "INFO" "📊 QC Statistics:        $STATS_DIR"
log "INFO" "🗂️ Run log:              $LOG_FILE"
log "INFO" "📈 Performance metrics:  $PERF_LOG"
log "INFO" "📊 Processing summary:"
log "INFO" "🧮 Total groups:          ${#groups[@]}"
log "INFO" "🔄 Processed:             $processed_groups"
log "INFO" "⏭️ Skipped:               $skipped_groups"
log "INFO" "🧊 Failed:                $failed_groups"
log "INFO" "⏱️ Elapsed time:         $(get_elapsed_time)"
log "INFO" "🕒 End time:             $(date '+%F %T')"
log "INFO" "=============================================="
print_header "END OF MODULE"

record_metrics "PROCESS COMPLETE: $processed_groups groups processed"
} | tee -a "$PERF_LOG" "$MODULE_LOG"
