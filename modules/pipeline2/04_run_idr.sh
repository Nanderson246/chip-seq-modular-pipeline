#!/bin/bash

# Module: run_idr.sh
# Description: Combined script to run IDR for Replicate vs Replicate, Pseudo vs Pseudo, and Pooled/Pseudo vs Replicate
# Author: Merged by ChatGPT based on Nancy Anderson's original scripts

################################################################################
# SOFTWARE REQUIREMENTS (for Docker container build):
#
# Required:
# - bash >= 4.2                 # Required for associative arrays and robust CLI
# - HOMER (>= 4.11)             # Provides findPeaks, makeTagDirectory
# - yq >= 4.0                   # YAML parser used for target styles
# - samtools >= 1.10            # For FRiP scoring (total read counting)
# - bedtools >= 2.29            # For FRiP scoring (peak read overlap)
# - bc                          # Required for floating point operations
# - Rscript                     # Required for format conversion (via custom R script)
#
# Notes:
# - No reference genome required
# - Script assumes presence of R script for peak file conversion
################################################################################


################################################################################
# USAGE 
# bash modules/pipeline2/04_run_idr.sh --source [macs/homer] 
# Dry run: 
# bash modules/pipeline2/04_run_idr.sh --source macs3 --dry-run
# Usage (standalone):
#   bash run_idr.sh \
#     --source macs3 \
#     --mapping /absolute/path/mapping.tsv \
#     --replicate-dir /data/peaks \
#     --pseudo-dir /data/pseudo \
#     --idr-dir /data/idr_results
# NOTE:
#--replicate-dir   Path to replicate peak directory (narrowPeak or broadPeak)

#   This module run IDR analysis for Replicate vs Replicate, Pseudo vs Pseudo, and Pooled/Pseudo vs Replicate peak call done with homer or macs3
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

# === Draw progress bar ===
draw_progress() {
    local name="$1"
    ((current++))

    # Cap current to avoid going over 100%
    local capped_current=$(( current > total_calls ? total_calls : current ))

    local percent=$((capped_current * 100 / total_calls))
    local filled=$((percent / 5))
    ((filled > 20)) && filled=20  # Prevent overflowing the bar

    local feet_bar
    feet_bar="$(printf '👣%.0s' $(seq 1 "$filled"))$(printf '▫️%.0s' $(seq 1 $((20 - filled))))"

    printf "\r📦 Progress: [${BLUE}%-20s${NC}] ${GREEN}%3d%%${NC} (%d/%d) %s" "$feet_bar" "$percent" "$capped_current" "$total_calls" "$name"
}



# === Configuration ===
DRY_RUN=false
SCALED_METADATA="${SCALED_METADATA:-$PROJECT_ROOT/metadata/mapping_scaled.tsv}"
RESULTS_DIR="${RESULTS_DIR:-$PROJECT_ROOT/results}"
BASE_ANALYSIS_DIR="${BASE_ANALYSIS_DIR:-$PROJECT_ROOT/analysis}"
MAPPING="${MAPPING:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
IDR_DIR="${IDR_DIR:-${BASE_ANALYSIS_DIR}/IDR_Results}"
LOG_DIR="${LOG_DIR:-./logs/idr}"
PERF_LOG="$LOG_DIR/IDR_performance.log"
THREADS="${THREADS:-4}"
RANK="${RANK:-p.value}"
SOURCE=""



#=== Performance metrics functions ===
record_metrics() {
local msg="$1"
local timestamp=$(date +%Y-%m-%dT%H:%M:%S)
echo "[$timestamp] $msg" >> "$PERF_LOG"

echo "  Resource Snapshot (from /usr/bin/time):" >> "$PERF_LOG"
grep -E 'Maximum resident set size|User time|System time|Elapsed' /tmp/idr_last_time.log | sed 's/^/    /' >> "$PERF_LOG"

local disk_usage=$(df -h . | awk 'NR==2 {print $4 " free of "$2}')
echo "  Disk: $disk_usage" >> "$PERF_LOG"

}

run_profiled() {
    local label="$1"
    shift
    if $DRY_RUN; then
        echo "[DRY-RUN] 🔍 Would run: $*" | tee -a "$LOG_FILE"
    else
        /usr/bin/time -v "$@" 2> >(tee /tmp/idr_last_time.log >&2)
        record_metrics "$label"
    fi
}

run_and_sort_idr() {
    local base="$1"

    local bed_file="${base}.bed"
    local txt_file="${base}.idr.txt"

    [[ -s "$bed_file" ]] && sort -k1,1 -k2,2n "$bed_file" > "${base}.sorted.bed"
    [[ -s "$txt_file" ]] && sort -k1,1 -k2,2n "$txt_file" > "${base}.sorted.txt"

    log "INFO" "✅ Sorted outputs for $base"
}

# === Safeguard avoidance process repetition ===

idr_exists() {
    local base="$1"
    [[ -s "${base}.idr.txt" && -s "${base}.bed" ]]
}

# === CLI Parser ===
usage() {
    echo ""
    echo "Usage:"
    echo "  $0 --source [homer|macs3] [options]"
    echo ""
    echo "Required:"
    echo "  --source           Peak calling source: 'homer' or 'macs3'"
    echo ""
    echo "Optional:"
    echo "  -m, --mapping FILE         Path to metadata mapping file (default: metadata/mapping_filtered.tsv)"
    echo "  --replicate-dir DIR        Directory with replicate peak files (narrowPeak or broadPeak)"
    echo "  --pseudo-dir DIR           Directory with pseudo-replicate peak files"
    echo "  --idr-dir DIR              Output directory for IDR results (default: ./analysis/IDR_Results)"
    echo "  --log-dir DIR              Log directory (default: ./logs/idr)"
    echo "  --rank FIELD               Ranking field for IDR (default: p.value). Options: p.value, q.value, signal.value"
    echo "  --threads INT              Number of threads to use (default: 4)"
    echo "  --dry-run                  Simulate execution without running commands"
    echo "  -h, --help                 Show this help message and exit"
    echo ""
    echo "Examples:"
    echo "  bash $0 --source macs3 --mapping metadata.tsv"
    echo "  bash $0 --source homer --replicate-dir /path/to/peaks --pseudo-dir /path/to/pseudo --dry-run"
    echo ""
    exit 1
}


# === Parse Args ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m|--mapping) 
           MAPPING="$2"; 
           shift 2
           ;;
        --source) 
          SOURCE="$2"; 
          shift 2
          ;;
        --rank) 
          RANK="$2"; 
          shift 2
          ;;
        --threads) 
          THREADS="$2"; 
          shift 2
          ;;
        --replicate-dir) 
          NARROW_DIR="$2"; 
          shift 2 
          ;;
        --pseudo-dir) 
          PSEUDO_DIR="$2"; 
          shift 2 
          ;;
        --idr-dir) 
          IDR_DIR="$2"; 
          shift 2 
          ;;
        --log-dir) 
          LOG_DIR="$2"; 
          shift 2 ;;
        --dry-run)
          DRY_RUN=true
          shift
          ;;
        -h|--help) 
           usage;;
        *) echo "[ERROR] Unknown option: $1"; usage;;
    esac
done

# === Validate source before proceeding ===
[[ -z "$SOURCE" || ! "$SOURCE" =~ ^(homer|macs3)$ ]] && {
    log "ERROR" "❌ Invalid or missing --source argument (must be 'homer' or 'macs3')"
    exit 1
}
# === HOMER paths ===
if [[ "$SOURCE" == "homer" ]]; then
    HOMER_DIR_BASE="$BASE_ANALYSIS_DIR/PeakCalling_HOMER"
    if [[ -d "$HOMER_DIR_BASE/narrow_PEAKS" ]]; then
        HOMER_DIR="$HOMER_DIR_BASE/narrow_PEAKS"
        log "INFO" "📁 Using HOMER narrow_PEAKS output"
    elif [[ -d "$HOMER_DIR_BASE/broad_PEAKS" ]]; then
        HOMER_DIR="$HOMER_DIR_BASE/broad_PEAKS"
        log "INFO" "📁 Using HOMER broad_PEAKS output"
    else
        log "ERROR" "❌ No HOMER peak folder found (expected narrow_PEAKS or broad_PEAKS under $HOMER_DIR_BASE)"
        exit 1
    fi

    HOMER_PSEUDO_DIR_BASE="$BASE_ANALYSIS_DIR/PeakCalling_HOMER_pool_pseudo"
    if [[ -d "$HOMER_PSEUDO_DIR_BASE/narrow_PEAKS" ]]; then
        HOMER_PSEUDO_DIR="$HOMER_PSEUDO_DIR_BASE/narrow_PEAKS"
        log "INFO" "📁 Using HOMER pseudo narrow_PEAKS output"
    elif [[ -d "$HOMER_PSEUDO_DIR_BASE/broad_PEAKS" ]]; then
        HOMER_PSEUDO_DIR="$HOMER_PSEUDO_DIR_BASE/broad_PEAKS"
        log "INFO" "📁 Using HOMER pseudo broad_PEAKS output"
    else
        log "ERROR" "❌ No HOMER pseudo peak folder found (expected narrow_PEAKS or broad_PEAKS under $HOMER_PSEUDO_DIR_BASE)"
        exit 1
    fi
fi

# === MACS3 paths ===
if [[ "$SOURCE" == "macs3" ]]; then
    MACS3_DIR_BASE="$BASE_ANALYSIS_DIR/PeakCalling_MACS3"
    if [[ -d "$MACS3_DIR_BASE/narrowPeak" ]]; then
        MACS3_DIR="$MACS3_DIR_BASE/narrowPeak"
        log "INFO" "📁 Using MACS3 narrowPeak output"
    elif [[ -d "$MACS3_DIR_BASE/broadPeak" ]]; then
        MACS3_DIR="$MACS3_DIR_BASE/broadPeak"
        log "INFO" "📁 Using MACS3 broadPeak output"
    else
        log "ERROR" "❌ No MACS3 peak folder found (expected narrowPeak or broadPeak under $MACS3_DIR_BASE)"
        exit 1
    fi

    MACS3_PSEUDO_DIR_BASE="$BASE_ANALYSIS_DIR/PeakCalling_MACS3_pool_pseudo"
    if [[ -d "$MACS3_PSEUDO_DIR_BASE/narrowPeak" ]]; then
        MACS3_PSEUDO_DIR="$MACS3_PSEUDO_DIR_BASE/narrowPeak"
        log "INFO" "📁 Using MACS3 pseudo narrowPeak output"
    elif [[ -d "$MACS3_PSEUDO_DIR_BASE/broadPeak" ]]; then
        MACS3_PSEUDO_DIR="$MACS3_PSEUDO_DIR_BASE/broadPeak"
        log "INFO" "📁 Using MACS3 pseudo broadPeak output"
    else
        log "ERROR" "❌ No MACS3 pseudo peak folder found (expected narrowPeak or broadPeak under $MACS3_PSEUDO_DIR_BASE)"
        exit 1
    fi
fi

# === Optional informational banner for user ===
print_header "🦠 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Purpose: run IDR for Replicate vs Replicate, Pseudo vs Pseudo, and Pooled/Pseudo vs Replicate"
log "INFO" "📥 Input: Peak directories "
log "INFO" "📁 Output: $IDR_DIR"
log "INFO" "🕒 Start time: $(date '+%F %T')"
print_header "INITIALIZATION COMPLETE"

# === DRY-RUN announcement ===
if $DRY_RUN; then

detect_peak_folder() {
    local base_dir="$1"
    local narrow_name="$2"
    local broad_name="$3"

    if [[ -d "$base_dir/$narrow_name" ]]; then
        echo "$base_dir/$narrow_name"
    elif [[ -d "$base_dir/$broad_name" ]]; then
        echo "$base_dir/$broad_name"
    else
        log "ERROR" "❌ No peak folder found in $base_dir (expected $narrow_name or $broad_name)"
        exit 1
    fi
}
fi
    log "INFO" "💤 DRY-RUN MODE: Commands will not be executed"
# === Initial Validation ===
[[ -z "$SOURCE" || ! "$SOURCE" =~ ^(homer|macs3)$ ]] && usage
VALID_RANKS=("p.value" "q.value" "signal.value")
[[ ! " ${VALID_RANKS[*]} " =~ " $RANK " ]] && echo "[ERROR] Invalid rank: $RANK" && exit 1

NARROW_DIR="$([ "$SOURCE" == "homer" ] && echo "$HOMER_DIR" || echo "$MACS3_DIR")"
PSEUDO_DIR="$([ "$SOURCE" == "homer" ] && echo "$HOMER_PSEUDO_DIR" || echo "$MACS3_PSEUDO_DIR")"

OUT_DIR="$IDR_DIR/$SOURCE"
LOG_FILE="$LOG_DIR/idr_run_${SOURCE}_$(date +%Y%m%d_%H%M%S).log"
SUMMARY_FILE="$OUT_DIR/idr_summary.tsv"
REPLICATE_OUT="$OUT_DIR/replicate_vs_replicate"
PSEUDO_OUT="$OUT_DIR/pseudo_vs_pseudo"
PSEUDO_REP_OUT="$OUT_DIR/pseudo_vs_replicate"
POOLED_OUT="$OUT_DIR/pooled_vs_replicate"
mkdir -p "$REPLICATE_OUT" "$PSEUDO_OUT" "$PSEUDO_REP_OUT" "$POOLED_OUT"

mkdir -p "$OUT_DIR" "$LOG_DIR"


# === Begin real processing ===
SECONDS=0
run_profiled "PROCESS START" echo "Starting IDR..."

# === Metadata Parsing ===
log "INFO" "🔢 Parsing mapping file for experiment groups..."
HEADER=$(head -n 1 "$MAPPING")
IFS=$'\t' read -ra COLS <<< "$HEADER"

tail -n +2 "$MAPPING" | while IFS=$'\t' read -ra FIELDS; do
    declare -A meta
    for i in "${!COLS[@]}"; do
        meta["${COLS[$i]}"]="${FIELDS[$i]}"
    done

    Sample_ID="${meta[Sample_ID]}"
    Sample_Type="${meta[Sample_Type]}"
    Condition="${meta[Condition]}"
    [[ -z "$Sample_ID" || -z "$Sample_Type" || -z "$Condition" ]] && continue
    [[ "$Sample_Type" != IP* && "$Sample_Type" != "ChIP" ]] && continue

    parts=("${meta[Instrument]}" "${meta[Condition]}" "${meta[Target]}")
    group_key=$(IFS=_; echo "${parts[*]}")

    name_parts=("$Sample_ID")
    for key in Instrument Condition Replicate Sample_Type Target; do
        val="${meta[$key]}"
        val_clean=$(echo "$val" | sed 's/[^a-zA-Z0-9]/_/g')
        [[ -n "$val_clean" ]] && name_parts+=("$val_clean")
    done
    ipname=$(IFS=_; echo "${name_parts[*]}")

    npk=$(find "$NARROW_DIR" -type f \( -name "${ipname}_vs_*_peaks.narrowPeak" -o -name "${ipname}_vs_*_peaks.broadPeak" \) | head -n 1)
   # [[ -n "$npk" ]] && echo "${group_key}|${npk}" >> "$OUT_DIR/group_file_links.txt"
    
    if [[ -f "$npk" ]]; then
        echo "${group_key}|${npk}" >> "$OUT_DIR/group_file_links.txt"
    else
        log "WARN" "⚠️ No peak found for $ipname" | tee -a "$LOG_FILE"
    fi
done

sort -u "$OUT_DIR/group_file_links.txt" -o "$OUT_DIR/group_file_links.txt"
cut -d'|' -f1 "$OUT_DIR/group_file_links.txt" | sort | uniq > "$OUT_DIR/unique_groups.txt"

# === Progress bar tracking ===
#total_calls=$(cut -d'|' -f1 "$OUT_DIR/group_file_links.txt" | sort | uniq | wc -l)

total_calls=$(grep -c '|' "$OUT_DIR/group_file_links.txt")  # count total peak file entries (replicates)
current=0
# === Replicate vs Replicate ===
log "INFO" "🏃‍♂️💨 Running IDR for replicate vs replicate" | tee -a "$LOG_FILE"
while read -r grp; do
    log "INFO" "🎯 IDR for $grp" | tee -a "$LOG_FILE"

    # Extract files linked to group
    files=( $(grep -F "${grp}|" "$OUT_DIR/group_file_links.txt" | cut -d'|' -f2) )
    echo "[DEBUG] Found ${#files[@]} files for group $grp" | tee -a "$LOG_FILE"

    # Warn if not enough replicates
    if [[ ${#files[@]} -lt 2 ]]; then
        log "WARN" "⚠️ Skipping $grp: Only ${#files[@]} replicate(s) found" | tee -a "$LOG_FILE"
        continue
    fi

# Pairwise IDR with canonical ordering
for ((i = 0; i < ${#files[@]} - 1; i++)); do
    for ((j = i + 1; j < ${#files[@]}; j++)); do
        f1="${files[i]}"
        f2="${files[j]}"
        ext="${f1##*.}"
        input_type="$ext"

        # Get tag base names
        tag1=$(basename "$f1" ."$ext")
        tag2=$(basename "$f2" ."$ext")

        # Sort tag names to create canonical pair
        sorted_tags=($(printf "%s\n%s" "$tag1" "$tag2" | sort))
        sorted_tag="${sorted_tags[0]}__vs__${sorted_tags[1]}"

        #base="$OUT_DIR/${sorted_tag}"
        base="$REPLICATE_OUT/${sorted_tag}"
        tmp_log="$base.idr.stderr.log"

        log "INFO" "\n 🔖 Running IDR: ${sorted_tag}" | tee -a "$LOG_FILE"

        log "INFO" "🧪 Checking for existing outputs: $base.idr.txt and $base.bed"
        if idr_exists "$base"; then
            log "INFO" "✅ Skipping ${sorted_tag} – already done"
            continue
        fi

        run_profiled "IDR: $sorted_tag" \
        idr --samples "$f1" "$f2" \
            --input-file-type "$input_type" \
            --output-file "$base.idr.txt" \
            --rank "$RANK" \
            --soft-idr-threshold 0.05 \
            --plot 2> "$base.idr.stats" >"$tmp_log" | tee -a "$LOG_FILE"

        if grep -q "ValueError" "$tmp_log"; then
            log "ERROR" "❌ IDR failed for ${sorted_tag}" | tee -a "$LOG_FILE"
            grep -A2 "ValueError" "$tmp_log" | tee -a "$LOG_FILE"
            continue
        fi

        run_profiled "IDR bed output: ${sorted_tag}" \
        idr --samples "$f1" "$f2" \
            --input-file-type "$input_type" \
            --output-file "$base.bed" \
            --output-file-type bed \
            --rank "$RANK" \
            --soft-idr-threshold 0.05 >> "$LOG_FILE" 2>&1

        echo "[DONE] ${sorted_tag} → $base.bed" | tee -a "$LOG_FILE"
        draw_progress "${sorted_tag}"
        run_and_sort_idr "$base"


        if ! $DRY_RUN; then
            rm -f "$tmp_log"
        fi
    done
done

done < "$OUT_DIR/unique_groups.txt"

# === Pseudo vs Pseudo ===
log "INFO" "🏃‍♂️💨 Running IDR for pseudo1 vs pseudo2" | tee -a "$LOG_FILE"
find "$PSEUDO_DIR" -name "*pseudo1*_peaks.*" | while read -r pseudo1; do
    pseudo2="${pseudo1/pseudo1/pseudo2}"
    [[ ! -f "$pseudo2" ]] && {
        log "WARN" "⚠️ Missing pseudo2 file for $pseudo1" | tee -a "$LOG_FILE"
        continue
    }

    ext="${pseudo1##*.}"
    input_type="$ext"
    tag1=$(basename "$pseudo1" ."$ext")
    tag2=$(basename "$pseudo2" ."$ext")
    #base="$OUT_DIR/${tag1}__vs__${tag2}"
    base="$PSEUDO_OUT/${tag1}__vs__${tag2}"
    tmp_log="$base.idr.stderr.log"

    log "INFO" "\n 🔖 IDR for pseudo_pair: $tag1 vs $tag2" | tee -a "$LOG_FILE"
    
    log "INFO" "🧪 Checking for existing outputs: $base.idr.txt and $base.bed"
    if idr_exists "$base"; then
        log "INFO" "✅ Skipping $tag1 vs $tag2 – already done"
        continue
    fi

    # Run IDR and trap errors
    run_profiled "IDR: $tag1 vs $tag2" \
    idr --samples "$pseudo1" "$pseudo2" \
        --input-file-type "$input_type" \
        --output-file "$base.idr.txt" \
        --rank "$RANK" \
        --soft-idr-threshold 0.05 \
        --plot 2> "$base.idr.stats" >"$tmp_log" | tee a  "$LOG_FILE"

    if grep -q "ValueError" "$tmp_log"; then
        log "ERROR" "❌ Pseudo IDR failed for $tag1 vs $tag2" | tee -a "$LOG_FILE"
        grep -A2 "ValueError" "$tmp_log" | tee -a "$LOG_FILE"
        rm -f "$tmp_log"
        continue
    fi

    # Generate final BED output
    run_profiled "IDR bed output: $tag1 vs $tag2" \
    idr --samples "$pseudo1" "$pseudo2" \
        --input-file-type "$input_type" \
        --output-file "$base.bed" \
        --output-file-type bed \
        --rank "$RANK" \
        --soft-idr-threshold 0.05 >> "$LOG_FILE" 2>&1

    echo "[DONE] Pseudorep IDR: $tag1 vs $tag2 → $base.bed" | tee -a "$LOG_FILE"
    draw_progress "$tag1 vs $tag2"
    run_and_sort_idr "$base"

    if ! $DRY_RUN; then
        rm -f "$tmp_log"
    fi

done

# === Pseudo vs Replicate ===
log "INFO" "🏃‍♂️📊 Running IDR for pseudo vs replicate" | tee -a "$LOG_FILE"

find "$PSEUDO_DIR" -name "*_pseudo1_vs_*_peaks.*" | while read -r pseudo1; do
    ext="${pseudo1##*.}"
    input_type="$ext"
    pseudo_base=$(basename "$pseudo1")
    
    # Extract the input part after "_vs_"
    rep_input_id=$(echo "$pseudo_base" | sed -E 's/.*_vs_([^_]+_.*)_peaks\..*/\1/')

    echo "[DEBUG] Pseudo: $pseudo_base → Rep Input ID: $rep_input_id" | tee -a "$LOG_FILE"

    replicate_peaks=$(find "$NARROW_DIR" -type f -name "*${rep_input_id}_peaks.*")

    echo "[DEBUG] Matching replicates for $rep_input_id: $replicate_peaks" | tee -a "$LOG_FILE"

    for replicate_peak in $replicate_peaks; do
        tag1=$(basename "$pseudo1" ."$ext")
        tag2=$(basename "$replicate_peak" ."$ext")
        base="$PSEUDO_REP_OUT/${tag1}__vs__${tag2}"
        tmp_log="$base.idr.stderr.log"

        log "INFO" "\n 🔖 IDR for pseudo vs replicate: $tag1 vs $tag2" | tee -a "$LOG_FILE"

        if idr_exists "$base"; then
            log "INFO" "✅ Skipping $tag1 vs $tag2 – already done"
            continue
        fi

        run_profiled "IDR: $tag1 vs $tag2" \
        idr --samples "$replicate_peak" "$pseudo1" \
            --input-file-type "$input_type" \
            --output-file "$base.idr.txt" \
            --rank "$RANK" \
            --soft-idr-threshold 0.05 \
            --plot 2> "$base.idr.stats" >"$tmp_log" | tee -a "$LOG_FILE"

        if grep -q "ValueError" "$tmp_log"; then
            log "ERROR" "❌ IDR failed for $tag1 vs $tag2" | tee -a "$LOG_FILE"
            grep -A2 "ValueError" "$tmp_log" | tee -a "$LOG_FILE"
            rm -f "$tmp_log"
            continue
        fi

        run_profiled "IDR bed output: $tag1 vs $tag2" \
        idr --samples "$replicate_peak" "$pseudo1" \
            --input-file-type "$input_type" \
            --output-file "$base.bed" \
            --output-file-type bed \
            --rank "$RANK" \
            --soft-idr-threshold 0.05 >> "$LOG_FILE" 2>&1

        echo "[DONE] Pseudo vs Replicate IDR: $tag1 vs $tag2 → $base.bed" | tee -a "$LOG_FILE"
        draw_progress "$tag1 vs $tag2"
        run_and_sort_idr "$base"

        if ! $DRY_RUN; then
            rm -f "$tmp_log"
        fi
    done
done


# === Pooled/ Replicate ===
log "INFO" "🏃‍♂️💨 Running IDR for pooled/pseudo vs replicate" | tee -a "$LOG_FILE"

find "$PSEUDO_DIR" -name "pooled_*_vs_*_peaks.*" | while read -r pooled; do
    ext="${pooled##*.}"
    input_type="$ext"
     
    #old:greedy and fuzzy
    #rep_id=$(basename "$pooled" | sed -E 's/^pooled_(.*)_vs_(.*)_peaks.*/\2/')
    # SAFE: captures replicate id cleanly
    rep_id=$(basename "$pooled" | sed -E 's/^pooled_.*_vs_([a-zA-Z0-9_]+)_peaks\..*/\1/')
    
    replicate_peaks=$(find "$NARROW_DIR" -name "*${rep_id}_peaks.narrowPeak" -o -name "*${rep_id}_peaks.broadPeak")

    for replicate_peak in $replicate_peaks; do
        tag1=$(basename "$pooled" ."$ext")
        tag2=$(basename "$replicate_peak" ."$ext")
        #base="$OUT_DIR/${tag1}__vs__${tag2}"
        base="$POOLED_OUT/${tag1}__vs__${tag2}"
        tmp_log="$base.idr.stderr.log"

        log "INFO" "\n 🔖 IDR for pooled vs replicate: $tag1 vs $tag2" | tee -a "$LOG_FILE"

        log "INFO" "🧪 Checking for existing outputs: $base.idr.txt and $base.bed"
        if idr_exists "$base"; then
            log "INFO" "✅ Skipping $tag1 vs $tag2 – already done"
            continue
        fi

        run_profiled "IDR: $tag1 vs $tag2" \
        idr --samples "$pooled" "$replicate_peak" \
            --input-file-type "$input_type" \
            --output-file "$base.idr.txt" \
            --rank "$RANK" \
            --soft-idr-threshold 0.05 \
            --plot 2> "$base.idr.stats" >"$tmp_log" | tee -a "$LOG_FILE"

        if grep -q "ValueError" "$tmp_log"; then
            log "ERROR" "❌ IDR failed for $tag1 vs $tag2" | tee -a "$LOG_FILE"
            grep -A2 "ValueError" "$tmp_log" | tee -a "$LOG_FILE"
            rm -f "$tmp_log"
            continue
        fi

        run_profiled "IDR bed output: $tag1 vs $tag2" \
        idr --samples "$pooled" "$replicate_peak" \
            --input-file-type "$input_type" \
            --output-file "$base.bed" \
            --output-file-type bed \
            --rank "$RANK" \
            --soft-idr-threshold 0.05 >> "$LOG_FILE" 2>&1

        echo "[DONE] Pooled vs Replicate IDR: $tag1 vs $tag2 → $base.bed" | tee -a "$LOG_FILE"
        draw_progress "$tag1 vs $tag2"
        run_and_sort_idr "$base"

        if ! $DRY_RUN; then
            rm -f "$tmp_log"
        fi
    done
done


# === DONE ===
log "INFO" "🟢 IDR complete for $SOURCE at $(date)" | tee -a "$LOG_FILE"


log "INFO" "⚙️ Generating IDR summary table..."

SUMMARY_FILE="$OUT_DIR/idr_summary.tsv"
echo -e "Group\tPair\tParams\tTotal\tPass\tPct\tFile\tIDR_Type" > "$SUMMARY_FILE"

find "$OUT_DIR" -type f -name "*.idr.stats" | while IFS= read -r stats; do
    [[ ! -s "$stats" ]] && continue

    filename=$(basename "$stats")
    pair=$(echo "$filename" | sed 's/\.idr\.stats$//' | sed 's/__vs__/ vs /')

    sample1=$(echo "$pair" | cut -d' ' -f1)
    sample2=$(echo "$pair" | cut -d' ' -f3)

    # Default group
    group="Unassigned"

    sample1_base=$(basename "$sample1")
    sample2_base=$(basename "$sample2")

    group=$(grep -F "$sample1_base" "$OUT_DIR/group_file_links.txt" | cut -d'|' -f1 | head -n 1)
    [[ -z "$group" ]] && group=$(grep -F "$sample2_base" "$OUT_DIR/group_file_links.txt" | cut -d'|' -f1 | head -n 1)

    if [[ -z "$group" || "$group" == "Unassigned" ]]; then
        while IFS= read -r g; do
            if [[ "$pair" == *"$g"* ]]; then
                group="$g"
                break
            fi
        done < "$OUT_DIR/unique_groups.txt"
    fi

    # === Determine IDR type from filename ===
    # === Determine IDR type from path ===
case "$stats" in
  */pseudo_vs_pseudo/*)       idr_type="pseudo_vs_pseudo" ;;
  */pooled_vs_replicate/*)    idr_type="pooled_vs_replicate" ;;
  */pseudo_vs_replicate/*)    idr_type="pseudo_vs_replicate" ;;
  */replicate_vs_replicate/*) idr_type="replicate_vs_replicate" ;;
  *)                          idr_type="unknown" ;;
esac


    # Extract stats
    param=$(grep "Final parameter values" "$stats" | grep -o '\[.*\]')
    total=$(grep "Number of reported peaks" "$stats" | awk '{print $NF}' | cut -d'/' -f2)
    pass_line=$(grep "Number of peaks passing IDR cutoff" "$stats")
    pass=$(echo "$pass_line" | awk -F'- ' '{print $2}' | cut -d'/' -f1 | awk '{print $1}')
    pct=$(echo "$pass_line" | grep -o '[0-9.]\+%' | head -1)

    echo -e "${group}\t${pair}\t${param}\t${total}\t${pass}\t${pct}\t$stats\t${idr_type}" >> "$SUMMARY_FILE"
done

#sort -u "$SUMMARY_FILE" -o "$SUMMARY_FILE"


#++++++++++++++++++++++++++++++++


log "INFO" "📊 IDR summary written to $SUMMARY_FILE"

# Copy IDR BED files to annotation-ready folder based on source
if ! $DRY_RUN; then
    if [[ "$SOURCE" == "homer" ]]; then
        mkdir -p "${BASE_ANALYSIS_DIR}/PeakCalling_HOMER/IDR_annotation_ready"

        # Copy individual IDR .bed files
        cp "${IDR_DIR}/homer/"*.bed "${BASE_ANALYSIS_DIR}/PeakCalling_HOMER/IDR_annotation_ready/" 2>/dev/null || true

        # Copy merged replicate IDR if it exists
        if compgen -G "${IDR_DIR}/homer/merged_replicate_idr/*.bed" > /dev/null; then
            cp "${IDR_DIR}/homer/merged_replicate_idr/"*.bed "${BASE_ANALYSIS_DIR}/PeakCalling_HOMER/IDR_annotation_ready/"
        else
            log "WARN" "⚠️ No merged replicate IDR BEDs found for HOMER."
        fi

    elif [[ "$SOURCE" == "macs3" ]]; then
        mkdir -p "${BASE_ANALYSIS_DIR}/PeakCalling_MACS3/IDR_annotation_ready"

        # Copy individual IDR .bed files
        cp "${IDR_DIR}/macs3/"*.bed "${BASE_ANALYSIS_DIR}/PeakCalling_MACS3/IDR_annotation_ready/" 2>/dev/null || true

        # Copy merged replicate IDR if it exists
        if compgen -G "${IDR_DIR}/macs3/merged_replicate_idr/*.bed" > /dev/null; then
            cp "${IDR_DIR}/macs3/merged_replicate_idr/"*.bed "${BASE_ANALYSIS_DIR}/PeakCalling_MACS3/IDR_annotation_ready/"
        else
            log "WARN" "⚠️ No merged replicate IDR BEDs found for MACS3."
        fi
    fi
fi

# Optional plot
if ! $DRY_RUN; then
Rscript modules/pipeline2/04_plot_idr_summary.R "$SUMMARY_FILE" "$SOURCE" "$OUT_DIR/unique_groups.txt" "$OUT_DIR" || echo "[WARN] R plotting script failed"
fi
elapsed="${SECONDS}s"
record_metrics "PROCESS COMPLETE (Elapsed: $elapsed)"


log "INFO" "📑 Sorting final IDR BED and IDR TXT files for downstream compatibility..."

for file in "$OUT_DIR"/*.bed; do
    [[ -s "$file" ]] || continue
    base=$(basename "$file" .bed)
    sort -k1,1 -k2,2n "$file" > "$OUT_DIR/${base}.sorted.bed"
done

log "INFO" "✅ All IDR BEDs sorted."
log "INFO" "📁 Organizing files by condition and merging replicate-vs-replicate IDRs..."

MERGED_DIR="$OUT_DIR/merged_replicate_idr"
mkdir -p "$MERGED_DIR"

while read -r group; do
    group_safe=$(echo "$group" | sed 's/[^a-zA-Z0-9]/_/g')
    group_dir="$OUT_DIR/grouped_sort/${group_safe}"
    mkdir -p "$group_dir"

    group_tokens=($(echo "$group" | tr '_' ' '))  # tokens from group name

    matching_rep_beds=()
    while IFS= read -r file; do
        bed_base=$(basename "$file")
        match_all=true
        for token in "${group_tokens[@]}"; do
            if [[ "$bed_base" != *"$token"* ]]; then
                match_all=false
                break
            fi
        done
        if $match_all; then
            matching_rep_beds+=("$file")
        fi
    done < <(find "$REPLICATE_OUT" -type f -name "*.sorted.bed")
    
     # === DEBUG PRINT: show matching files ===
    log "INFO" "🔍 Matching files for $group: ${#matching_rep_beds[@]}"
    for f in "${matching_rep_beds[@]}"; do
        log "INFO" "   → $(basename "$f")"
    done

    merged="$MERGED_DIR/${group_safe}_merged_idr.bed"

    if [[ ${#matching_rep_beds[@]} -gt 1 ]]; then
        cat "${matching_rep_beds[@]}" | sort -k1,1 -k2,2n | bedtools merge > "$merged"
        
         #  Avoid self-copy
        if [[ "$merged" != "$group_dir/$(basename "$merged")" ]]; then
            cp "$merged" "$group_dir/"
        fi

        log "INFO" "✅ Merged replicate IDR BED for $group → $merged"

    elif [[ ${#matching_rep_beds[@]} -eq 1 ]]; then
        if [[ "${matching_rep_beds[0]}" != "$merged" ]]; then
            cp "${matching_rep_beds[0]}" "$merged"
        fi

        if [[ "$merged" != "$group_dir/$(basename "$merged")" ]]; then
            cp "$merged" "$group_dir/"
        fi

        log "INFO" "ℹ️ Only one replicate-vs-replicate IDR for $group — copying as merged IDR → $merged"

    else
        log "WARN" "⚠️ Not enough replicate IDR .bed files for $group to merge or use"
    fi

    grep -F "$group" "$SUMMARY_FILE" > "$group_dir/idr_summary.tsv"
done < "$OUT_DIR/unique_groups.txt"


log "INFO" "✅ All merged IDR BEDs generated from replicate comparisons."
log "INFO" "✅ All files grouped, summaries written, and replicates merged."
log "INFO" "🎉 Completed organization, categorization, and merging of IDR results."


# === SECTION: Intersect merged IDR peaks with raw replicate peaks ===
print_header "🔬 IDR-RAW INTERSECTION"
log "INFO" "📊 Starting intersection of merged IDR peaks with raw replicate peaks"

# Create output directory for intersection results
INTERSECTION_DIR="$OUT_DIR/IDR_replicate_intersection"
mkdir -p "$INTERSECTION_DIR"

# Initialize status tracking file
STATUS_FILE="$INTERSECTION_DIR/intersection_status.tsv"
echo -e "Group\tStatus" > "$STATUS_FILE"

# Loop over each unique experimental group
while read -r group; do
    # Normalize group name for filenames
    group_safe=$(echo "$group" | sed 's/[^a-zA-Z0-9]/_/g')

    # Path to the merged IDR BED file
    merged_bed="$MERGED_DIR/${group_safe}_merged_idr.bed"

    # If merged file doesn't exist, skip with warning
    if [[ ! -f "$merged_bed" ]]; then
        log "WARN" "⚠️ No merged IDR BED for $group — skipping intersection"
        echo -e "${group}\tSKIPPED_NO_MERGED_IDR" >> "$STATUS_FILE"
        continue
    fi

    # Track if we did any intersections for this group
    any_done=false

    # Find raw peak files that belong to this group using mapping file
    grep -F "$group|" "$OUT_DIR/group_file_links.txt" | cut -d'|' -f2 | while read -r raw_peak; do
        # Extract replicate name (cleaned base name)
        rep_name=$(basename "$raw_peak" | sed 's/_peaks\..*//')

        # Define output file name
        out_file="$INTERSECTION_DIR/${group_safe}_idr_vs_${rep_name}.txt"

        # Perform the intersection
        bedtools intersect -a "$merged_bed" -b "$raw_peak" -wa -wb > "$out_file"

        # Log the completed intersection
        log "INFO" "🔗 Intersected $group (IDR) vs $rep_name (raw) → $out_file"

        any_done=true
    done

    # Record status
    if $any_done; then
        echo -e "${group}\tOK" >> "$STATUS_FILE"
    else
        echo -e "${group}\tNO_RAW_REPLICATES_FOUND" >> "$STATUS_FILE"
    fi

done < "$OUT_DIR/unique_groups.txt"

log "INFO" "✅ Intersection status summary written to → $STATUS_FILE"

# === Enriched intersection summary with explanations ===
EXPLAINED_FILE="$INTERSECTION_DIR/intersection_status_explained.tsv"
echo -e "Group\tStatus\tExplanation" > "$EXPLAINED_FILE"

while read -r group; do
    group_safe=$(echo "$group" | sed 's/[^a-zA-Z0-9]/_/g')
    merged_file="$MERGED_DIR/${group_safe}_merged_idr.bed"

    if [[ ! -f "$merged_file" ]]; then
        status="SKIPPED_NO_MERGED_IDR"
        explanation="Merged IDR BED file not found — likely due to IDR failure or too few replicates."
    else
        raw_count=$(grep -F "$group|" "$OUT_DIR/group_file_links.txt" | cut -d'|' -f2 | wc -l)
        if [[ $raw_count -eq 0 ]]; then
            status="NO_RAW_REPLICATES_FOUND"
            explanation="Group was defined, but no matching raw replicate peak files were found."
        else
            status="OK"
            explanation="Merged IDR and raw replicates found — intersection performed."
        fi
    fi

    echo -e "$group\t$status\t$explanation" >> "$EXPLAINED_FILE"
done < "$OUT_DIR/unique_groups.txt"

log "INFO" "📄 Explained intersection status saved to → $EXPLAINED_FILE"



# === SECTION: Evaluate IDR Peak Quality ===
print_header "📈 IDR QUALITY SUMMARY"
log "INFO" "📊 Evaluating number and quality of IDR peaks per group..."

IDR_EVAL_DIR="$OUT_DIR/IDR_qc_summary"
mkdir -p "$IDR_EVAL_DIR"

EVAL_SUMMARY_FILE="$IDR_EVAL_DIR/idr_quality_metrics.tsv"
echo -e "Group\tMerged_IDR_Peaks\tHigh_Confidence_Peaks(IDR≤0.01)\tComment" > "$EVAL_SUMMARY_FILE"

while read -r group; do
    group_safe=$(echo "$group" | sed 's/[^a-zA-Z0-9]/_/g')

    # Count merged peaks
    merged_file="$MERGED_DIR/${group_safe}_merged_idr.bed"
    if [[ ! -f "$merged_file" ]]; then
        log "WARN" "⚠️ No merged IDR file for $group — skipping evaluation"
        echo -e "${group}\tNA\tNA\tNo merged IDR file" >> "$EVAL_SUMMARY_FILE"
        continue
    fi

    merged_count=$(wc -l < "$merged_file")

    # Look for the best replicate-vs-replicate IDR result (choose first available for now)
    best_idr_file=$(find "$REPLICATE_OUT" -type f -name "*${group_safe}*.idr.txt" | head -n 1)

    if [[ -f "$best_idr_file" ]]; then
        high_conf_count=$(awk '$7 <= 0.01' "$best_idr_file" | wc -l)
        comment="OK"

        if [[ $merged_count -lt 1000 ]]; then
            comment="Too few peaks"
        elif [[ $high_conf_count -lt 500 ]]; then
            comment="Low confidence peaks"
        fi

        echo -e "${group}\t${merged_count}\t${high_conf_count}\t${comment}" >> "$EVAL_SUMMARY_FILE"
        log "INFO" "🧾 $group: ${merged_count} merged, ${high_conf_count} high-confidence → $comment"
    else
        echo -e "${group}\t${merged_count}\tNA\tNo replicate IDR file found" >> "$EVAL_SUMMARY_FILE"
        log "WARN" "⚠️ No IDR .txt file found for $group to evaluate confidence"
    fi
done < "$OUT_DIR/unique_groups.txt"

log "INFO" "✅ IDR quality evaluation summary written to → $EVAL_SUMMARY_FILE"

# === SECTION: Rescue & Self-Consistency Ratios ===
print_header "🔎 IDR REPLICATE AGREEMENT METRICS"
log "INFO" "📊 Calculating Rescue Ratio and Self-Consistency Ratio..."

RATIO_SUMMARY_FILE="$IDR_EVAL_DIR/idr_ratios.tsv"
echo -e "Group\tReplicate_vs_Replicate\tPooled_vs_Pseudo\tPseudo_vs_Pseudo\tRescue_Ratio\tSelf_Consistency_Ratio\tComment" > "$RATIO_SUMMARY_FILE"

while read -r group; do
    group_safe=$(echo "$group" | sed 's/[^a-zA-Z0-9]/_/g')

    # Locate relevant IDR result files
    rep_idr=$(find "$REPLICATE_OUT" -type f -name "*${group_safe}*.idr.txt" | wc -l)
    pool_idr=$(find "$POOLED_OUT" -type f -name "*${group_safe}*.idr.txt" | wc -l)
    pseudo_idr=$(find "$PSEUDO_OUT" -type f -name "*${group_safe}*.idr.txt" | wc -l)

    rep_count=$(find "$REPLICATE_OUT" -type f -name "*${group_safe}*.idr.txt" -exec awk 'END{print NR}' {} + | awk '{sum+=$1} END{print sum}')
    pool_count=$(find "$POOLED_OUT" -type f -name "*${group_safe}*.idr.txt" -exec awk 'END{print NR}' {} + | awk '{sum+=$1} END{print sum}')
    pseudo_count=$(find "$PSEUDO_OUT" -type f -name "*${group_safe}*.idr.txt" -exec awk 'END{print NR}' {} + | awk '{sum+=$1} END{print sum}')

    # Avoid divide-by-zero
    if [[ $rep_count -gt 0 && $pool_count -gt 0 && $pseudo_count -gt 0 ]]; then
        rescue_ratio=$(awk -v a="$pool_count" -v b="$rep_count" 'BEGIN {printf "%.2f", a/b}')
        self_consistency_ratio=$(awk -v a="$rep_count" -v b="$pseudo_count" 'BEGIN {printf "%.2f", a/b}')

        comment="OK"
        if (( $(echo "$rescue_ratio > 2.0" | bc -l) )) || (( $(echo "$self_consistency_ratio > 2.0" | bc -l) )); then
            comment="Low reproducibility"
        fi

        echo -e "$group\t$rep_count\t$pool_count\t$pseudo_count\t$rescue_ratio\t$self_consistency_ratio\t$comment" >> "$RATIO_SUMMARY_FILE"
        log "INFO" "🔍 $group: Rescue=$rescue_ratio, Self-consistency=$self_consistency_ratio → $comment"
    else
        echo -e "$group\t$rep_count\t$pool_count\t$pseudo_count\tNA\tNA\tMissing data" >> "$RATIO_SUMMARY_FILE"
        log "WARN" "⚠️ Incomplete IDR counts for $group — cannot compute ratios"
    fi

done < "$OUT_DIR/unique_groups.txt"

log "INFO" "✅ IDR agreement metrics summary written to → $RATIO_SUMMARY_FILE"
# === SECTION: Jaccard Index Evaluation ===
print_header "🧮 JACCARD INDEX (Overlap % between Replicate Peaks)"
JACCARD_SUMMARY="$IDR_EVAL_DIR/idr_jaccard.tsv"
echo -e "Group\tRep1\tRep2\tJaccard_Index\tIntersection\tUnion\tInterpretation" > "$JACCARD_SUMMARY"

while read -r group; do
    group_safe=$(echo "$group" | sed 's/[^a-zA-Z0-9]/_/g')
    peaks=( $(grep -F "$group|" "$OUT_DIR/group_file_links.txt" | cut -d'|' -f2) )

    if [[ ${#peaks[@]} -lt 2 ]]; then
        continue
    fi

    for ((i = 0; i < ${#peaks[@]} - 1; i++)); do
        for ((j = i + 1; j < ${#peaks[@]}; j++)); do
            a="${peaks[i]}"
            b="${peaks[j]}"
            rep1=$(basename "$a" | sed 's/_peaks\..*//')
            rep2=$(basename "$b" | sed 's/_peaks\..*//')

            result=$(bedtools jaccard -a "$a" -b "$b" | tail -n1)
            index=$(echo "$result" | awk '{print $3}' | tr -d '[:space:]')
            intersect=$(echo "$result" | awk '{print $2}' | tr -d '[:space:]')
            union=$(echo "$result" | awk '{print $1}' | tr -d '[:space:]')

            # Protect against non-numeric or blank index
            if [[ "$index" =~ ^[0-9.]+$ ]]; then
                if (( $(echo "$index >= 0.5" | bc -l) )); then
                    interp="🔥 High similarity"
                elif (( $(echo "$index >= 0.3" | bc -l) )); then
                    interp="⚠️ Moderate similarity"
                else
                    interp="❌ Low similarity"
                fi
            else
                interp="❓ Invalid Jaccard"
                index="NA"
            fi

            echo -e "$group\t$rep1\t$rep2\t$index\t$intersect\t$union\t$interp" >> "$JACCARD_SUMMARY"
            log "INFO" "🔬 Jaccard for $group: $rep1 vs $rep2 = $index ($interp)"
        done
    done
done < "$OUT_DIR/unique_groups.txt"


log "INFO" "✅ Jaccard overlap summary written to → $JACCARD_SUMMARY"

# === SECTION: Final HTML Report (Enhanced) ===
print_header "🧾 GENERATING FINAL HTML QC REPORT"
HTML_REPORT="$IDR_EVAL_DIR/idr_report.html"
: "${TIMESTAMP:=$(date +'%Y-%m-%d %H:%M:%S')}"

cat <<EOF > "$HTML_REPORT"
<html>
<head>
  <title>IDR QC Report</title>
  <p style="font-size: 10px; margin-top: -10px; margin-bottom: 20px;">
  📑 This report summarizes the quality control and reproducibility assessment of ChIP-seq replicates using the Irreproducible Discovery Rate (IDR) framework. It includes analysis of normalized and cleaned BAM files, detection of spike-in usage (when applicable), correlation heatmaps, PCA plots, FRiP scores, and PCR Bottleneck Coefficient (PBC) metrics. BAM files that failed replicate correlation were identified and excluded from downstream analysis. The goal of this report is to identify high-confidence peaks from biologically consistent replicates and provide transparency on sample quality, filtering decisions, and the reproducibility of ChIP enrichment across conditions.
</p>
 
  <style>
body {
  font-family: Georgia, serif;
  font-size: 10px;
  line-height: 1.6;
  padding: 10px;
  max-width: 1100px;
  margin: auto;
  color: #333;
}

    h1, h2 {
      color: #333;
    }
    table {
  width: 100%;
  border-collapse: collapse;
  table-layout: fixed;
  word-wrap: break-word;
  margin-bottom: 40px;
  font-size: 9px;            /* NEW: Smaller text */
  line-height: 1.4;          /* NEW: Tighter line spacing */
  font-family: Georgia, serif;    /* OPTIONAL: Good for numeric data */
}

    th, td {
      border: 1px solid #ccc;
      padding: 8px;
      text-align: left;
    }
    th {
      background-color: #f4f4f4;
    }
    tr:nth-child(even) {
      background-color: #fafafa;
    }
    .timestamp {
      font-size: 0.9em;
      color: #666;
      margin-top: -10px;
      margin-bottom: 30px;
    }
    .logo {
      width: 120px;
      float: right;
    }
    .section {
      clear: both;
    }
    .image-grid {
  display: flex;
  flex-wrap: wrap;
  justify-content: center;
  gap: 20px;
}
.image-grid div {
  flex: 1 0 45%;
  text-align: center;
}
.image-grid img {
  max-width: 90%;
  margin: 10px auto;
}
  </style>
</head>
<body>

<img src="../../../../assets/ninja_8.png" alt="Logo" class="logo" />
<h1>IDR QC Summary</h1>
<p class="timestamp">Generated on: $TIMESTAMP</p>
EOF
# === 1. Spike-In Detection Summary ===
if [[ -f "metadata/mapping_scaled.tsv" ]]; then
  echo '<div class="section">
  <h2>🧬 Spike-In Detection Summary</h2>
  <p> 🧪Spike detection is performed to determine if exogenous spike was added and information was not available </p>
  <table>
    <tr><th>Sample</th><th>Spike Type</th><th>Spike Reads</th><th>Host Reads</th></tr>' >> "$HTML_REPORT"
  tail -n +2 metadata/mapping_scaled.tsv | while IFS=$'\t' read -r sample cell cond rep type inst target spike spikereads hostreads scaling; do
    echo "<tr><td>$sample</td><td>$spike</td><td>$spikereads</td><td>$hostreads</td></tr>" >> "$HTML_REPORT"
  done
  echo '</table></div>' >> "$HTML_REPORT"
fi

# === 2. Replicate QC Plots and 3. Failed Replicates Summary ===
# === Add dynamic section for Replicate QC plots (heatmaps + PCA) ===

echo '<div class="section">
  <h2 style="text-align:center;">📊 Replicate QC Plots</h2>
  <p style="text-align:center;">Correlation heatmaps and PCA plots based on replicate BAM files.</p>
  <div class="image-grid">' >> "$HTML_REPORT"

find "$BASE_ANALYSIS_DIR/Replicate_QC/deeptools" -type f \( -name "*_correlation_heatmap.png" -o -name "*_PCA_plot.png" \) | sort | while read -r plot; do
    plot_name=$(basename "$plot")
    relative_path=$(realpath --relative-to="$(dirname "$HTML_REPORT")" "$plot")
    echo "    <div>
              <img src=\"$relative_path\" alt=\"$plot_name\">
              <div style=\"word-wrap: break-word; max-width: 90%; margin: 5px auto; font-family: monospace; font-size: 0.85em;\">
                $plot_name
              </div>
          </div>" >> "$HTML_REPORT"
           
done 

echo '</div>' >> "$HTML_REPORT"

# === 2b. Fingerprint plots (added) ===
echo '<div class="section">
  <h2 style="text-align:center;">🖋️ Fingerprint Curves (Pre-IDR)</h2>
  <p style="text-align:center;">Signal-to-background profiles for all cleaned BAMs and one representative per group.</p>
  <div class="image-grid">' >> "$HTML_REPORT"

for fp in "$BASE_ANALYSIS_DIR/Replicate_QC/deeptools/"*preIDR_fingerprint.png; do
    [[ -f "$fp" ]] || continue
    fp_name=$(basename "$fp")
    rel_fp=$(realpath --relative-to="$(dirname "$HTML_REPORT")" "$fp")
    echo "  <div><img src=\"$rel_fp\" alt=\"$fp_name\" style=\"max-width:500px;\"><br><code>$fp_name</code></div>" \
      >> "$HTML_REPORT"
done

echo '  </div></div>' >> "$HTML_REPORT"
# === 3. Failed Replicates Summary ===
FAILED_REPLICATE_FILE="${BASE_ANALYSIS_DIR}/Replicate_QC/tmp_failed_ids.txt"
FAILED_BAM_DIR="${BASE_ANALYSIS_DIR}/BAM_replicate_fail"

if [[ -s "$FAILED_REPLICATE_FILE" ]]; then
  echo '<div class="section">
    <h2>🧬 Rejected Samples</h2>
    <p>BAM files that failed correlation with other replicates. Moved to <code>'"$FAILED_BAM_DIR"'</code></p>
    <table>
      <tr><th>Original ID</th><th>Renamed BAM(s)</th></tr>' >> "$HTML_REPORT"

  while IFS= read -r failed_id; do
    [[ -z "$failed_id" ]] && continue
    # Find all BAM/BAI files matching the failed ID
    matched_files=$(find "$FAILED_BAM_DIR" -type f \( -name "*${failed_id}*.bam" -o -name "*${failed_id}*.bai" \) | sort)

    # Format matched file names
    file_list=""
    while IFS= read -r file; do
      rel_path=$(realpath --relative-to="$(dirname "$HTML_REPORT")" "$file")
      file_list+="<code>$rel_path</code><br>"
    done <<< "$matched_files"

    echo "<tr><td><code>$failed_id</code></td><td>$file_list</td></tr>" >> "$HTML_REPORT"
  done < "$FAILED_REPLICATE_FILE"

  echo '</table>
</div>' >> "$HTML_REPORT"
fi

# === 4. PBC Summary ===
echo '<div class="section">
  <h2>🧬 PCR Bottleneck Coefficient (PBC)</h2>' >> "$HTML_REPORT"

PBC_REPORT="${BASE_ANALYSIS_DIR}/Replicate_QC/pbc_metrics.tsv"

if [[ -f "$PBC_REPORT" ]]; then
    echo "<table><tr><th>Sample</th><th>PBC</th><th>N1</th><th>Nd</th><th>Cutoff</th><th>Interpretation</th></tr>" >> "$HTML_REPORT"
    tail -n +2 "$PBC_REPORT" | while IFS=$'\t' read -r sample pbc n1 nd cutoff interp; do
        if (( $(echo "$pbc >= 0.9" | bc -l) )); then
            color="style=\"background-color:#d4edda\""
        elif (( $(echo "$pbc >= 0.7" | bc -l) )); then
            color="style=\"background-color:#fff3cd\""
        elif (( $(echo "$pbc >= 0.5" | bc -l) )); then
            color="style=\"background-color:#ffeeba\""
        else
            color="style=\"background-color:#f8d7da\""
        fi

        echo "<tr>
                <td>$sample</td>
                <td $color><strong>$pbc</strong></td>
                <td>$n1</td>
                <td>$nd</td>
                <td>$cutoff</td>
                <td>$interp</td>
              </tr>" >> "$HTML_REPORT"
    done

    echo "</table>" >> "$HTML_REPORT"
    echo "<p><strong>Legend:</strong> 
      <span style='background-color:#d4edda'>&nbsp;High (≥ 0.9)&nbsp;</span> 
      <span style='background-color:#fff3cd'>&nbsp;Moderate (0.7–0.89)&nbsp;</span> 
      <span style='background-color:#ffeeba'>&nbsp;Low (0.5–0.69)&nbsp;</span> 
      <span style='background-color:#f8d7da'>&nbsp;Very Low (&lt; 0.5)&nbsp;</span></p>" >> "$HTML_REPORT"
else
    echo "<p><em>⚠️ PBC summary not found at $PBC_REPORT</em></p>" >> "$HTML_REPORT"
fi

echo '</div>' >> "$HTML_REPORT"

# === 5. IDR Evaluation Plots ===
echo '<div class="section">
  <h2>📈 IDR Evaluation Plots</h2>
  <p>These plots summarize IDR pass rates and parameter distributions for <code>'"$SOURCE"'</code> peaks.</p>
  <div class="image-grid">' >> "$HTML_REPORT"

IDR_PLOT_DIR="$IDR_DIR/$SOURCE"
if [[ -d "$IDR_PLOT_DIR" ]]; then
    find "$IDR_PLOT_DIR" -type f -name "*.png" | sort | while read -r img; do
        img_name=$(basename "$img")
        formatted_name=$(echo "$img_name" | sed 's/__vs__/__vs__<br>/g')
        relative_path=$(realpath --relative-to="$(dirname "$HTML_REPORT")" "$img")
        echo "    <div><img src=\"$relative_path\" alt=\"$img_name\" style=\"max-width: 500px;\"><br><code>$formatted_name</code></div>" >> "$HTML_REPORT"
    done
else
    echo "<p><em>No IDR plots found for source: $SOURCE</em></p>" >> "$HTML_REPORT"
fi

echo '  </div>
</div>' >> "$HTML_REPORT"


# === 6. FRiP Summary ===
# === SECTION: Add FRiP summary based on peak caller ===
echo '<div class="section">
  <h2>📊 FRiP Scores from Initial Peak Calling</h2>' >> "$HTML_REPORT"

if [[ "$SOURCE" == "homer" ]]; then
    for frip_path in "${BASE_ANALYSIS_DIR}/PeakCalling_HOMER/frip_summary.tsv" "${BASE_ANALYSIS_DIR}/PeakCalling_HOMER_pool_pseudo/frip_summary.tsv"; do
        if [[ -f "$frip_path" ]]; then
            caller_name=$(basename "$(dirname "$frip_path")")
            echo "<h3>📦 $caller_name</h3>" >> "$HTML_REPORT"
            echo "<table><tr><th>Sample</th><th>Tag Count</th><th>Peak Count</th><th>FRiP</th></tr>" >> "$HTML_REPORT"
            tail -n +2 "$frip_path" | while IFS=$'\t' read -r sample tags peaks frip; do
                echo "<tr><td>$sample</td><td>$tags</td><td>$peaks</td><td>$frip</td></tr>" >> "$HTML_REPORT"
            done
            echo "</table>" >> "$HTML_REPORT"
        fi
    done
elif [[ "$SOURCE" == "macs3" ]]; then
    for frip_path in "${BASE_ANALYSIS_DIR}/PeakCalling_MACS3/frip_summary.tsv" "${BASE_ANALYSIS_DIR}/PeakCalling_MACS3_pool_pseudo/frip_summary.tsv"; do
        if [[ -f "$frip_path" ]]; then
            caller_name=$(basename "$(dirname "$frip_path")")
            echo "<h3>📦 $caller_name</h3>" >> "$HTML_REPORT"
            echo "<table><tr><th>Sample</th><th>Tag Count</th><th>Peak Count</th><th>FRiP</th></tr>" >> "$HTML_REPORT"
            tail -n +2 "$frip_path" | while IFS=$'\t' read -r sample tags peaks frip; do
                echo "<tr><td>$sample</td><td>$tags</td><td>$peaks</td><td>$frip</td></tr>" >> "$HTML_REPORT"
            done
            echo "</table>" >> "$HTML_REPORT"
        fi
    done
else
    echo "<p><em>⚠️ Unknown peak caller source: $SOURCE — FRiP section skipped</em></p>" >> "$HTML_REPORT"
fi

echo '</div>' >> "$HTML_REPORT"


cat <<EOF >> "$HTML_REPORT"
<div class="section">
  <h2>🔹 High-confidence Peak Metrics</h2>
  <table>
    <tr><th>Group</th><th>Merged_IDR_Peaks</th><th>High_Confidence_Peaks(IDR≤0.01)</th><th>Comment</th></tr>
EOF

# High-confidence rows
tail -n +2 "$EVAL_SUMMARY_FILE" | while IFS=$'\t' read -r group merged high_conf comment; do
  echo "<tr><td>$group</td><td>$merged</td><td>$high_conf</td><td>$comment</td></tr>" >> "$HTML_REPORT"
done

cat <<EOF >> "$HTML_REPORT"
  </table>
</div>

<div class="section">
  <h2>🔹 Replicate Agreement (Rescue & Self-Consistency)</h2>
  <table>
    <tr><th>Group</th><th>Replicate_vs_Replicate</th><th>Pooled_vs_Pseudo</th><th>Pseudo_vs_Pseudo</th><th>Rescue_Ratio</th><th>Self_Consistency_Ratio</th><th>Comment</th></tr>
EOF

tail -n +2 "$RATIO_SUMMARY_FILE" | while IFS=$'\t' read -r group rep rep_pseudo pseudo_pseudo rescue self comment; do
  echo "<tr><td>$group</td><td>$rep</td><td>$rep_pseudo</td><td>$pseudo_pseudo</td><td>$rescue</td><td>$self</td><td>$comment</td></tr>" >> "$HTML_REPORT"
done

cat <<EOF >> "$HTML_REPORT"
  </table>
</div>

<div class="section">
  <h2>🔹 Jaccard Overlap (Replicate Peak Similarity)</h2>
  <table>
    <tr>
      <th>Group</th>
      <th>Replicate 1</th>
      <th>Replicate 2</th>
      <th>Jaccard Index</th>
      <th>Intersection</th>
      <th>Union</th>
      <th>Interpretation</th>
    </tr>
EOF

tail -n +2 "$JACCARD_SUMMARY" | while IFS=$'\t' read -r group rep1 rep2 jaccard_index intersection union interpretation; do
  echo "<tr>
      <td>$group</td>
      <td style=\"word-break: break-word;\">$rep1</td>
      <td style=\"word-break: break-word;\">$rep2</td>
      <td>$jaccard_index</td>
      <td>$intersection</td>
      <td>$union</td>
      <td>$interpretation</td>
  </tr>" >> "$HTML_REPORT"
done

cat <<EOF >> "$HTML_REPORT"
  </table>
  <p><strong>Legend:</strong></p>
  <ul>
    <li>🔥 High similarity (≥ 0.5)</li>
    <li>⚠️ Moderate similarity (0.3–0.49)</li>
    <li>❌ Low similarity (&lt; 0.3)</li>
  </ul>
  <p><strong>Value range:</strong><br>
    • 0 → No overlap at all<br>
    • 1 → Perfect overlap (identical peak sets)
  </p>
  <p> </p>
  <p> </p>
  <p style="margin-top: 20px;"></p>
</div>

</body>
</html>
EOF

log "INFO" "✅ Final HTML report generated → $HTML_REPORT"

# === Save a PDF copy of the HTML report ===
HTML_PDF="${HTML_REPORT%.html}.pdf"

if [[ ! -f "$HTML_REPORT" ]]; then
    log "ERROR" "❌ HTML report not found: $HTML_REPORT"
elif command -v chromium >/dev/null 2>&1; then
    ABS_HTML="$(realpath "$HTML_REPORT")"
chromium --headless --disable-gpu --no-sandbox \
         --virtual-time-budget=10000 \
         --print-to-pdf="$HTML_PDF" "file://$ABS_HTML" \
         2>/dev/null

    log "INFO" "✅ PDF generated: $HTML_PDF"
elif command -v google-chrome >/dev/null 2>&1; then
     ABS_HTML="$(realpath "$HTML_REPORT")"
    google-chrome --headless --disable-gpu --no-sandbox \
                   --virtual-time-budget=10000 \
                  --print-to-pdf="$HTML_PDF" "file://$ABS_HTML" \
                  2>/dev/null
    log "INFO" "✅ PDF generated: $HTML_PDF"
else
    log "WARN" "⚠️ Neither Chromium nor Google Chrome is installed. Skipping PDF export."
fi

# === Final message ===
print_header "📦 MODULE COMPLETE"
log "INFO" "🧾 IDR analysis complete: $(date)"
print_header "✅ [MODULE: run_idr.sh] Finished at $(date '+%F_%H-%M-%S')"
