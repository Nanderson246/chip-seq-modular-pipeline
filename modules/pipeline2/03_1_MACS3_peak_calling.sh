#!/usr/bin/env bash

# Module: MACS3_peak_calling with MACS3 (flexible metadata support)
# Author: Nancy Anderson
# Description: Util script to call peaks with MACS3

################################################################################
# SOFTWARE REQUIREMENTS:
# - macs3==3.0.3        # For peak calling
# - samtools >= 1.10    # For read counting in BAM
# - bedtools >= 2.29    # For computing FRiP scores
# - yq >= 4.0           # For parsing YAML (must be installed as `yq` from  
# https://github.com/mikefarah/yq)
# Optional:
# - bash >= 4.2         # For associative arrays and better CLI parsing
# - bc                  # For floating point arithmetic in FRiP scoring
################################################################################

################################################################################
# USAGE (Standalone):
#    bash modules/pipeline2/03_1_MACS3_peak_calling.sh --mapping metadata/mapping_filtered.tsv
#  or
#    bash modules/pipeline2/03_1_MACS3_peak_calling.sh --mapping metadata/mapping_filtered.tsv
# NOTE:
#   This module renamed BAM files from the results/Filtered/Cleaned directory to adjust name to experiment metadata and conditions.
#   Once the instrumentation, sample type, target, and Input is assigned, peaks can be called using MACS3 or HOMER.
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

# === Default paths (can be overridden by CLI or env) ===
BAM_DIR="${BAM_DIR:-$PROJECT_ROOT/analysis/Renamed_Cleaned}"
PEAK_DIR="${PEAK_DIR:-$PROJECT_ROOT/analysis/PeakCalling_MACS3}"
MAPPING="${MAPPING:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
TARGET_YAML="${TARGET_YAML:-$PROJECT_ROOT/templates/target_peak_rules.yaml}"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/logs/pipeline2}"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_performance_${TIMESTAMP}.log"
REF_PREFIX="${REF_PREFIX:-hg38}"
mkdir -p "$PEAK_DIR"
: > "$MODULE_LOG"
: > "$PERF_LOG"

# === Performance metrics ===
record_metrics() {
    local msg="$1"
    local timestamp
    timestamp=$(date +%Y-%m-%dT%H:%M:%S)
    local cpu
    cpu=$(ps -p $$ -o %cpu= | xargs)
    local mem
    mem=$(ps -p $$ -o %mem= | xargs)
    local disk
    disk=$(df -h "$PEAK_DIR" | tail -1)

    echo "[$timestamp] $msg" >> "$PERF_LOG"
    echo "  CPU: $cpu% | MEM: $mem%" >> "$PERF_LOG"
    echo "  Disk: $disk" >> "$PERF_LOG"
}

# === Draw progress bar ===
draw_progress() {
    local name="$1"
    ((current++))
    local percent=$((current * 100 / total_calls))
    local filled=$((percent / 5))
    local feet_bar
    feet_bar="$(printf '👣%.0s' $(seq 1 "$filled"))$(printf '▫️%.0s' $(seq 1 $((20 - filled))))"
    printf "\r📦 Progress: [${BLUE}%-20s] ${GREEN}%3d%% (%d/%d) %s" "$feet_bar" "$percent" "$current" "$total_calls" "$name"
}

# === Optional informational banner for user ===

print_header "🦠 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Purpose: Peak calling with MACS3"
log "INFO" "📥 Input: $BAM_DIR"
log "INFO" "📁 Output:  $PEAK_DIR"
log "INFO" "🕒 Start time: $(date '+%F %T')"
print_header "INITIALIZATION COMPLETE"

# === CLI Parser ===
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Optional arguments:"
    echo "  -m, --mapping FILE        Metadata mapping file (default: metadata/mapping.tsv)"
    echo "  -b, --bam-dir DIR         Directory containing input BAM files"
    echo "  -o, --out-dir DIR         Output directory for MACS3 results"
    echo "  -y, --yaml FILE           YAML file with broad peak target list"
    echo "  --force-broad             Force MACS3 to use --broad mode"
    echo "  -d, --dry-run             Simulate execution without modifying files"
    echo "  -h, --help                Show this help message"
    exit 1
}

FORCE_BROAD=false
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d | --dry-run)
            DRY_RUN=true
            shift
            ;;
        -m | --mapping)
            MAPPING="$2"
            shift 2
            ;;
        -b | --bam-dir)
            BAM_DIR="$2"
            shift 2
            ;;
        -o | --out-dir)
            PEAK_DIR="$2"
            shift 2
            ;;
        -y | --yaml)
            TARGET_YAML="$2"
            shift 2
            ;;
        --force-broad | --broad)
            FORCE_BROAD=true
            shift
            ;;
        -h | --help) usage ;;
        *)
            echo "[ERROR] 🤷‍♀️ Unknown option: $1"
            usage
            ;;
    esac
done

# === Progress bar tracking ===
total_calls=$(grep -v "^#" "$MAPPING" | grep -E "ChIP|^IP" | wc -l)
current=0

if (( total_calls == 0 )); then
    log "WARN" "⚠️ No ChIP/IP samples found — exiting."
    exit 0
fi


# === Initial Validation ===
if [[ ! -f "$MAPPING" ]]; then
    log "ERROR" "❌ Mapping file not found: $MAPPING"
    exit 1
fi

HEADER=$(head -n 1 "$MAPPING")
IFS=$'\t' read -ra REQUIRED_COLS <<< "$HEADER"

for field in Sample_ID Sample_Type Condition Target; do
    log "INFO" "🗃️ Parsed header columns: ${REQUIRED_COLS[*]}"
    if ! [[ " ${REQUIRED_COLS[@]} " =~ " ${field} " ]]; then
        log "ERROR" "❌ Missing required field in mapping file: $field"
        exit 1
    fi
done

mkdir -p "$PEAK_DIR"

# === Genome size determination ===
case "$REF_PREFIX" in
  hg38|GRCh38) GENOME_SIZE="hs" ;;
  mm10|GRCm38) GENOME_SIZE="mm" ;;
  dm6)         GENOME_SIZE="dm" ;;
  ce11|WBcel235) GENOME_SIZE="ce" ;;
  *) 
    log "ERROR" "🤷‍♀️ Unknown genome 🧬 prefix: $REF_PREFIX"
    exit 1
    ;;
esac

# === MACS3 availability check ===
if ! command -v macs3 &> /dev/null; then
    log "ERROR" "❌ macs3 tools 🔧 not found in PATH."
    exit 1
fi
REQUIRED_VERSION="3.0.3"
current_version=$(macs3 --version | awk '{print $2}')
if [[ "$current_version" != "$REQUIRED_VERSION" ]]; then
    log "WARN" "MACS3 version mismatch: expected $REQUIRED_VERSION, got $current_version"
fi

# === Samtools and bedtools availability check ===
if ! command -v samtools &>/dev/null || ! command -v bedtools &>/dev/null; then
    log "ERROR" "❌ samtools or bedtools not found — required for FRiP calculation"
    exit 1
fi

# === Setup logging ===
log "INFO" "⛰️ Peak calling started: $(date)" | tee -a "$MODULE_LOG"

# === Setup Dry run ===
if $DRY_RUN; then
    log "INFO" " ✅DRY-RUN mode enabled — no files will be modified"
    echo "[INFO] DRY-RUN mode enabled" >> "$MODULE_LOG"
fi

# === Progress bar tracking ===
total_calls=$(grep -v "^#" "$MAPPING" | grep -E "ChIP|^IP" | wc -l)
current=0

if (( total_calls == 0 )); then
    log "WARN" "⚠️ No ChIP/IP samples found — exiting."
    exit 0
fi

# === Begin real processing ===
SECONDS=0
record_metrics "PROCESS START"

# === Load broad peak targets from YAML ===
if [[ -f "$TARGET_YAML" ]]; then
    BROAD_TARGETS=$(yq e '.broad_peak_targets[]' "$TARGET_YAML")
else
    log "ERROR" "❌ target.yaml not found 👻 at $TARGET_YAML"
    exit 1
fi

# === Parse metadata and call peaks ===
{
  read -r _  # Skip header
  while IFS=$'\t' read -r line; do
    [[ -z "$line" ]] && continue
    
    IFS=$'\t' read -ra FIELDS <<< "$line"
    declare -A meta
    for i in "${!REQUIRED_COLS[@]}"; do
        meta["${REQUIRED_COLS[$i]}"]="${FIELDS[$i]:-}"
    done

    Sample_ID="${meta[Sample_ID]}"
    Sample_Type="${meta[Sample_Type]}"
    Condition="${meta[Condition]}"
    Target="${meta[Target]}"

    [[ -z "$Sample_ID" || -z "$Sample_Type" || -z "$Condition" || -z "$Target" ]] && continue
    [[ "$Sample_Type" != IP* && "$Sample_Type" != "ChIP" ]] && continue

 # Build base name safely
    parts=("$Sample_ID")
    for key in Instrument Condition Replicate Sample_Type Target; do
        val="${meta[$key]:-}"
        val_clean=$(echo "$val" | sed 's/[^a-zA-Z0-9]/_/g')
        [[ -n "$val_clean" ]] && parts+=("$val_clean")
    done
    base=$(IFS=_; echo "${parts[*]}")
    bam_path="${BAM_DIR}/${base}.bam"

# 🔧 Assign target for clarity/logging
   Target="${meta[Target]}"

  # Decide broad or narrow peak mode
    USE_BROAD=""
    REASON=""

    if echo "$BROAD_TARGETS" | grep -qw "${meta[Target]}"; then
        USE_BROAD="--broad"
        REASON="Target=${meta[Target]} matched broad list"
    else
        REASON="Target=${meta[Target]} did not match broad list"
    fi

    $FORCE_BROAD && {
        USE_BROAD="--broad"
        REASON="User override via --force-broad"
    }

    # === Find matching control sample ===
    input_line=""
    {
      read -r _  # Skip header again
      while IFS=$'\t' read -r control_line; do
        [[ -z "$control_line" ]] && continue
        
        IFS=$'\t' read -ra control_fields <<< "$control_line"
        declare -A control_meta
        for i in "${!REQUIRED_COLS[@]}"; do
            control_meta["${REQUIRED_COLS[$i]}"]="${control_fields[$i]:-}"
        done
        if [[ "${control_meta[Sample_Type]}" =~ ^(Input|IgG|Mock)$ ]] && \
           [[ "${control_meta[Condition]}" == "$Condition" ]] && \
           [[ "${control_meta[Instrument]:-}" == "${meta[Instrument]:-}" ]] && \
           [[ "${control_meta[Target]}" == "$Target" ]]; then
            input_line="$control_line"
            break
        fi
      done
    } < "$MAPPING"

    if [[ -z "$input_line" ]]; then
        log "WARN" "⚠️ No matching Input/IgG/Mock for $base — Calling peaks without control!"
        # Define output directory and name
        name="${base}_noInput"
        macs_out_dir="$PEAK_DIR/$name"
        
 # Skip if ZERO_PEAKS_MARKER was previously set
        if [[ -f "$macs_out_dir/ZERO_PEAKS_MARKER" ]]; then
            log "INFO" "⚠️ Skipping $name: previously produced zero peaks"
            continue
        fi       
        
  # === SAFEGUARD: skip if all expected outputs are present ===
        expected_xls_file="$macs_out_dir/${name}_peaks.xls"
        expected_peak_file="$macs_out_dir/${name}_peaks.narrowPeak"
        expected_summit_file="$macs_out_dir/${name}_summits.bed"

        if [[ "$USE_BROAD" == "--broad" ]]; then
            expected_peak_file="$macs_out_dir/${name}_peaks.broadPeak"
            expected_summit_file=""
        fi

# === Evaluate final file check ===
        skip_call=false
        [[ -s "$expected_peak_file" && -s "$expected_xls_file" ]] && skip_call=true
        if [[ -n "$expected_summit_file" && ! -s "$expected_summit_file" ]]; then
            skip_call=false
        fi

        if $skip_call; then
            log "INFO" "✅ Skipping $name: All expected peak outputs exist"
            continue
        fi
       
        if $DRY_RUN; then
            log "INFO" "⚡ [DRY-RUN] Would call MACS3: $name (no input), $USE_BROAD"
        else
            mkdir -p "$macs_out_dir"
            macs3 callpeak \
              -t "$bam_path" \
              -f BAM \
              -g "$GENOME_SIZE" \
              -n "$name" \
              --outdir "$macs_out_dir" \
              --keep-dup all \
              --nomodel \
              --extsize 200 \
              --qvalue 0.01 \
              $USE_BROAD 2>&1 | tee -a "$MODULE_LOG"
              
            # === Check for zero peaks ===
             if [[ -f "$expected_peak_file" ]]; then
                 peak_count=$(wc -l < "$expected_peak_file")
                 if [[ "$peak_count" -eq 0 ]]; then
                    log "WARN" "⚠️ MACS3 call for $name produced ZERO peaks"
                    touch "$macs_out_dir/ZERO_PEAKS_MARKER"
                 fi
             else
                 log "ERROR" "❌ Expected peak file missing for $name"
                 touch "$macs_out_dir/MISSING_PEAK_FILE_MARKER"
             fi

        log "INFO" "📦 Finished MACS3 call (no input) for $name"
        draw_progress "$name"

        fi
                        
    else
        IFS=$'\t' read -ra input_fields <<< "$input_line"
        declare -A input_meta
        for i in "${!REQUIRED_COLS[@]}"; do
            input_meta["${REQUIRED_COLS[$i]}"]="${input_fields[$i]:-}"
        done

        Input_ID="${input_meta[Sample_ID]}"
        
        parts_input=("$Input_ID")
        for key in Instrument Condition Replicate Sample_Type Target; do
            val="${input_meta[$key]:-}"
            val_clean=$(echo "$val" | sed 's/[^a-zA-Z0-9]/_/g')
            [[ -n "$val_clean" ]] && parts_input+=("$val_clean")
        done
        Input_base=$(
            IFS=_; 
            echo "${parts_input[*]}"
            )
        input_path="${BAM_DIR}/${Input_base}.bam"

        if [[ ! -f "$bam_path" ]]; then
            log "ERROR" "❌ Missing treatment BAM: $bam_path"
            continue
        fi
        if [[ ! -f "$input_path" ]]; then
            log "ERROR" "❌ Missing control 🕳️ BAM: $input_path"
            continue
        fi

        # Define output directory and name
        name="${base}_vs_${Input_base}"
        macs_out_dir="$PEAK_DIR/$name"
        
        # Skip if ZERO_PEAKS_MARKER 
        if [[ -f "$macs_out_dir/ZERO_PEAKS_MARKER" ]]; then
            log "INFO" "⚠️ Skipping $name: previously produced zero peaks"
            continue
        fi
        
        
        # === SAFEGUARD: skip if all expected outputs are present ===
        expected_xls_file="$macs_out_dir/${name}_peaks.xls"
        expected_peak_file="$macs_out_dir/${name}_peaks.narrowPeak"
        expected_summit_file="$macs_out_dir/${name}_summits.bed"

        if [[ "$USE_BROAD" == "--broad" ]]; then
            expected_peak_file="$macs_out_dir/${name}_peaks.broadPeak"
            expected_summit_file=""
        fi

        # Evaluate final file check
        skip_call=false
        [[ -s "$expected_peak_file" && -s "$expected_xls_file" ]] && skip_call=true
        if [[ -n "$expected_summit_file" && ! -s "$expected_summit_file" ]]; then
    skip_call=false
        fi

        if $skip_call; then
            log "INFO" "✅ Skipping $name: All expected peak outputs exist"
            continue
        fi

        if $DRY_RUN; then
            log "INFO" "⚡ [DRY-RUN] Would call MACS3: $name (with input), $USE_BROAD"
        else
            mkdir -p "$macs_out_dir"
            log "INFO" " ⛰️ Calling MACS3 for $name"
            log "INFO" "$base: Applying MACS3 with ${USE_BROAD:-narrow peak mode} — $REASON"
            macs3 callpeak \
              -t "$bam_path" \
              -c "$input_path" \
              -f BAM \
              -g "$GENOME_SIZE" \
              -n "$name" \
              --outdir "$macs_out_dir" \
              --keep-dup all \
              --nomodel \
              --extsize 200 \
              --qvalue 0.01 \
              $USE_BROAD 2>&1 | tee -a "$MODULE_LOG"
              
         # === Check for zero peaks ===
            if [[ -f "$expected_peak_file" ]]; then
                peak_count=$(wc -l < "$expected_peak_file")
                if [[ "$peak_count" -eq 0 ]]; then
                    log "WARN" "⚠️ MACS3 call for $name produced ZERO peaks"
                    touch "$macs_out_dir/ZERO_PEAKS_MARKER"
                fi
            else
                log "ERROR" "❌ Expected peak file missing for $name"
                touch "$macs_out_dir/MISSING_PEAK_FILE_MARKER"
            fi                 
              
            log "INFO" "📦 Finished MACS3 call (with input) for $name"
             draw_progress "$name"
       
        fi
    fi
  done
} < "$MAPPING"

printf "\n" >&2
# === Summarize called peaks ===
record_metrics "MACS3 call complete ☎️. Starting summary 🧾 ..."

# === Generate summary ===
print_header  "📑 Generating MACS3 peak summary table..."

SUMMARY_OUT="$PEAK_DIR/macs3_summary.tsv"
echo -e "Sample\tTotal_Peaks\tAvg_Peak_Length\tOutput_Dir\tWarning" > "$SUMMARY_OUT"

log "INFO" "📄 Initializing fresh macs3_summary.tsv"
for xls in "$PEAK_DIR"/*/*_peaks.xls; do
    [[ ! -f "$xls" ]] && continue
    sample=$(basename "$xls" _peaks.xls)
    output_dir=$(dirname "$xls")

    total_peaks=$(grep -v "^#" "$xls" | wc -l)
#avg_len=$(grep -v "^#" "$xls" | awk '{sum += $3 - $2} END {if (NR > 0) printf "%.1f", sum/NR; else print "NA"}')
    avg_len=$(awk '!/^#/ {sum += $3 - $2; n++} END {print (n > 0 ? sprintf("%.1f", sum/n) : "NA")}' "$xls")

    warning="OK"
    if [[ "$total_peaks" -lt 100 ]]; then
        warning="VERY_LOW_PEAKS"
    elif [[ "$total_peaks" -lt 500 ]]; then
        warning="LOW_PEAKS"
    fi

    echo -e "${sample}\t${total_peaks}\t${avg_len}\t${output_dir}\t${warning}" >> "$SUMMARY_OUT"
done

# === Organize output by peak type ===
narrow_files=$(find "$PEAK_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_peaks.narrowPeak" -print)
broad_files=$(find "$PEAK_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_peaks.broadPeak" -print)

if ! $DRY_RUN; then
    if [[ -n "$narrow_files" ]]; then
        mkdir -p "$PEAK_DIR/narrowPeak"
        while IFS= read -r file; do
            [[ -f "$file" ]] || continue
            dest="$PEAK_DIR/narrowPeak/$(basename "$file")"
            [[ "$file" != "$dest" ]] && cp "$file" "$dest"
        done <<< "$narrow_files"
    fi

    if [[ -n "$broad_files" ]]; then
        mkdir -p "$PEAK_DIR/broadPeak"
        while IFS= read -r file; do
            [[ -f "$file" ]] || continue
            dest="$PEAK_DIR/broadPeak/$(basename "$file")"
            [[ "$file" != "$dest" ]] && cp "$file" "$dest"
        done <<< "$broad_files"
    fi
else
    log "INFO" "🧯 [DRY-RUN] Would organize peak files into narrowPeak/ and broadPeak/"
fi

# === Compute FRiP scores ===

FRIP_SUMMARY="$PEAK_DIR/frip_summary.tsv"

# Initialize new temporary file
TMP_FRIP="$FRIP_SUMMARY.tmp"
echo -e "Sample\tFRiP\tReads_in_Peaks\tTotal_Reads\tPeak_File\tConclusion" > "$TMP_FRIP"

# Find peak files from both narrow and broad output dirs
for peak_file in "$PEAK_DIR"/narrowPeak/*_peaks.narrowPeak "$PEAK_DIR"/broadPeak/*_peaks.broadPeak; do
    [[ ! -f "$peak_file" ]] && continue

    # Grab everything BEFORE "_vs_" — that's the IP BAM sample
    ip_sample=$(basename "$peak_file" | sed 's/_vs_.*//; s/_peaks\..*//')
    bam_file="$BAM_DIR/${ip_sample}.bam"


    if [[ ! -f "$bam_file" ]]; then
        log "WARN" "⚠️ BAM file not found for FRiP: $bam_file"
        continue
    fi

    reads_in_peaks=$(bedtools intersect -a "$bam_file" -b "$peak_file" -bed | wc -l)
    total_reads=$(samtools view -c -F 260 "$bam_file")
    
    if [[ "$total_reads" -eq 0 ]]; then
        frip_score="NA"
        conclusion="NO_READS"
    else
        frip_score=$(awk -v r="$reads_in_peaks" -v t="$total_reads" 'BEGIN {printf "%.4f", r/t}')
        
        # Add interpretation
        if (( $(echo "$frip_score < 0.01" | bc -l) )); then
            conclusion="FAIL_FRIP: Poor enrichment"
        elif (( $(echo "$frip_score < 0.05" | bc -l) )); then
            conclusion="LOW_FRIP: Weak enrichment"
        elif (( $(echo "$frip_score < 0.2" | bc -l) )); then
            conclusion="OK_FRIP: Moderate signal"
        else
            conclusion="HIGH_FRIP: Strong signal"
        fi
    fi

    # Only add if not already in TMP
    if ! grep -q "^${ip_sample}\t" "$TMP_FRIP"; then
        echo -e "${ip_sample}\t${frip_score}\t${reads_in_peaks}\t${total_reads}\t${peak_file}\t${conclusion}" >> "$TMP_FRIP"
    fi
done

# Move temp to final summary
mv "$TMP_FRIP" "$FRIP_SUMMARY"

log "INFO" "📈 Annotated FRiP scores written to: $FRIP_SUMMARY"
log "INFO" "📂 Using BAM for FRiP: $bam_file"  # Optional for debugging

elapsed="${SECONDS}s"
record_metrics "PROCESS COMPLETE (Elapsed: $elapsed)"

# === Final message ===
converted=$(find "$PEAK_DIR"/{narrowPeak,broadPeak} -name '*_peaks.*' 2> /dev/null | wc -l)

if [[ -z "$converted" || "$converted" -eq 0 ]]; then
    log "WARN" "⚠️ No converted peak files were found in $PEAK_DIR/narrowPeak or broadPeak"
else
    log "INFO" "🧾 Final organized peak files: $converted"
fi
print_header "📦 MODULE COMPLETE"
log "INFO" "🧾 Final converted peak files: $converted"
log "INFO" "⛰️ 🏁 Peak calling finished: $(date)"
print_header "✅ [MODULE: ${SCRIPT_BASE_NAME}] Finished at $(date '+%F_%H-%M-%S')"
