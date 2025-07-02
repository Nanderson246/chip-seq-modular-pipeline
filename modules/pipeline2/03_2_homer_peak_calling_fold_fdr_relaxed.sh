#!/usr/bin/env bash
# Module: 03_2_homer_peak_calling_fold_fdr_relaxed.sh (flexible metadata support)
# Author: Nancy Anderson (updated with input fallback and sample_type_rules YAML support)

################################################################################
# SOFTWARE REQUIREMENTS (for Docker container build):
#
# Core:
# - HOMER >= 4.11             # For peak calling (findPeaks, makeTagDirectory)
# - Rscript                   # For HOMER-to-BED format conversion
# - bedtools >= 2.29          # For FRiP score calculation
# - samtools >= 1.10          # For BAM read counting
# - yq >= 4.0                 # For parsing YAML rules (Go version: https://github.com/mikefarah/yq)
#
# Optional:
# - bash >= 4.2               # For associative arrays and advanced parsing
# - bc                        # For floating point arithmetic in FRiP scoring
#
# - R >= 4.0                 # Required to run this script
# - Rscript                  # To execute the script from shell
#
# Base R packages used:
# - read.delim
# - write.table
# - readLines / writeLines
# - p.adjust
# - tempfile, sprintf, etc.
#
# Notes:
# - No external CRAN or Bioconductor packages required
# - No internet access or additional installs needed inside the container
################################################################################

################################################################################
# USAGE:
#    bash modules/pipeline2/03_2_homer_peak_calling_fold_fdr_relaxed.sh --style 
#    bash modules/pipeline2/03_2_homer_peak_calling_fold_fdr_relaxed.sh --mapping metadata/mapping_filtered.tsv
#
#Override Style:
#bash modules/pipeline2/03_2_homer_peak_calling_fold_fdr_relaxed.sh \
#  --mapping metadata/mapping_filtered.tsv \
#  --style histone
# 
#USAGE STANDALONE:
# bash 03_2_homer_peak_calling_fold_fdr_relaxed.sh \
#  --mapping metadata/mapping_filtered.tsv \
#  --bam-dir analysis/Renamed_Cleaned \
#  --out-dir analysis/PeakCalling_HOMER \
#  --yaml templates/target_peak_rules_HOMER.yaml 
# You can also override the style if needed:
#  --style histone

################################################################################

set -uo pipefail

# === SCRIPT ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}" 
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
DRY_RUN=false

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



# Default paths (can be overridden by CLI or env)
STYLE="${STYLE:-}"
FDR="0.001"
FOLD_INPUT=4
FOLD_LOCAL=4
BAM_DIR="${BAM_DIR:-$PROJECT_ROOT/analysis/Renamed_Cleaned}"
PEAK_DIR="${PEAK_DIR:-$PROJECT_ROOT/analysis/PeakCalling_HOMER}"
MAPPING="${MAPPING:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
TARGET_YAML="${TARGET_YAML:-$PROJECT_ROOT/templates/target_peak_rules_HOMER.yaml}"
RULES_YAML="${RULES_YAML:-$PROJECT_ROOT/templates/mapping_schema.yaml}"
LOG_FILE="$PEAK_DIR/homer_peak_calling.log"
PERF_LOG="$PEAK_DIR/homer_peak_calling_metrics.log"

mkdir -p "$PEAK_DIR"/{tags,peaks}
 
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
print_header "🧬 Module: ${SCRIPT_NAME}"
log "INFO" "📌 Purpose: Peak calling with HOMER"
log "INFO" "📁 Output:  analysis/PeakCalling_HOMER"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "🗂️ Log file: $LOG_FILE"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "📦 Script version: ${VERSION} (${SCRIPT_NAME})"
print_header "INITIALIZATION COMPLETE"

# === CLI Parser ===
usage() {
  echo "Usage: $0 [OPTIONS]"
  echo ""
  echo "Options:"
  echo "  -s, --style STYLE           Override HOMER style manually (e.g., factor, histone)"
  echo "  -m, --mapping FILE          Metadata file (TSV)"
  echo "  -b, --bam-dir DIR           Directory with input BAMs [default: analysis/Renamed_Cleaned]"
  echo "  -o, --out-dir DIR           Output directory for HOMER results [default: analysis/PeakCalling_HOMER]"
  echo "  -y, --yaml FILE             YAML file for target-to-style mapping [default: templates/target_peak_rules_HOMER.yaml]"
  echo "  --fdr FLOAT                 FDR cutoff for HOMER (default: 0.001)"
  echo "  -F, --fold-input INT        Fold over input cutoff (default: 4)"
  echo "  -L, --fold-local INT        Fold over local background cutoff (default: 4)"
  echo "  -d, --dry-run               Simulate commands without executing"
  echo "  -h, --help                  Show this help message and exit"
  exit 1
}


# === Parse CLI arguments ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dry-run) DRY_RUN=true; shift ;;
        -s|--style) STYLE="$2"; shift 2;;
        -m|--mapping) MAPPING="$2"; shift 2;;
        -b|--bam-dir) BAM_DIR="$2"; shift 2 ;;
        -o|--out-dir) PEAK_DIR="$2"; shift 2 ;;
        -y|--yaml) TARGET_YAML="$2"; shift 2 ;;
        --fdr) FDR="$2"; shift 2;;
        -F|--fold-input) FOLD_INPUT="$2"; shift 2;;
        -L|--fold-local) FOLD_LOCAL="$2"; shift 2;;
        -h|--help) usage;;
        *) log "ERROR" "Unknown option: $1"; usage;;
    esac
done

# === Initial Validation ===
if [[ ! -f "$MAPPING" ]]; then
    log "ERROR" "Mapping file not found: $MAPPING"
    exit 1
fi

HEADER=$(head -n 1 "$MAPPING")
IFS=$'\t' read -ra COLS <<< "$HEADER"

for field in Sample_ID Sample_Type Condition Target; do
    if ! [[ " ${COLS[@]} " =~ " ${field} " ]]; then
        log "ERROR" "Missing required field in mapping file: $field"
        exit 1
    fi
done

# === Validate HOMER style ===
valid_styles=("factor" "atac" "histone" "groseq" "tsr" "dnase" "super" "superhistone" "mC" "damid" "clip")

if [[ -n "$STYLE" && ! " ${valid_styles[*]} " =~ " $STYLE " ]]; then
    log "ERROR" "Invalid HOMER style: $STYLE. Valid styles: ${valid_styles[*]}"
    exit 1
fi

# === Check HOMER, R and YQ availability ===
if ! command -v findPeaks &>/dev/null || ! command -v makeTagDirectory &>/dev/null; then
    log "ERROR" "HOMER tools not found in PATH."
    exit 1
fi

if ! command -v yq &>/dev/null; then
    log "ERROR" "yq YAML parser not found in PATH."
    exit 1
fi

if ! command -v Rscript &>/dev/null; then
    log "ERROR" "Rscript is required for HOMER peak format conversion"
    exit 1
fi

# === Read allowed sample types from YAML ===
IP_TYPES=($(yq '.sample_type_rules.ip_types[]' "$RULES_YAML"))
CONTROL_TYPES=($(yq '.sample_type_rules.control_types[]' "$RULES_YAML"))

log "INFO" "Accepted IP types: ${IP_TYPES[*]}"
log "INFO" "Accepted Control types: ${CONTROL_TYPES[*]}"

# === Setup logging ===
echo "[INFO] Peak calling started: $(date)" > "$LOG_FILE"
echo "[INFO] HOMER style: $STYLE | FDR: $FDR | Fold Input: $FOLD_INPUT | Fold Local: $FOLD_LOCAL" >> "$LOG_FILE"

# === Setup Dry run ===
if $DRY_RUN; then
    log "INFO" "DRY-RUN mode enabled — no files will be modified"
    echo "[INFO] DRY-RUN mode enabled" >> "$LOG_FILE"
fi

# === Begin real processing ===
SECONDS=0
record_metrics "PROCESS START"

# === Load target → style rules from YAML ===
if [[ ! -f "$TARGET_YAML" ]]; then
    log "ERROR" "target_peak_rules_HOMER.yaml not found at $TARGET_YAML"
    exit 1
fi

log "INFO" "Loading HOMER styles from: $TARGET_YAML"

# Function to get style for a target
get_homer_style() {
    local target="$1"
    local style
    
    # Check factor style targets
    if yq e '.factor_style_targets[]' "$TARGET_YAML" | grep -qw "$target"; then
        echo "factor"
        return
    fi
    
    # Check other style categories
    for category in histone_style_targets dnase_style_targets groseq_style_targets \
                   tss_style_targets super_style_targets mc_style_targets; do
        if yq e ".${category}[]" "$TARGET_YAML" | grep -qw "$target"; then
            case "$category" in
                histone_style_targets) echo "histone" ;;
                dnase_style_targets) echo "dnase" ;;
                groseq_style_targets) echo "groseq" ;;
                tss_style_targets) echo "tss" ;;
                super_style_targets) echo "super" ;;
                mc_style_targets) echo "mC" ;;
                *) echo "factor" ;;
            esac
            return
        fi
    done
    
    # Fallback to factor style
    echo "factor"
}

# === Count how many IP samples will be processed ===
total_calls=$(awk -F'\t' -v ip_types="${IP_TYPES[*]}" '
    BEGIN { split(ip_types, a, " "); count = 0 }
    NR > 1 {
        for (i in a) {
            if ($5 ~ "^" a[i]) { count++ }
        }
    }
    END { print count }
' "$MAPPING")
current=0

# === Parse metadata and call peaks ===
tail -n +2 "$MAPPING" | while IFS=$'\t' read -r line; do
    [[ -z "$line" ]] && continue
    
    IFS=$'\t' read -ra FIELDS <<< "$line"
    declare -A meta
    for i in "${!COLS[@]}"; do
        meta["${COLS[$i]}"]="${FIELDS[$i]:-}"
    done

    Sample_ID="${meta[Sample_ID]}"
    Sample_Type="${meta[Sample_Type]}"
    Condition="${meta[Condition]}"
    Target="${meta[Target]}"
    
 # === VALIDATION Sample_Type found in sample ===
    if [[ ! "${Sample_Type}" =~ ^(Input|IgG|Mock|IP(_rep[0-9]+)?|ChIP)$ ]]; then
        log "ERROR" "Invalid Sample_Type '${Sample_Type}' for $Sample_ID (Allowed: Input, IgG, Mock, IP, ChIP)"
        continue
    fi
    # === END VALIDATION ===


    [[ -z "$Sample_ID" || -z "$Sample_Type" || -z "$Condition" || -z "$Target" ]] && {
        log "WARN" "Missing required fields for sample - skipping"
        continue
    }

    # === Flexible Sample Type Checking ===
    is_ip_sample=false
    for type in "${IP_TYPES[@]}"; do
        if [[ "$Sample_Type" == "$type"* ]]; then
            is_ip_sample=true
            break
        fi
    done

    if [[ "$is_ip_sample" != true ]]; then
        log "INFO" "Skipping non-IP sample: $Sample_ID (Type: $Sample_Type)"
        continue
    fi

    # === Build base name ===
    parts=("$Sample_ID")
    for key in Instrument Condition Replicate Sample_Type Target; do
        val="${meta[$key]:-}"
        [[ -n "$val" ]] && parts+=("$(echo "$val" | sed 's/[^a-zA-Z0-9]/_/g')")
    done
    base=$(IFS=_; echo "${parts[*]}")
    bam_path="$BAM_DIR/$base.bam"
    
    if [[ ! -f "$bam_path" ]]; then
        log "ERROR" "BAM file not found: $bam_path"
        continue
    fi

    # === Determine style dynamically if not set by CLI ===
    local_style="$STYLE"
    if [[ -z "$local_style" ]]; then
        local_style=$(get_homer_style "$Target")
        log "INFO" "$base: Detected style '$local_style' from target '$Target'"
    fi

    # === Determine output subfolder and file extension based on style ===
    case "$local_style" in
        factor|dnase|tss|clip)
            SUBDIR="narrow_PEAKS"
            EXT="narrowPeak"
            ;;
        histone|super|superhistone|groseq)
            SUBDIR="broad_PEAKS"
            EXT="broadPeak"
            ;;
        mc)
            SUBDIR="methylation_PEAKS"
            EXT="bed"
            ;;
        *)
            SUBDIR="narrow_PEAKS"
            EXT="narrowPeak"
            ;;
    esac

    mkdir -p "$PEAK_DIR/$SUBDIR"
   # === Try to find best-matching control (Input/IgG/Mock) ===
    input_line=$(awk -F'\t' -v id="$Sample_ID" \
        -v cond="$Condition" \
        -v target="$Target" \
        -v inst="${meta[Instrument]:-}" \
        -v rep="${meta[Replicate]:-}" \
        'BEGIN {OFS=FS}
     $1 != id && ($5 ~ /^(Input|IgG|Mock)$/) && $3 == cond && $7 == target && $6 == inst && $4 == rep {print $0; exit}' "$MAPPING")

# Relax to instrument match only
    if [[ -z "$input_line" ]]; then
      input_line=$(awk -F'\t' -v id="$Sample_ID" \
          -v cond="$Condition" \
          -v target="$Target" \
          -v inst="${meta[Instrument]:-}" \
          'BEGIN {OFS=FS}
       $1 != id && ($5 ~ /^(Input|IgG|Mock)$/) && $3 == cond && $7 == target && $6 == inst {print $0; exit}' "$MAPPING")
    fi

# Relax to just condition + target
    if [[ -z "$input_line" ]]; then
      input_line=$(awk -F'\t' -v id="$Sample_ID" \
          -v cond="$Condition" \
          -v target="$Target" \
          'BEGIN {OFS=FS}
       $1 != id && ($5 ~ /^(Input|IgG|Mock)$/) && $3 == cond && $7 == target {print $0; exit}' "$MAPPING")
    fi

    has_input=false
    input_path=""
    input_tag=""

    if [[ -n "$input_line" ]]; then
        IFS=$'\t' read -ra input_fields <<< "$input_line"
        declare -A input_meta
        for i in "${!COLS[@]}"; do
            input_meta["${COLS[$i]}"]="${input_fields[$i]:-}"
        done

        Input_ID="${input_meta[Sample_ID]}"
        log "INFO" "Found matching control sample: $Input_ID (Type: ${input_meta[Sample_Type]})"  
        
        parts_input=("$Input_ID")
        for key in Instrument Condition Replicate Sample_Type Target; do
            val="${input_meta[$key]:-}"
            [[ -n "$val" ]] && parts_input+=("$(echo "$val" | sed 's/[^a-zA-Z0-9]/_/g')")
        done
        Input_base=$(IFS=_; echo "${parts_input[*]}")
        input_path="$BAM_DIR/$Input_base.bam"

        if [[ -f "$input_path" ]]; then
            has_input=true
        else
            log "WARN" "Control BAM not found: $input_path"
        fi
    fi

    # === Process IP sample ===
    ip_tag="$PEAK_DIR/tags/$base"
    peak_file="$PEAK_DIR/peaks/${base}_peaks.txt"
    output_prefix="$base"
    
    # === SAFEGUARD: Skip if all outputs already exist (peak call + final converted files) ===
    peak_out_base="$PEAK_DIR/$SUBDIR/${base}"
    peak_file_out="${peak_out_base}_peaks.${EXT}"
    bed_out="${peak_out_base}_summits.bed"
    xml_out="${peak_out_base}_peaks.xml"
    if [[ -d "$ip_tag" && -f "$peak_file" && -f "$peak_file_out" && -f "$bed_out" && -f "$xml_out" ]]; then
        log "INFO" "✅ Skipping HOMER: All outputs exist (including XML)  for $base → $peak_file_out and $bed_out"
        continue
    fi

    
    if $DRY_RUN; then
        log "INFO" "[DRY-RUN] Would process IP sample: $base"
        [[ "$has_input" == true ]] && {
            output_prefix="${base}_vs_${Input_base}"
            log "INFO" "[DRY-RUN] Would use control: $Input_base"
        }
        continue
    fi
 
    # Create tag directory if needed
    if [[ ! -d "$ip_tag" ]]; then
        mkdir -p "$ip_tag"
        if ! makeTagDirectory "$ip_tag" "$bam_path" 2>> "$LOG_FILE"; then
            log "ERROR" "Failed to create tag directory for $base"
            continue
        fi
   else
    log "INFO" "🗂️ Tag directory already exists for $base: skipping makeTagDirectory"
   fi     
        

    if [[ "$has_input" == true ]]; then
        input_tag="$PEAK_DIR/tags/$Input_base"
        
        if [[ ! -d "$input_tag" ]]; then
            mkdir -p "$input_tag"
            if ! makeTagDirectory "$input_tag" "$input_path" 2>> "$LOG_FILE"; then
                log "ERROR" "Failed to create tag directory for control $Input_base"
                continue
            fi
        else
            log "INFO" "🗂️ Control tag directory already exists for $Input_base: skipping makeTagDirectory"
        fi     

        log "INFO" "Calling HOMER peaks with Control for $base (Style: $local_style)"
        if ! findPeaks "$ip_tag" -style "$local_style" -i "$input_tag" \
            -F "$FOLD_INPUT" -L "$FOLD_LOCAL" -fdr "$FDR" -o "$peak_file" 2>> "$LOG_FILE"; then
            log "ERROR" "HOMER peak calling failed for $base with control"
            continue
        fi
        output_prefix="${base}_vs_${Input_base}"
    else
        log "WARN" "No Control found for $base — calling HOMER peaks without control"
        if ! findPeaks "$ip_tag" -style "$local_style" \
            -F "$FOLD_INPUT" -L "$FOLD_LOCAL" -fdr "$FDR" -o "$peak_file" 2>> "$LOG_FILE"; then
            log "ERROR" "HOMER peak calling failed for $base without control"
            continue
        fi
        log "INFO" "Finished HOMER call (with input) for $base"
    fi

    # === Post-process peaks ===
    mkdir -p "$PEAK_DIR/$SUBDIR" 
    peak_out_base="$PEAK_DIR/$SUBDIR/$output_prefix"
    bed_out="${peak_out_base}_summits.bed"
    peak_file_out="${peak_out_base}_peaks.${EXT}"
    xml_out="${peak_out_base}_peaks.xml"

    if [[ ! -f "$peak_file_out" || ! -f "$bed_out" ]]; then
        log "INFO" "Converting HOMER peaks to BED/$EXT for $base"
        if ! Rscript modules/pipeline2/03_2_homer_to_peakFormat.R "$peak_file" "$peak_out_base" "$base" "$local_style" 2>> "$LOG_FILE"; then
            log "ERROR" "Failed to convert peaks for $base"
            continue
        fi

        if [[ -f "$peak_file_out" ]]; then
            sort -k1,1 -k2,2n "$peak_file_out" -o "${peak_file_out}.tmp" && \
            mv "${peak_file_out}.tmp" "$peak_file_out"
        fi
        draw_progress "$base"
    fi
done

echo "" >&2
echo -e "\n${GREEN}✅ All peak calls completed${NC}"

# === Summary ===
log "INFO" "Generating HOMER peak summary table..."
SUMMARY_OUT="$PEAK_DIR/homer_summary.tsv"
echo -e "Sample\tTotal_Peaks\tPeak_Size\tTotal_Tags\tTags_in_Peaks\tIP_Efficiency\tWarning" > "$SUMMARY_OUT"

for f in "$PEAK_DIR"/peaks/*_peaks.txt; do
    [[ ! -f "$f" ]] && continue
    sample=$(basename "$f" _peaks.txt)

    total_peaks=$(awk '/# total peaks/ {print $5}' "$f")
    peak_size=$(awk '/# peak size/ {print $5}' "$f")
    total_tags=$(awk '/# Total tags =/ {print $5}' "$f")
    tags_in_peaks=$(awk '/# Total tags in peaks/ {print $7}' "$f")
    ip_eff=$(awk '/# Approximate IP efficiency/ {print $6}' "$f" | tr -d '%')

    warning="OK"
    if [[ -z "$ip_eff" ]]; then
        warning="MISSING_DATA"
    elif (( $(echo "$ip_eff < 0.3" | bc -l) )); then
        warning="CRITICAL_IP"
    elif (( $(echo "$ip_eff < 1.0" | bc -l) )); then
        warning="LOW_IP"
    fi

    echo -e "${sample}\t${total_peaks}\t${peak_size}\t${total_tags}\t${tags_in_peaks}\t${ip_eff}%\t${warning}" >> "$SUMMARY_OUT"
done

{ head -n 1 "$SUMMARY_OUT"; tail -n +2 "$SUMMARY_OUT" | sort -k1,1; } > "${SUMMARY_OUT}.tmp" && \
mv "${SUMMARY_OUT}.tmp" "$SUMMARY_OUT"

log "INFO" "HOMER summary saved to: $SUMMARY_OUT"

# === Compute FRiP scores for HOMER peaks ===
FRIP_SUMMARY="$PEAK_DIR/frip_summary.tsv"
echo -e "Sample\tFRiP\tReads_in_Peaks\tTotal_Reads\tPeak_File\tConclusion" > "$FRIP_SUMMARY"

for peak_file in "$PEAK_DIR"/narrow_PEAKS/*_peaks.narrowPeak "$PEAK_DIR"/broad_PEAKS/*_peaks.broadPeak; do
    [[ ! -f "$peak_file" ]] && continue

     # Extract base sample name from peak file (before _vs_ and _peaks)
    sample=$(basename "$peak_file" | sed 's/_vs_.*//; s/_peaks\..*//')
    bam_file="$BAM_DIR/${sample}.bam"
    
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

    echo -e "${sample}\t${frip_score}\t${reads_in_peaks}\t${total_reads}\t${peak_file}\t${conclusion}" >> "$FRIP_SUMMARY"
done

log "INFO" "📈 Annotated HOMER FRiP scores written to: $FRIP_SUMMARY"



elapsed="${SECONDS}s"
record_metrics "PROCESS COMPLETE (Elapsed: $elapsed)"


converted=$(find "$PEAK_DIR"/narrow_PEAKS "$PEAK_DIR"/broad_PEAKS -name '*_peaks.*' 2>/dev/null | wc -l)

if [[ -z "$converted" || "$converted" -eq 0 ]]; then
    log "WARN" "⚠️ No converted peak files were found in $PEAK_DIR/narrowPeak or broadPeak"
else
    log "INFO" "🧾 Final organized peak files: $converted"
fi
raw=$(find "$PEAK_DIR"/peaks -name '*_peaks.txt' 2>/dev/null | wc -l)

# === Final message ===
print_header "📦 MODULE COMPLETE"
log "INFO" "🧾 Raw HOMER peak summary files: $raw"
log "INFO" "📦 Final converted peak files: $converted"
log "INFO" "⛰️ 🏁 Peak calling finished: $(date)"
print_header "✅ [MODULE: ${SCRIPT_BASE_NAME}] Finished at $(date '+%F_%H-%M-%S')"
