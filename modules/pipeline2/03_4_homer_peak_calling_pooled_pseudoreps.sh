#!/usr/bin/env bash
# Module: 03_4_homer_peak_calling_pooled_pseudoreps.sh
# Author: Nancy Anderson
# Module: HOMER_peak_calling (Pooled/Pseudoreplicates)
# Description: Call peaks with HOMER findPeaks on pooled ChIP-seq BAMs and pseudoreplicates.

################################################################################
# SOFTWARE REQUIREMENTS (for Docker container build):
#
# Required Tools:
# - bash >= 4.2                 # For associative arrays and CLI parsing
# - HOMER (>= 4.11)             # Includes findPeaks, makeTagDirectory
# - yq >= 4.0                   # For parsing YAML target/control styles
# - samtools >= 1.10            # For FRiP scoring (read counting from BAM)
# - bedtools >= 2.29            # For FRiP scoring (peak overlaps)
# - bc                          # For floating point math in FRiP and IP efficiency
# - Rscript                     # For HOMER peak format conversion (uses custom R script)
#
# Notes:
# - HOMER must be installed and available in PATH
# - Assumes a custom R script for HOMER-to-peak format conversion
# - No reference genome or annotation required
################################################################################

#############################################################################
#              Identifies groups by Condition+Target+Instrument and uses matching controls if available.
# Usage in pipeline and inside package:
# bash modules/pipeline2/03_4_homer_peak_calling_pooled_pseudoreps.sh -m metadata/mapping_filtered.tsv
# Output:      analysis/PeakCalling_HOMER (tag directories, raw peak files, converted .narrowPeak/.broadPeak in subfolders, and summary)
# STANDALONE USAGE:
#   bash <pathway>/03_4_HOMER_peak_calling_pooled_pseudoreps.sh [--mapping <metadata.tsv>] [--pooled-dir <dir>] [--pseudo-dir <dir>] [--out-dir <dir>] [--yaml <target_styles.yaml>] [--style <override>] [--dry-run]
#############################################################################
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
    local level="$1"; local msg="$2"; local ts=$(date '+%Y-%m-%d %H:%M:%S')
    case "$level" in
        ERROR) echo -e "${RED}[${ts}] ERROR: ${msg}${NC}" >&2 ;;
        WARN)  echo -e "${YELLOW}[${ts}] WARNING: ${msg}${NC}" >&2 ;;
        INFO)  echo -e "${BLUE}[${ts}] INFO: ${msg}${NC}" ;;
        *)     echo "[${ts}] ${msg}" ;;
    esac
}

print_header() {
    echo "========================================"
    echo "$1"
    echo "========================================"
}


# === Default parameters ===
STYLE_OVERRIDE="${STYLE_OVERRIDE:-}"
FDR="0.001"
FOLD_INPUT=4
FOLD_LOCAL=4
BAM_DIR="${BAM_DIR:-$PROJECT_ROOT/analysis/Renamed_Cleaned}"
POOLED_DIR="${POOLED_DIR:-$PROJECT_ROOT/analysis/Pooled_BAMs}"
PSEUDO_DIR="${PSEUDO_DIR:-$PROJECT_ROOT/analysis/Pseudoreplicates}"
PEAK_DIR="${PEAK_DIR:-$PROJECT_ROOT/analysis/PeakCalling_HOMER_pool_pseudo}"
MAPPING="${MAPPING:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
TARGET_YAML="${TARGET_YAML:-$PROJECT_ROOT/templates/target_peak_rules_HOMER.yaml}"
RULES_YAML="${RULES_YAML:-$PROJECT_ROOT/templates/mapping_schema.yaml}"
LOG_FILE="$PEAK_DIR/homer_peak_calling_pooled.log"
PERF_LOG="$PEAK_DIR/homer_peak_calling_pooled_metrics.log"

# Ensure output directories exist
mkdir -p "$PEAK_DIR"/{tags,peaks}

# === Performance logging function ===
record_metrics() {
    local msg="$1"; local ts=$(date +%Y-%m-%dT%H:%M:%S)
    local cpu=$(ps -p $$ -o %cpu= | xargs)
    local mem=$(ps -p $$ -o %mem= | xargs)
    local disk=$(df -h "$PEAK_DIR" | tail -1)
    echo "[$ts] $msg" >> "$PERF_LOG"
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
print_header "🔬  Module: ${SCRIPT_NAME}"
log "INFO" "📌 Purpose: Peak calling with HOMER on pooled and pseudoreplicate BAMs"
log "INFO" "📁 Output:  $PEAK_DIR"
log "INFO" "🚀 Dry-run mode: $DRY_RUN"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "🗂️ Log file: $LOG_FILE"
log "INFO" "📦 Script version: ${VERSION} (${SCRIPT_NAME})"
print_header "INITIALIZATION COMPLETE"

# === CLI Parser ===

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d | --dry-run)
            DRY_RUN=true
            shift
            ;;
        -s | --style)
            STYLE_OVERRIDE="$2"
            shift 2
            ;;
        -m | --mapping)
            MAPPING="$2"
            shift 2
            ;;
        -o | --out-dir)
            PEAK_DIR="$2"
            LOG_FILE="$PEAK_DIR/homer_peak_calling_pooled.log"
            PERF_LOG="$PEAK_DIR/homer_peak_calling_pooled_metrics.log"
            shift 2
            ;;
        -p | --pooled-dir)
            POOLED_DIR="$2"
            shift 2
            ;;
        -p | --pseudo-dir)
            PSEUDO_DIR="$2"
            shift 2
            ;;
        --yaml)
            TARGET_YAML="$2"
            shift 2
            ;;
        --rules)
            RULES_YAML="$2"
            shift 2
            ;;

        -h | --help)
            echo "Usage: $0 [options]"
            echo "  -m, --mapping FILE       Metadata TSV file (default: metadata/mapping_filtered.tsv)"
            echo "  -o, --out-dir DIR        Output directory for HOMER results (default: analysis/PeakCalling_HOMER)"
            echo "  -p, --pooled-dir DIR     Directory of pooled BAMs (default: analysis/Pooled_BAMs)"
            echo "  -s, --pseudo-dir DIR     Directory of pseudoreplicates (default: analysis/Pseudoreplicates)"
            echo "      --yaml FILE          YAML file for target→HOMER style mapping (default: templates/target_peak_rules_HOMER.yaml)"
            echo "      --rules FILE         YAML file for sample type rules (default: templates/mapping_schema.yaml)"
            echo "      --style NAME         Override HOMER peak calling style (e.g., histone, factor)"
            echo "  -d, --dry-run            Simulate commands without executing"
            echo "  -h, --help               Show this help message"
            exit 0
            ;;
        *)
            log "ERROR" "🤷‍♀️ Unknown option: $1"
            exit 1
            ;;
    esac
done

log "[INFO]" "⛰️ Peak calling started: $(date)" > "$LOG_FILE"
if $DRY_RUN; then
    log "[INFO]""🧯 DRY-RUN mode enabled — no files will be modified" | tee -a "$LOG_FILE"
fi

# === Check HOMER, R script and YQ availability ==
if ! command -v findPeaks &>/dev/null || ! command -v makeTagDirectory &>/dev/null; then
    log "ERROR" "❌ HOMER tools 🛠️ (findPeaks/makeTagDirectory) not found 😢 in PATH." | tee -a "$LOG_FILE"
    exit 1
fi

if ! command -v yq &>/dev/null; then
    log "ERROR" "❌ yq YAML parser not found 😢 in PATH." | tee -a "$LOG_FILE"
    exit 1
fi

if ! command -v Rscript &>/dev/null; then
    log "ERROR" "❌ Rscript is required for HOMER peak ⛰️ =>🏔️format conversion" | tee -a "$LOG_FILE"
    exit 1
fi


# Load sample type rules from schema
IP_TYPES=($(yq e '.sample_type_rules.ip_types[]' "$RULES_YAML"))
CONTROL_TYPES=($(yq e '.sample_type_rules.control_types[]' "$RULES_YAML"))

# Function to get HOMER style from target mapping YAML
get_homer_style() {
    local target="$1"
    local style
    
    # First check factor style targets
    if yq e '.factor_style_targets[]' "$TARGET_YAML" | grep -qw "$target"; then
        echo "factor"
        return
    fi
    
    # Then check other categories
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
    
    # Fallback based on target name
    if [[ "$target" =~ H3K|H4K ]]; then
        echo "histone"
    else
        echo "factor"
    fi
}

# === Group samples by Condition+Target+Instrument ===
declare -A group_count
declare -A group_cond group_target group_instrument
IFS=$'\t' read -r -a HEADER < <(head -n1 "$MAPPING")

# Find column indices dynamically
for i in "${!HEADER[@]}"; do
    case "${HEADER[$i]}" in
        "Sample_Type") idxSampleType=$i ;;
        "Condition") idxCondition=$i ;;
        "Instrument") idxInstrument=$i ;;
        "Target") idxTarget=$i ;;
    esac
done

while IFS=$'\t' read -r -a cols || [[ -n "${cols[*]}" ]]; do
    [[ "${cols[0]}" == "" || "${cols[0]}" == "Sample_ID" ]] && continue

    # Check if IP sample
    is_ip=false
    for t in "${IP_TYPES[@]}"; do
        if [[ "${cols[idxSampleType]}" == "$t"* ]]; then is_ip=true; break; fi
    done
    $is_ip || continue

    cond="${cols[idxCondition]}"; target="${cols[idxTarget]}"; inst="${cols[idxInstrument]}"
    group="${cond}"
    [[ -n "$target" ]] && group+="_${target}"
    [[ -n "$inst" ]] && group+="_${inst}"
    group_count["$group"]=$(( ${group_count["$group"]:-0} + 1 ))
    group_cond["$group"]="$cond"; group_target["$group"]="$target"; group_instrument["$group"]="$inst"
done < <(tail -n +2 "$MAPPING")


SECONDS=0
record_metrics "PROCESS START"

total_calls=0
for g in "${!group_count[@]}"; do
    (( group_count["$g"] >= 2 )) && (( total_calls += 3 )) # pooled + 2 pseudos
done
current=0

for group in "${!group_count[@]}"; do
    count=${group_count["$group"]}
    if [[ $count -lt 2 ]]; then
        log "INFO" "⏭️Skipping group '$group' (only $count sample)" | tee -a "$LOG_FILE"
        continue
    fi
    cond="${group_cond[$group]}"; inst="${group_instrument[$group]}"; targ="${group_target[$group]}"
    log "INFO" "🧠 Processing group: $group (n=$count samples)" | tee -a "$LOG_FILE"

    # Determine HOMER peak calling style for this group's target
    if [[ -n "$STYLE_OVERRIDE" ]]; then
        STYLE="$STYLE_OVERRIDE"
        log "INFO" "$group: Using user-specified HOMER style '$STYLE_OVERRIDE'" | tee -a "$LOG_FILE"
    else
        STYLE=$(get_homer_style "$targ")
        log "INFO" "$group: 👀 Detected HOMER peak calling style '$STYLE'🎨 for target 🎯 '$targ'" | tee -a "$LOG_FILE"
    fi

    # Decide output subdirectory and file extension based on style
    case "$STYLE" in
        factor | dnase | tss | clip)
            SUBDIR="narrow_PEAKS"
            EXT="narrowPeak"
            ;;
        histone | super* | groseq)
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

    # Locate matching control (Input/IgG/Mock with same Condition, Instrument, Target)
    input_line=$(awk -F'\t' -v c="$cond" -v i="$inst" -v t="$targ" \
                 '$5 ~ /^(Input|IgG|Mock)$/ && $3 == c && $6 == i && $7 == t {print; exit}' "$MAPPING")
    has_input=false; Input_base=""; input_path=""
    if [[ -n "$input_line" ]]; then
        IFS=$'\t' read -r -a input_fields <<< "$input_line"
        declare -A input_meta
        for idx in "${!HEADER[@]}"; do input_meta["${HEADER[$idx]}"]="${input_fields[$idx]}"; done
        Input_ID="${input_meta[Sample_ID]}"
        parts_input=("$Input_ID")
        for key in Instrument Condition Replicate Sample_Type Target; do
            val="${input_meta[$key]}"; [[ -n "$val" ]] && parts_input+=("$(echo "$val" | sed 's/[^a-zA-Z0-9]/_/g')")
        done
        Input_base=$(IFS=_; echo "${parts_input[*]}")
        input_path="$BAM_DIR/${Input_base}.bam"
        if [[ -f "$input_path" ]]; then
            has_input=true
            log "INFO" "📍 Found matching control: $Input_base" | tee -a "$LOG_FILE"
        else
            log "WARN" "💨 Control BAM not found: $input_path" | tee -a "$LOG_FILE"
        fi
    else
        log "WARN" "🧩 No matching control found for group $group" | tee -a "$LOG_FILE"
    fi

    # Define base names for pooled and pseudorep outputs
    base_pooled="pooled_${group}"
    base_pseudo1="${group}_pseudo1"
    base_pseudo2="${group}_pseudo2"

    # Prepare input tag directory if control exists
    if [[ "$has_input" == true ]]; then
        input_tag="$PEAK_DIR/tags/$Input_base"
        if $DRY_RUN; then
            log "INFO" "🧯 [DRY-RUN] Would create tag 📌 directory for control: $input_tag" | tee -a "$LOG_FILE"
        else
            mkdir -p "$input_tag"
            if [[ -z "$(ls -A "$input_tag")" ]]; then
                makeTagDirectory "$input_tag" "$input_path" 2>> "$LOG_FILE" || {
                    log "ERROR" "❌ Failed to create tag 📌 directory for control $Input_base" | tee -a "$LOG_FILE"
                    continue
                }
            fi
        fi
    fi

   # Process pooled and pseudoreplicate BAMs
    for type in pooled pseudo1 pseudo2; do
        if [[ "$type" == "pooled" ]]; then
            base="$base_pooled"
            bam_path="$POOLED_DIR/${base}.bam"
        elif [[ "$type" == "pseudo1" ]]; then
            base="$base_pseudo1"
            bam_path="$PSEUDO_DIR/${group}/${base}.bam"
        elif [[ "$type" == "pseudo2" ]]; then
            base="$base_pseudo2"
            bam_path="$PSEUDO_DIR/${group}/${base}.bam"
        else
            log "ERROR" "Unknown replicate type: $type"
            continue
        fi


        # Skip if BAM doesn't exist
        if [[ ! -f "$bam_path" ]]; then
            log "ERROR" "Missing BAM for $base ($bam_path)" | tee -a "$LOG_FILE"
            continue
        fi




        # Prepare tag directory for this IP
        ip_tag="$PEAK_DIR/tags/$base"
        if $DRY_RUN; then
            log "INFO" "🧯[DRY-RUN] Would create tag 📌 directory: $ip_tag" | tee -a "$LOG_FILE"
        else
            mkdir -p "$ip_tag"
            makeTagDirectory "$ip_tag" "$bam_path" 2>> "$LOG_FILE" || {
                log "ERROR" "❌ Failed to create tag 📌 directory for $base" | tee -a "$LOG_FILE"
                continue
            }
        fi

        # Set output peak file path
        peak_txt="$PEAK_DIR/peaks/${base}_peaks.txt"
        output_prefix="$base"
        [[ "$has_input" == true ]] && output_prefix="${base}_vs_${Input_base}"
        
       # === Skip if all outputs already exist (peak call + final converted files) ===
        peak_out_base="$PEAK_DIR/$SUBDIR/$output_prefix"
        peak_out_file="${peak_out_base}_peaks.${EXT}"
        bed_out="${peak_out_base}_summits.bed"
        xml_out="${peak_out_base}_peaks.xml"

  # === Skip if all outputs already exist ===

        if [[ -d "$ip_tag" && -f "$peak_txt" && -f "$peak_out_file" && -f "$bed_out" && -f "$xml_out" ]]; then
            log "INFO" "✅ Skipping $base — All outputs exist for peak + conversion" | tee -a "$LOG_FILE"
            continue
        fi
            
        # Call peaks with or without control
        if $DRY_RUN; then
            if [[ "$has_input" == true ]]; then
                log "INFO" "[DRY-RUN] Would call peaks: findPeaks $ip_tag -style $STYLE -i $input_tag -F $FOLD_INPUT -L $FOLD_LOCAL -fdr $FDR -o $peak_txt" | tee -a "$LOG_FILE"
            else
                log "INFO" "[DRY-RUN] Would call peaks: findPeaks $ip_tag -style $STYLE -F $FOLD_INPUT -L $FOLD_LOCAL -fdr $FDR -o $peak_txt" | tee -a "$LOG_FILE"
            fi
        else
            if [[ "$has_input" == true ]]; then
                log "INFO" "⛰️ Calling peaks with control for $base" | tee -a "$LOG_FILE"
                findPeaks "$ip_tag" -style "$STYLE" -i "$input_tag" -F "$FOLD_INPUT" -L "$FOLD_LOCAL" -fdr "$FDR" -o "$peak_txt" 2>> "$LOG_FILE" || {
                    log "ERROR" "❌ Peak calling ⛰️failed for $base with control" | tee -a "$LOG_FILE"
                    continue
                }
            else
                log "WARN" "⚠️ Calling peaks ⛰️ without control for $base" | tee -a "$LOG_FILE"
                findPeaks "$ip_tag" -style "$STYLE" -F "$FOLD_INPUT" -L "$FOLD_LOCAL" -fdr "$FDR" -o "$peak_txt" 2>> "$LOG_FILE" || {
                    log "ERROR" "❌ Peak calling ⛰️ failed for $base without control" | tee -a "$LOG_FILE"
                    continue
                }
            fi
        fi

        # Convert peak files
        peak_out_base="$PEAK_DIR/$SUBDIR/$output_prefix"
        bed_out="${peak_out_base}_summits.bed"
        peak_out_file="${peak_out_base}_peaks.${EXT}"

        # Ensure output subdirectory exists *before* conversion
        mkdir -p "$PEAK_DIR/$SUBDIR"

        if $DRY_RUN; then
            log "INFO" "🧯[DRY-RUN] Would convert peaks ⛰️ to 🏔️ $EXT and BED format for $base" | tee -a "$LOG_FILE"
       else
            log "INFO" "♻️ Converting peaks ⛰️ for $base to  🏔️ $EXT format" | tee -a "$LOG_FILE"
            draw_progress "$base"
            printf "\n" >&2

            Rscript modules/pipeline2/03_2_homer_to_peakFormat.R "$peak_txt" "$peak_out_base" "$base" "$STYLE" 2>> "$LOG_FILE" || {
                 log "ERROR" "❌ Failed to convert peaks ⛰️ for $base" | tee -a "$LOG_FILE"
                 continue
    }

           if [[ -f "$peak_out_file" ]]; then
               sort -k1,1 -k2,2n "$peak_out_file" -o "${peak_out_file}.tmp" && mv "${peak_out_file}.tmp" "$peak_out_file"
           fi
       fi

    done
done

# Generate summary
if ! $DRY_RUN; then
    log "INFO" "🧾 Generating peak calling summary" | tee -a "$LOG_FILE"
    SUMMARY_FILE="$PEAK_DIR/homer_pooled_summary.tsv"
    echo -e "Sample\tTotal_Peaks\tType\tPeak_Size\tTotal_Tags\tTags_in_Peaks\tIP_Efficiency\tWarning" > "$SUMMARY_FILE"
    
    for f in "$PEAK_DIR"/peaks/pooled_*_peaks.txt; do
        [[ -f "$f" ]] || continue
        sample=$(basename "$f" _peaks.txt)
        
        # Classify sample type
if [[ "$sample" == *pseudo1* ]]; then
    type="Pseudorep1"
elif [[ "$sample" == *pseudo2* ]]; then
    type="Pseudorep2"
elif [[ "$sample" == *pooled_* && "$sample" != *pseudo* ]]; then
    type="Pooled"
else
    type="Unknown"
fi

               
  # === Summary ===
log "INFO" "Generating HOMER peak summary table..."
SUMMARY_OUT="$PEAK_DIR/homer_summary.tsv"             
               
        # Extract metrics from HOMER output
        total_peaks=$(awk '/# total peaks/ {print $5}' "$f")
        peak_size=$(awk '/# peak size/ {print $5}' "$f")
        total_tags=$(awk '/# Total tags =/ {print $5}' "$f")
        tags_in_peaks=$(awk '/# Total tags in peaks/ {print $7}' "$f")
        ip_eff=$(awk '/# Approximate IP efficiency/ {print $6}' "$f" | tr -d '%')
        
  
        # Determine warning level
        warning="OK"
        if [[ -z "$ip_eff" || "$ip_eff" == "NaN" ]]; then
            warning="MISSING_DATA"
        elif (( $(echo "$ip_eff < 0.3" | bc -l) )); then
            warning="CRITICAL_IP"
        elif (( $(echo "$ip_eff < 1.0" | bc -l) )); then
            warning="LOW_IP"
        fi
        
        echo -e "${sample}\t${type}\t${total_peaks:-0}\t${peak_size:-NA}\t${total_tags:-0}\t${tags_in_peaks:-0}\t${ip_eff:-NA}%\t${warning}" >> "$SUMMARY_FILE"
    done
    
    # Sort summary
    { head -n1 "$SUMMARY_FILE"; tail -n +2 "$SUMMARY_FILE" | sort -k1,1; } > "${SUMMARY_FILE}.tmp" && mv "${SUMMARY_FILE}.tmp" "$SUMMARY_FILE"
    log "INFO" "Summary saved to $SUMMARY_FILE" | tee -a "$LOG_FILE"
fi
        
 
# === Compute FRiP scores for HOMER peaks ===
FRIP_SUMMARY="$PEAK_DIR/frip_summary.tsv"
echo -e "Sample\tType\tFRiP\tReads_in_Peaks\tTotal_Reads\tPeak_File\tConclusion" > "$FRIP_SUMMARY"

for peak_file in "$PEAK_DIR"/narrow_PEAKS/*_peaks.narrowPeak "$PEAK_DIR"/broad_PEAKS/*_peaks.broadPeak; do
    [[ ! -f "$peak_file" ]] && continue

    # Full peak name (e.g., pooled_WT_G4_HiSeq_3000_pseudo1_vs_SRR..._peaks.narrowPeak)
    sample=$(basename "$peak_file" | sed 's/_peaks\..*//')

    # Remove _vs_* to isolate BAM base name
    bam_base=$(echo "$sample" | sed 's/_vs_.*//')

    # Determine BAM path and type
    if [[ "$bam_base" == *_pseudo1 || "$bam_base" == *_pseudo2 ]]; then
        group="${bam_base%_pseudo*}"
        bam_file="$PSEUDO_DIR/$group/${bam_base}.bam"
        type=$( [[ "$bam_base" == *_pseudo1 ]] && echo "Pseudorep1" || echo "Pseudorep2" )
    else
        bam_file="$POOLED_DIR/${bam_base}.bam"
        type="Pooled"
    fi

    if [[ ! -f "$bam_file" ]]; then
        echo -e "${sample}\t${type}\tNA\tNA\tNA\t$peak_file\tMISSING_BAM" >> "$FRIP_SUMMARY"
        continue
    fi

    log "INFO" "📊 Calculating FRiP for $sample using $bam_file"
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

    echo -e "${sample}\t${type}\t${frip_score}\t${reads_in_peaks}\t${total_reads}\t${peak_file}\t${conclusion}" >> "$FRIP_SUMMARY"
done

log "INFO" "📈 Annotated HOMER FRiP scores written to: $FRIP_SUMMARY"
log "INFO" "📂 Using BAM for FRiP: $bam_file"  # Optional for debugging
                   
elapsed="${SECONDS}"
record_metrics "PROCESS COMPLETE (Elapsed: $elapsed)"

# === Final message ===
print_header "📦 MODULE COMPLETE"
log "INFO" "🏁 Peak calling completed at $(date) (Elapsed: $elapsed)" | tee -a "$LOG_FILE"
log "INFO" "🧮 Final output: $(find "$PEAK_DIR" -name '*_peaks.*' | wc -l) peak files" | tee -a "$LOG_FILE"
log "INFO" "📊 FRiP summary saved to: $PEAK_DIR/frip_summary.tsv" | tee -a "$LOG_FILE"
print_header "✅ [MODULE: ${SCRIPT_BASE_NAME}] Finished at $(date '+%F_%H-%M-%S')"
