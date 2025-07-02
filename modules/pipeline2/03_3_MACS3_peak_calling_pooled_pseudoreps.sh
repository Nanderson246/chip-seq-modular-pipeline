#!/usr/bin/env bash
# Module: 03_3_MACS3_peak_calling_pooled_pseudoreps.sh
# Author: Nancy Anderson
#MACS3_peak_calling (Pooled/Pseudoreplicates)
# Description: Call peaks with MACS3 on merged (pooled) ChIP-seq BAMs and their pseudoreplicates.

################################################################################
# SOFTWARE REQUIREMENTS (for Docker container build):
#
# Required Tools:
# - bash >= 4.2                # For associative arrays and reliable parsing
# - MACS3 == 3.0.3             # For peak calling (callpeak)
# - bedtools >= 2.29           # For FRiP score calculations (intersect)
# - samtools >= 1.10           # For read counting from BAMs
# - yq >= 4.0                  # For YAML parsing (Go version: https://github.com/mikefarah/yq)
# - bc                         # For floating point math in FRiP and thresholds
#
# Notes:
# - No R or Python dependencies required
# - No external genome FASTA file required — only genome size keyword (e.g. "hs", "mm")
# - MACS3 must be installed with `pip install macs3` (Python 3.6+)
################################################################################

#############################################################################
# Identifies groups by Condition+Target+Instrument and uses matching Input/IgG/Mock controls if available.
# Output:      analysis/PeakCalling_MACS3 (subdirectories per peak set, with broad/narrow peak files copied to summary folders)
#
# Usage in pipeline and inside package:
# bash modules/pipeline2/03_3_MACS3_peak_calling_pooled_pseudoreps.sh -m metadata/mapping_filtered.tsv
#
#STANDALONE USAGE outside package:
#   bash modules/pipeline2/03_3_MACS3_peak_calling_pooled_pseudoreps.sh [--mapping <metadata.tsv>] [--pooled-dir <dir>] [--pseudo-dir <dir>] [--out-dir <dir>] [--yaml <target_peak_rules.yaml>] [--force-broad] [--dry-run]
#############################################################################

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
        WARN) echo -e "${YELLOW}[${timestamp}] WARNING: ${message}${NC}" >&2 ;;
        INFO) echo -e "${BLUE}[${timestamp}] INFO: ${message}${NC}" ;;
        *) echo "[${timestamp}] ${message}" ;;
    esac
}

print_header() {
    echo "========================================"
    echo "$1"
    echo "========================================"
}

# === Default parameters ===
BAM_DIR="${BAM_DIR:-$PROJECT_ROOT/analysis/Renamed_Cleaned}"
POOLED_DIR="${POOLED_DIR:-$PROJECT_ROOT/analysis/Pooled_BAMs}"
PSEUDO_DIR="${PSEUDO_DIR:-$PROJECT_ROOT/analysis/Pseudoreplicates}"
PEAK_DIR="${PEAK_DIR:-$PROJECT_ROOT/analysis/PeakCalling_MACS3_pool_pseudo}"
MAPPING="${MAPPING:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
TARGET_YAML="${TARGET_YAML:-$PROJECT_ROOT/templates/target_peak_rules.yaml}"
REF_PREFIX="${REF_PREFIX:-hg38}"
LOG_FILE="$PEAK_DIR/macs3_peak_calling_pooled.log"
PERF_LOG="$PEAK_DIR/macs3_peak_calling_pooled_metrics.log"
mkdir -p "$PEAK_DIR"
: > "$LOG_FILE"
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
print_header "🧬 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Purpose: Peak calling with MACS3 on pooled and pseudoreplicate BAMs"
log "INFO" "📥 Input folder: $BAM_DIR ,$POOLED_DIR,$PSEUDO_DIR "
log "INFO" "📁 Output:  $PEAK_DIR"
log "INFO" "🕒 Start time: $(date '+%F %T')"
print_header "INITIALIZATION COMPLETE"

# === CLI Parser ===
DRY_RUN=false
FORCE_BROAD=false
while [[ $# -gt 0 ]]; do
    case "$1" in
        -m | --mapping)
            MAPPING="$2"
            shift 2
            ;;
        -o | --out-dir)
            PEAK_DIR="$2"
            LOG_FILE="$PEAK_DIR/macs3_peak_calling_pooled.log"
            PERF_LOG="$PEAK_DIR/macs3_peak_calling_pooled_metrics.log"
            shift 2
            ;;
        -p | --pooled-dir)
            POOLED_DIR="$2"
            shift 2
            ;;
        -s | --pseudo-dir)
            PSEUDO_DIR="$2"
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
        -d | --dry-run)
            DRY_RUN=true
            shift
            ;;
        -h | --help)
            echo "Usage: $0 [options]"
            echo "  -m, --mapping FILE       Metadata TSV file (default: metadata/mapping_filtered.tsv)"
            echo "  -o, --out-dir DIR        Output directory for MACS3 results (default: analysis/PeakCalling_MACS3)"
            echo "  -p, --pooled-dir DIR     Directory of pooled BAMs (default: analysis/Pooled_BAMs)"
            echo "  -s, --pseudo-dir DIR     Directory of pseudoreplicate BAMs (default: analysis/Pseudoreplicates)"
            echo "  -y, --yaml FILE          YAML file with list of broad-peak targets (default: templates/target_peak_rules.yaml)"
            echo "      --force-broad        Force broad-peak calling for all targets"
            echo "  -d, --dry-run            Print commands without executing (simulation mode)"
            echo "  -h, --help               Show this help message"
            exit 0
            ;;
        *)
            log "ERROR" "🤷‍♀️ Unknown option: $1"
            exit 1
            ;;
    esac
done

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

# === Genome size determination ===
case "$REF_PREFIX" in
    hg38 | GRCh38) GENOME_SIZE="hs" ;;
    mm10 | GRCm38) GENOME_SIZE="mm" ;;
    dm6) GENOME_SIZE="dm" ;;
    ce11 | WBcel235) GENOME_SIZE="ce" ;;
    *)
        log "ERROR" "🤷‍♀️ Unknown genome 🧬 prefix: $REF_PREFIX"
        exit 1
        ;;
esac

# === MACS3 availability check ===
if ! command -v macs3 &> /dev/null; then
    log "ERROR" "❌ MACS3 😢 not found in PATH"
    exit 1
fi

REQUIRED_VERSION="3.0.3"
current_version=$(macs3 --version | awk '{print $2}')
if [[ "$current_version" != "$REQUIRED_VERSION" ]]; then
    log "WARN" "MACS3 version mismatch: expected $REQUIRED_VERSION, got $current_version"
fi

# === Samtools and bedtools availability check ===
if ! command -v samtools &> /dev/null || ! command -v bedtools &> /dev/null; then
    log "ERROR" "❌ samtools or bedtools not found — required for FRiP calculation"
    exit 1
fi

# === Load broad peak targets ===
if [[ -f "$TARGET_YAML" ]]; then
    BROAD_TARGETS=$(yq e '.broad_peak_targets[]' "$TARGET_YAML")
else
    log "ERROR" "❌ YAML 😢 not found: $TARGET_YAML"
    exit 1
fi

# === Parse metadata: Group samples by Condition+Target+Instrument ===
declare -A group_count
declare -A group_cond group_target group_instrument

# Find column indices
for i in "${!REQUIRED_COLS[@]}"; do
    case "${REQUIRED_COLS[$i]}" in
        "Sample_Type") idxSampleType=$i ;;
        "Condition") idxCondition=$i ;;
        "Instrument") idxInstrument=$i ;;
        "Target") idxTarget=$i ;;
    esac
done

while IFS=$'\t' read -r line || [[ -n "$line" ]]; do
    [[ "$line" =~ ^Sample_ID ]] && continue # skip header
    [[ -z "$line" ]] && continue
    IFS=$'\t' read -ra FIELDS <<< "$line"

    Sample_Type="${FIELDS[$idxSampleType]}"
    if [[ ! "$Sample_Type" =~ ^(IP|IP_rep[0-9]+|ChIP)$ ]]; then
        continue
    fi

    cond="${FIELDS[$idxCondition]}"
    target="${FIELDS[$idxTarget]}"
    inst="${FIELDS[$idxInstrument]}"

    group="${cond}_${target}_${inst}"
    group_count["$group"]=$((${group_count["$group"]:-0} + 1))
    group_cond["$group"]="$cond"
    group_target["$group"]="$target"
    group_instrument["$group"]="$inst"
done < "$MAPPING"

# === Process groups ===
SECONDS=0
record_metrics "PROCESS START"

# === Count total MACS3 calls for progress bar ===
total_calls=0
for group in "${!group_count[@]}"; do
    count=${group_count["$group"]}
    [[ $count -lt 2 ]] && continue
    total_calls=$((total_calls + 3)) # one for pooled, pseudo1, pseudo2
done
# === Check for zero work
if (( total_calls == 0 )); then
    log "WARN" "⚠️ No valid groups to process."
    exit 0
fi
current=0 # progress counter

for group in "${!group_count[@]}"; do
    log "INFO" "🔍 Found group: $group with ${group_count[$group]} IP samples"
    count=${group_count["$group"]}
    [[ $count -lt 2 ]] && continue

    cond="${group_cond[$group]}"
    target="${group_target[$group]}"
    inst="${group_instrument[$group]}"
    log "INFO" "🧠 Processing group: $group (n=$count samples)"

    # === Find matching control ===
    input_line=""
    while IFS=$'\t' read -r control_line; do
        [[ -z "$control_line" ]] && continue
        IFS=$'\t' read -ra control_fields <<< "$control_line"
        declare -A control_meta
        for i in "${!REQUIRED_COLS[@]}"; do
            control_meta["${REQUIRED_COLS[$i]}"]="${control_fields[$i]:-}"
        done

        if [[ "${control_meta[Sample_Type]}" =~ ^(Input|IgG|Mock)$ ]] &&
            [[ "${control_meta[Condition]}" == "$cond" ]] &&
            [[ "${control_meta[Instrument]}" == "$inst" ]] &&
            [[ "${control_meta[Target]}" == "$target" ]]; then
            input_line="$control_line"
            break
        fi
    done < <(tail -n +2 "$MAPPING")

    # === Process control ===
    has_input=false
    if [[ -n "$input_line" ]]; then
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
            IFS=_
            echo "${parts_input[*]}"
        )
        input_path="${BAM_DIR}/${Input_base}.bam"
        [[ -f "$input_path" ]] && has_input=true
    fi

    # === Determine peak type ===
    USE_BROAD=""
    if echo "$BROAD_TARGETS" | grep -qw "$target"; then
        USE_BROAD="--broad"
        REASON="Target=$target in broad list"
    else
        REASON="Target=$target not in broad list"
    fi
    $FORCE_BROAD && {
        USE_BROAD="--broad"
        REASON="User override"
    }

    # === Process pooled and pseudoreplicate BAMs ===
    for type in pooled pseudo1 pseudo2; do
        if [[ "$type" == "pooled" ]]; then
            bam_path="$POOLED_DIR/pooled_${group}.bam"
            base="pooled_${group}"
        else
            bam_path="$PSEUDO_DIR/${group}/${group}_${type}.bam"
            base="${group}_${type}"
        fi

        if [[ ! -f "$bam_path" ]]; then
            log "WARN" "⚠️ Missing BAM: $bam_path — ⏭️ skipping"
            continue
        fi

        # === Build output name ===
        if $has_input; then
            name="${base}_vs_${Input_base}"
        else
            name="${base}_noInput"
        fi
        macs_out_dir="$PEAK_DIR/$name"

        # === Dry Run ===
        if $DRY_RUN; then
            log "INFO" "🧯 [DRY-RUN] Would call MACS3: $name ($USE_BROAD)"
            continue
        fi

        # === Actual Processing ===
        mkdir -p "$macs_out_dir"

        # === Skip if output already exists (either broad or narrow peak file)

        expected_xls_file="$macs_out_dir/${name}_peaks.xls"
        expected_peak_file="$macs_out_dir/${name}_peaks.narrowPeak"
        expected_summit_file="$macs_out_dir/${name}_summits.bed"

        if [[ "$USE_BROAD" == "--broad" ]]; then
            expected_peak_file="$macs_out_dir/${name}_peaks.broadPeak"
            expected_summit_file=""
        fi

        # Skip if ZERO_PEAKS_MARKER was previously set
        if [[ -f "$macs_out_dir/ZERO_PEAKS_MARKER" ]]; then
            log "INFO" "⚠️ Skipping $name: previously produced zero peaks"
            continue
        fi

        # Optional debug print
        log "INFO" "🔍 Checking: $expected_peak_file, $expected_xls_file${expected_summit_file:+, $expected_summit_file}"

        # === Check if all expected outputs exist ===
        skip_call=false
        [[ -s "$expected_peak_file" && -s "$expected_xls_file" ]] && skip_call=true
        if [[ -n "$expected_summit_file" && ! -s "$expected_summit_file" ]]; then
            skip_call=false
        fi

        if $skip_call; then
            log "INFO" "✅ Skipping $name: All expected peak outputs exist"
            continue
        fi

        # === Run MACS3 ===
        if $has_input; then
            log "INFO" "📈 Calling peaks with control: $name"
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
                $USE_BROAD >> "$LOG_FILE" 2>&1

        else
            log "WARN" "🧩 No control found - calling ⛰️ without input: $name"
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
                $USE_BROAD >> "$LOG_FILE" 2>&1
        fi
        # === Check for zero or missing peak file ===
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

        # === Log progress ===
        log "INFO" "✅ Finished MACS3 call for $name"
        draw_progress "$name"
        printf "\n"

    done

done

# === Post-processing ===
record_metrics "🧾 MACS3 call complete.📊 Starting summary..."

# === Generate summary ===
print_header "📑 Generating MACS3 peak summary table..."

SUMMARY_OUT="$PEAK_DIR/macs3_pooled_summary.tsv"
echo -e "Sample\tType\tTotal_Peaks\tAvg_Peak_Length\tOutput_Dir\tWarning" > "$SUMMARY_OUT"

log "INFO" "📄 Initializing fresh macs3_summary.tsv"
for xls in "$PEAK_DIR"/*/*_peaks.xls; do
    [[ ! -f "$xls" ]] && continue
    sample=$(basename "$xls" _peaks.xls)
    output_dir=$(dirname "$xls")

    # === Determine sample type (pooled, pseudo, etc.) ===
    if [[ "$sample" == *pseudo1* ]]; then
        type="Pseudorep1"
        log "INFO" "🧪 Processing summary for pseudoreplicate 1: $sample"
    elif [[ "$sample" == *pseudo2* ]]; then
        type="Pseudorep2"
        log "INFO" "🧪 Processing summary for pseudoreplicate 2: $sample"
    elif [[ "$sample" == *pooled_* && "$sample" != *pseudo* ]]; then
        type="Pooled"
        log "INFO" "🌊 Processing summary for true pooled sample: $sample"
    else
        type="Other"
        log "INFO" "📦 Processing summary for unclassified sample: $sample"
    fi

    total_peaks=$(grep -v "^#" "$xls" | wc -l)
    avg_len=$(awk '!/^#/ {sum += $3 - $2; n++} END {print (n > 0 ? sprintf("%.1f", sum/n) : "NA")}' "$xls")

    warning="OK"
    if [[ "$total_peaks" -lt 100 ]]; then
        warning="VERY_LOW_PEAKS"
    elif [[ "$total_peaks" -lt 500 ]]; then
        warning="LOW_PEAKS"
    fi

    echo -e "${sample}\t${type}\t${total_peaks}\t${avg_len}\t${output_dir}\t${warning}" >> "$SUMMARY_OUT"
done

# === Organize output ===
narrow_files=$(find "$PEAK_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_peaks.narrowPeak")
broad_files=$(find "$PEAK_DIR" -mindepth 2 -maxdepth 2 -type f -name "*_peaks.broadPeak")

if ! $DRY_RUN; then
    if [[ -n "$narrow_files" ]]; then
        mkdir -p "$PEAK_DIR/narrowPeak"
        while IFS= read -r file; do
            dest="$PEAK_DIR/narrowPeak/$(basename "$file")"
            [[ "$file" != "$dest" && -f "$file" ]] && cp "$file" "$dest"
        done <<< "$narrow_files"
    fi

    if [[ -n "$broad_files" ]]; then
        mkdir -p "$PEAK_DIR/broadPeak"
        while IFS= read -r file; do
            dest="$PEAK_DIR/broadPeak/$(basename "$file")"
            [[ "$file" != "$dest" && -f "$file" ]] && cp "$file" "$dest"
        done <<< "$broad_files"
    fi
fi

# === Compute FRiP scores for pooled and pseudoreps ===

log "INFO" "📊 Calculating FRiP scores for all peak files..."

FRIP_SUMMARY="$PEAK_DIR/frip_summary.tsv"
TMP_FRIP="$FRIP_SUMMARY.tmp"
echo -e "Sample\tType\tFRiP\tReads_in_Peaks\tTotal_Reads\tPeak_File\tConclusion" > "$TMP_FRIP"

for peak_file in "$PEAK_DIR"/*/*_peaks.narrowPeak "$PEAK_DIR"/*/*_peaks.broadPeak; do
    [[ ! -f "$peak_file" ]] && continue

    sample=$(basename "$peak_file" | sed 's/_peaks\..*//')
    bam_base=$(echo "$sample" | sed 's/_vs_.*//')
    
    # Remove "pooled_" from pseudoreplicate names only
    if [[ "$bam_base" == pooled_*_pseudo* ]]; then
        bam_base="${bam_base#pooled_}"
    fi

    if [[ "$bam_base" == *_pseudo1 ]]; then
        group="${bam_base%_pseudo1}"
        bam_file="$PSEUDO_DIR/${group}/${bam_base}.bam"
        type="Pseudorep1"
    elif [[ "$bam_base" == *_pseudo2 ]]; then
        group="${bam_base%_pseudo2}"
        bam_file="$PSEUDO_DIR/${group}/${bam_base}.bam"
        type="Pseudorep2"
    else
        bam_file="$POOLED_DIR/${bam_base}.bam"
        type="Pooled"
    fi

    if [[ ! -f "$bam_file" ]]; then
        log "WARN" "⚠️ BAM not found for FRiP: $sample → looked for: $bam_file"
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

    echo -e "${sample}\t${type}\t${frip_score}\t${reads_in_peaks}\t${total_reads}\t${peak_file}\t${conclusion}" >> "$TMP_FRIP"
done

mv "$TMP_FRIP" "$FRIP_SUMMARY"

log "INFO" "📈 Pooled/Pseudorep FRiP scores written to: $FRIP_SUMMARY"
log "INFO" "📂 Using BAM for FRiP: $bam_file"

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
log "INFO" "🧾 Final organized peak files: $converted"
log "INFO" "🏁 Peak calling completed at $(date) (Elapsed: $elapsed)"
print_header "✅ [MODULE: ${SCRIPT_BASE_NAME}] Finished at $(date '+%F_%H-%M-%S')"

