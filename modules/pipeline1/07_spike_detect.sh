#!/usr/bin/env bash
# Module: 05_2_spike_detect.sh
# Author: Nancy Anderson
# Description: Detect exogenous spike

################################################################################
# SOFTWARE REQUIREMENTS:
#   • bash               >= 4.x
#   • seqtk              >= 1.3        # For read subsampling
#   • bwa                >= 0.7.x      # For read alignment to spike references
#   • samtools           >= 1.10       # For BAM/SAM read filtering and counting
#   • yq                 >= 4.x        # For YAML parsing (metadata and rules)
#   • gzip/zcat          (coreutils)   # For decompressing .fastq.gz files
#   • findutils          (find, xargs) # For file discovery
################################################################################

################################################################################
# Usage:
#default  bash modules/pipeline1/07_spike_detect.sh (hg38 default prefix)
#using mm10: bash modules/pipeline1/07_spike_detect.sh.sh --ref-prefix mm10
# bash modules/pipeline1/07_spike_detect.sh /path/to/fastq_dir /path/to/SpikeinReference target_peak_rules.yaml mapping_schema.yaml
#Standalone use
# bash modules/pipeline1/07_spike_detect.sh.sh \
# --fastq-dir results/Trimmed \
# --spike-dir SpikeinReference \
# --peak-rules templates/target_peak_rules.yaml \
# --meta-schema templates/mapping_schema.yaml \
# --meta-tsv metadata/mapping.tsv \
# --ref-prefix mm10
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
    local message="${2:-}"
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
FASTQ_DIR=""
SPIKE_REF_DIR=""
TARGET_RULES=""
META_FILE=""
META_TSV=""
REF_PREFIX="hg38"
READS=200000
LOG_DIR="${LOG_DIR:-logs/pipeline1}"
MODULE_LOG="${LOG_DIR}/${SCRIPT_BASE_NAME}_${TIMESTAMP}.log"
OUTPUT=${OUTPUT:-${PROJECT_ROOT}/results/spike_analysis}
mkdir -p "$LOG_DIR"
mkdir -p ${OUTPUT}

log "INFO" ""
log "INFO" "=============================================="
log "INFO" "🔍 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Purpose: Detecting exogenous spiked type"
log "INFO" "📁 Input: results/Trimmed"
log "INFO" "🗓️  Timestamp: $TIMESTAMP"
log "INFO" "📦 Script version: ${VERSION} (${SCRIPT_NAME})"
log "INFO" "=============================================="

HELP_MSG=$'\nUsage: detect_random_spikein.sh [options]\nOptions:\n  -d, --dry-run              Simulate execution\n  --fastq-dir                Input Fastqz.trimmed samples folder (default: results/Trimmed)\n  --spike-dir                Spike reference genomes directory (default: SpikeinReference)\n  --peak-rules               Path to peak rules YAML (default: templates/target_peak_rules.yaml)\n  --meta-schema              Path to metadata schema YAML (default: templates/mapping_schema.yaml)\n  --meta-tsv                 Metadata TSV file (default: metadata/mapping.tsv)\n  --ref-prefix               Genome prefix for host (default: hg38)\n  -h, --help                 Show help'

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dry-run) DRY_RUN=true; shift ;;
        --fastq-dir) FASTQ_DIR="$2"; shift 2 ;;
        --spike-dir) SPIKE_REF_DIR="$2"; shift 2 ;;
        --peak-rules) TARGET_RULES="$2"; shift 2 ;;
        --meta-schema) META_FILE="$2"; shift 2 ;;
        --meta-tsv) META_TSV="$2"; shift 2 ;;
        --ref-prefix) REF_PREFIX="$2"; shift 2 ;;
        -h|--help) echo "$HELP_MSG"; exit 0 ;;
        -*) echo "❌ Error: Unknown option $1"; echo "$HELP_MSG"; exit 1 ;;
        *) echo "❌ Error: Positional arguments not allowed"; echo "$HELP_MSG"; exit 1 ;;
    esac
done

FASTQ_DIR=${FASTQ_DIR:-${PROJECT_ROOT}/results/Trimmed}
SPIKE_REF_DIR=${SPIKE_REF_DIR:-${PROJECT_ROOT}/SpikeinReference}
TARGET_RULES=${TARGET_RULES:-${PROJECT_ROOT}/templates/target_peak_rules.yaml}
META_FILE=${META_FILE:-${PROJECT_ROOT}/templates/mapping_schema.yaml}
META_TSV=${META_TSV:-${PROJECT_ROOT}/metadata/mapping.tsv}

if [[ -f "$META_TSV" ]]; then
    if awk -F'\t' 'NR==1 {for (i=1; i<=NF; i++) if ($i=="Spike_Type") col=i} NR>1 && col && $col != "" {count++} END {exit count>0 ? 0 : 1}' "$META_TSV"; then
        log INFO "Spike_Type already defined in $META_TSV — skipping detection."
        exit 0
    fi
else
    log WARN "Metadata TSV not found at $META_TSV — continuing with detection."
fi


# === Tool Verification ===
for tool in seqtk bwa samtools yq; do
    command -v "$tool" >/dev/null 2>&1 || { 
        log "ERROR" "$tool is not installed. Aborting."
        exit 1
    }
done


# === Dry-run or real mode: Select 3 R1 samples and process ===

log INFO "🔁 Selecting 3 random paired samples for spike-in detection..."
R1_FILES=($(find "$FASTQ_DIR" -name "*_1_trimmed.fastq.gz" | shuf | head -n 3))

if [ "${#R1_FILES[@]}" -eq 0 ]; then
    log ERROR "❌ No *_1_trimmed.fastq.gz files found in $FASTQ_DIR"
    exit 1
fi

if [ "$DRY_RUN" = true ]; then
    echo -e "\n🛑 DRY-RUN MODE: No files will be modified"
    set -x
    log INFO "🔍 SIMULATION RESULTS:"
    log INFO "Total candidate R1 FASTQ files: ${#R1_FILES[@]}"
    for fq1 in "${R1_FILES[@]}"; do
        fq2="${fq1/_1/_2}"
        log INFO "Selected pair: $fq1 + $fq2"
    done
    log INFO "✅ [DRY-RUN] Simulation completed (no changes made)"
    exit 0
fi

# === Process 3 samples ===
SELECTED_SAMPLES=()
for fq1 in "${R1_FILES[@]}"; do
    fq2="${fq1/_1/_2}"
    sample_name=$(basename "$fq1" | sed -E 's/_1(_trimmed)?\.fastq\.gz$//')
    combined_fq="${OUTPUT}/${sample_name}_combined.fastq"
    subsampled_fq="${OUTPUT}/${sample_name}_subsample.fastq"

    log INFO "🎲 Selected pair: $fq1 + $fq2"
    log INFO "📎 Combining into: $combined_fq"

    zcat "$fq1" "$fq2" > "$combined_fq"
    
    log INFO "📉 Subsampling $READS reads to: $subsampled_fq"
    seqtk sample -s100 "$combined_fq" $READS > "$subsampled_fq"

    # Save for later usage
    SELECTED_SAMPLES+=("$subsampled_fq")
done

ALL_SPIKES=()

# === Detect spike-in for each subsampled FASTQ ===
for subsampled_fq in "${SELECTED_SAMPLES[@]}"; do
    base=$(basename "$subsampled_fq" .fastq)
    log INFO "\n🔬 Processing: $base"

    BEST_SPIKE="none"
    BEST_MAPPED=0
    SPIKE_COUNTS="${OUTPUT}/${base}_spike_counts.tsv"
    > "$SPIKE_COUNTS"

    for SPIKE_PATH in $(find "$SPIKE_REF_DIR" -type f -name "*.fa" ! -lname "*"); do
        SPIKE_NAME=$(basename "$SPIKE_PATH" .fa)
     
        if [[ "$SPIKE_NAME" == "$REF_PREFIX" ]]; then
            log INFO "⏭️  Skipping host genome $SPIKE_NAME as spike candidate"
            continue
        fi

        echo -n "Aligning to $SPIKE_NAME... "
        mapped=$(bwa mem "$SPIKE_PATH" "$subsampled_fq" 2>/dev/null | samtools view -F 4 -c -)
        echo "$mapped mapped reads"

        echo -e "$mapped\t$SPIKE_NAME" >> "$SPIKE_COUNTS"

        if (( mapped > BEST_MAPPED )); then
            BEST_MAPPED=$mapped
            BEST_SPIKE=$SPIKE_NAME
        fi
    done

    # Handle phiX special case
    if [[ "$BEST_SPIKE" == "phiX174" ]]; then
        SECOND_BEST=$(sort -nr "$SPIKE_COUNTS" | awk 'NR==2 {print $2}')
        SECOND_MAPPED=$(sort -nr "$SPIKE_COUNTS" | awk 'NR==2 {print $1}')
        if (( BEST_MAPPED < SECOND_MAPPED * 2 )); then
            log INFO "phiX174 likely sequencing control — ignoring."
            BEST_SPIKE="none"
            BEST_MAPPED=0
        fi
    fi

    echo "✅ Best spike-in: $BEST_SPIKE ($BEST_MAPPED reads)"
    
    ALL_SPIKES+=("$BEST_SPIKE")
    
    # Check if any target expects spike-in
    SCALING_NEEDED="unknown"
    WARNING=""

    if [[ -f "$TARGET_RULES" && -f "$META_FILE" ]]; then
        mapfile -t broad_targets < <(yq '.broad_peak_targets[]' "$TARGET_RULES")
        mapfile -t used_targets < <(yq '.samples[].Target' "$META_FILE" | sort -u)

        for target in "${used_targets[@]}"; do
            if printf '%s\n' "${broad_targets[@]}" | grep -q -x "$target"; then
                SCALING_NEEDED="true"
                if [[ "$BEST_SPIKE" == "none" ]]; then
                    WARNING="⚠️ Spike-in required for $target, but not detected"
                fi
                break
            fi
        done
        [[ "$SCALING_NEEDED" == "unknown" ]] && SCALING_NEEDED="false"
    fi

    # Save result
    SPIKE_TMP="${OUTPUT}/${base}.spike.tmp"
    echo -e "$subsampled_fq\t$BEST_SPIKE\t$BEST_MAPPED\t$SCALING_NEEDED\t$WARNING" > "$SPIKE_TMP"
    log INFO "📄 Result written: $SPIKE_TMP"
done

# === Update metadata TSV with detected Spike_Type per sample ===
# === Majority vote across subsampled results ===

declare -A SPIKE_COUNTS=()
MAJORITY_SPIKE="none"

for spike in "${ALL_SPIKES[@]}"; do
    ((SPIKE_COUNTS["$spike"]++))
done

for spike in "${!SPIKE_COUNTS[@]}"; do
    if (( SPIKE_COUNTS["$spike"] >= 2 )); then
        MAJORITY_SPIKE="$spike"
        break
    fi
done

log INFO "📊 Consensus spike-in across 3 samples: $MAJORITY_SPIKE"
if [[ "$MAJORITY_SPIKE" == "none" ]]; then
    log WARN "⚠️ No consensus spike-in found in the 3 samples — 'Spike_Type' will be 'none'"
fi


# === Update metadata.tsv with consensus ===
if [[ -f "$META_TSV" ]]; then
    TMP_META="tmp_metadata.tsv"
    HEADER=$(head -n 1 "$META_TSV")

    [[ ! -f "${META_TSV}.bak" ]] && cp "$META_TSV" "${META_TSV}.bak"

    if echo "$HEADER" | grep -q -w "Spike_Type"; then
        awk -v OFS='\t' -v spike="$MAJORITY_SPIKE" '
BEGIN { spike_col = 0 }
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "Spike_Type") spike_col = i;
    }
    print $0;
    next;
}
{
    if (spike_col > 0 && ($(spike_col) == "" || $(spike_col) == "none" || $(spike_col) == "NA")) {
        $(spike_col) = spike;
    }
    print $0;
}' "${META_TSV}.bak" > "$TMP_META"

    else
        awk -v OFS='\t' -v spike="$MAJORITY_SPIKE" '
        NR==1 { print $0, "Spike_Type"; next }
        { print $0, spike }' "${META_TSV}.bak" > "$TMP_META"
    fi

    mv "$TMP_META" "$META_TSV"
    log INFO "✅ mapping.tsv updated with consensus Spike_Type: $MAJORITY_SPIKE"
    UPDATED_COUNT=$(awk 'END{print NR-1}' "$META_TSV")
    log INFO "📄 Total metadata entries: $UPDATED_COUNT"
else
    log ERROR "❌ Metadata TSV not found at $META_TSV"
fi

# === Final summary rebuild ===
SUMMARY_OUT="${OUTPUT}/spike_summary.tsv"
echo -e "FASTQ\tDetected_Spike\tMapped_Reads\tScaling_Expected\tWarning" > "$SUMMARY_OUT"
cat "${OUTPUT}"/*.spike.tmp >> "$SUMMARY_OUT"
log INFO "📊 Final spike-in summary rebuilt: $SUMMARY_OUT"

# === Completion ===
echo ""
log "INFO" "=============================================="
log "INFO" "✅ MODULE COMPLETED"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "📂 Output location: $SUMMARY_OUT"
log "INFO" "=============================================="
