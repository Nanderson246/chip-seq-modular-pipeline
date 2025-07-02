#!/bin/bash
# Module: 01_replicate_qc.sh
# Author: Nancy Anderson
# Description: Replicate-level QC with DeepTools correlation analysis and resource monitoring
################################################################################
# SOFTWARE REQUIREMENTS:
#
# This script requires the following tools and packages to be installed:
#
# Command-line tools:
#   - bash              (>= 4.0)
#   - coreutils         (for basic commands like `mv`, `sort`, `awk`, `cut`)
#   - GNU parallel      (for concurrent processing of BAMs)
#   - bedtools          (for PBC calculation using `bamtobed`)
#   - deeptools         (>= 3.5) [provides multiBamSummary, plotCorrelation, plotPCA, plotFingerprint, bamCoverage]
#   - Rscript           (>= 3.5) [for summary report and filtering]
#
# R Packages (auto-installed if missing):
#   - base R (data.frame, read.table, write.table)
#
# Optional but recommended:
#   - multiqc           (to aggregate QC results visually)
#
# Install deeptools (e.g. via conda):
#   conda install -c bioconda deeptools bedtools parallel
#
# Run this script inside a containerized environment (e.g. Docker/Singularity)
# to ensure reproducibility and avoid system-specific issues.
################################################################################

################################################################################
# USAGE:
#   Normal run: bash modules/pipeline2/01_replicate_qc.sh [options]
#   Dry run:    bash modules/pipeline2/01_replicate_qc.sh --dry-run [options]
#
# OPTIONS:
#   -d, --dry-run       Simulate execution without changes
#   -g, --group-fields  Fields to group replicates by (default: Instrument Condition Sample_Type Target)
#   -h, --help          Show this help message
#USAGE STANDALONE
#bash 01_replicate_qc.sh \
#  --meta /my/data/mapping.tsv \
#  --bam-dir /my/results/renamed \
#  --out-dir /my/output/qc \
#  --log-dir /my/output/logs \
#  --threads 8
# EXAMPLES:
#   # Dry-run with custom grouping
#   bash modules/pipeline202_replicate_qc.sh --dry-run --group-fields Condition Sample_Type
################################################################################
set -uo pipefail

# === SCRIPT ===
readonly VERSION="3.0.0"
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
GROUP_FIELDS=()
META=""
BAM_DIR=""
QC_OUT=""
FAIL_DIR=""
FILTERED_MAP=""
LOG_DIR=""
THREADS=4
FORCE_PBC=false 

# === CLI Parser ===
HELP_MSG=$'\nUsage: 01_replicate_qc.sh [options]
Options:
  -d, --dry-run              Simulate execution
  -g, --group-fields         Grouping fields (space separated)
  --meta FILE                Metadata file (default: metadata/mapping.tsv)
  --bam-dir DIR              Input BAM directory (default: analysis/Renamed_Cleaned)
  --out-dir DIR              Output directory (default: analysis/Replicate_QC)
  --log-dir DIR              Log directory (default: logs/pipeline2_analysis)
  --threads N                Number of threads for DeepTools (default: 4)
  --correlation-method STR   Correlation method for plotCorrelation [default: pearson]Options: pearson, spearman
  --force-pbc               FORCE_PBC, in case PBC need to be repeat[default: false]
                             
  -h, --help                 Show help'

while [[ $# -gt 0 ]]; do
    case "$1" in
        -d|--dry-run) DRY_RUN=true; shift ;;
        -g|--group-fields)
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                GROUP_FIELDS+=("$1")
                shift
            done
            ;;
        --meta)
          if [[ -z "${2:-}" || "$2" =~ ^- ]]; then
             log "ERROR" "Missing value for --meta"
             exit 1
          fi
          META="$2"; shift 2 ;;
        --bam-dir) BAM_DIR="$2"; shift 2 ;;
        --out-dir) QC_OUT="$2"; shift 2 ;;
        --log-dir) LOG_DIR="$2"; shift 2 ;;
        --threads)
            if [[ -z "${2:-}" || "$2" =~ ^- ]]; then
                log "ERROR" "Missing value for --threads"
                exit 1
            fi
            THREADS="$2"  
            shift 2
            ;;
        --correlation-method)
            CORRELATION_METHOD="$2"
            shift 2
            ;; 
        --force-pbc) 
            FORCE_PBC=true
            shift 
            ;;    
        -h|--help) echo "$HELP_MSG"; exit 0 ;;
        -*)
            echo "❌ Error: Unknown option $1"; echo "$HELP_MSG"; exit 1 ;;
        *) echo "❌ Error: Positional arguments not allowed"; echo "$HELP_MSG"; exit 1 ;;
    esac
done

# === Default grouping if not specified ===
[[ ${#GROUP_FIELDS[@]} -eq 0 ]] && GROUP_FIELDS=(Instrument Condition Sample_Type Target)

# === Initialization ===
BAM_DIR="${BAM_DIR:-$PROJECT_ROOT/analysis/Renamed_Cleaned}"
META="${META:-$PROJECT_ROOT/metadata/mapping.tsv}"
QC_OUT="${QC_OUT:-$PROJECT_ROOT/analysis/Replicate_QC}"
FAIL_DIR="${FAIL_DIR:-$PROJECT_ROOT/analysis/BAM_replicate_fail}"
FILTERED_MAP="${FILTERED_MAP:-$PROJECT_ROOT/metadata/mapping_filtered.tsv}"
CORRELATION_METHOD="${CORRELATION_METHOD:-pearson}"
LOG_DIR="${LOG_DIR:-$PROJECT_ROOT/logs/pipeline2}"
MODULE_LOG="${LOG_DIR}/01_replicate_qc_$(date +%F_%H-%M-%S).log"
PERF_LOG="${LOG_DIR}/01_replicate_qcperformance_metrics.log"

mkdir -p "$QC_OUT/deeptools" "$QC_OUT/tmp_groups" "$FAIL_DIR" "$LOG_DIR"


###############################################################################
# ⬇️  ADD THIS BLOCK – restores any BAMs you previously moved to $FAIL_DIR
###############################################################################
log "INFO" "♻️  Restoring any previously failed BAMs (if present)..."
restored=0
shopt -s nullglob          # silence the loop when no matches
for bam in "$FAIL_DIR"/*.bam; do
    base=$(basename "$bam")
    log "INFO" "🧬 Restoring $base → $BAM_DIR/"
    mv "$bam" "$BAM_DIR/"
    [[ -f "${FAIL_DIR}/${base}.bai" ]] && mv "${FAIL_DIR}/${base}.bai" "$BAM_DIR/"
    ((restored++))
done
shopt -u nullglob
if (( restored > 0 )); then
    log "INFO" "✅ Restored $restored BAM(s) from $FAIL_DIR"
else
    log "INFO" "📂 No BAMs to restore"
fi

# ===Correlation validation

if [[ ! "$CORRELATION_METHOD" =~ ^(pearson|spearman)$ ]]; then
    echo "❌ Error: Invalid correlation method: '$CORRELATION_METHOD'"
    echo "Allowed values: pearson, spearman"
    exit 1
fi

log "INFO" "📊 Correlation method: $CORRELATION_METHOD"

# === Resource Monitoring Functions ===

# === Resource Monitoring Note ===
# This module uses a simplified `record_metrics` function instead of the full
# PROCESS_METRICS/SYSTEM_METRICS setup from pipeline1
# Rationale:
#   - This QC module is not CPU/IO intensive.
#   - The simplified function captures enough performance context (CPU, memory, disk, output count).
#   - Keeps the script clean and more portable for standalone use.
# For alignment-heavy modules or performance benchmarking, see pipeline1/module6.


record_metrics() {
    local message="$1"
    local timestamp=$(date +%s)
    local cpu_usage=$(ps -p $$ -o %cpu | tail -n 1 | awk '{print $1}')
    local mem_usage=$(ps -p $$ -o %mem | tail -n 1 | awk '{print $1}')
    local disk_usage=$(df -h "${QC_OUT}" | tail -n 1)
    local group_count=$(ls -1 "${QC_OUT}/deeptools"/*.npz 2>/dev/null | wc -l)
    
    echo "[${timestamp}] ${message}" >> "$PERF_LOG"
    echo "  CPU: ${cpu_usage}% | Memory: ${mem_usage}% | Groups Processed: ${group_count}" >> "$PERF_LOG"
    echo "  Disk: ${disk_usage}" >> "$PERF_LOG"
}

calculate_pbc() {
    local bam_file="$1"
    local out_file="$2"

    bedtools bamtobed -i "$bam_file" | \
    awk '{print $1"\t"$2"\t"$3"\t"$6}' | \
    sort -k1,1 -k2,2n -k3,3n | \
    uniq -c | \
    awk -v bam="$bam_file" '
        {
            total++;
            if ($1 == 1) n1++;
        }
        END {
            if (total > 0) {
                pbc = n1 / total;
                if (pbc >= 0.9) {
                    cutoff = "≥ 0.9";
                    interp = "✅ High complexity";
                } else if (pbc >= 0.7) {
                    cutoff = "0.7–0.89";
                    interp = "⚠️ Moderate complexity";
                } else if (pbc >= 0.5) {
                    cutoff = "0.5–0.69";
                    interp = "⚠️ Low complexity";
                } else {
                    cutoff = "< 0.5";
                    interp = "❌ Very low complexity";
                }
                printf "%s\t%.4f\t%d\t%d\t%s\t%s\n", bam, pbc, n1, total, cutoff, interp;
            } else {
                printf "%s\tNA\t0\t0\tNA\tInsufficient data\n", bam;
            }
        }' >> "$out_file"
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
    echo -e "\n🛑 DRY-RUN MODE: No files will be modified"
    set -x  # Enable command echo
fi

start_timer

{
echo "=============================================="
echo "🧬 MODULE: ${SCRIPT_BASE_NAME}"
echo "📌 Purpose: Replicate QC with correlation analysis"
echo "💡 Dry-run mode: $DRY_RUN"
echo "🕒 Start time: $(date '+%F %T')"
echo "=============================================="

record_metrics "PROCESS START"

# === Validate Inputs ===
if [[ ! -f "$META" ]]; then
    log "ERROR" "❌ Error: Metadata file not found: $META"
    record_metrics "PROCESS FAILED - MISSING METADATA"
    exit 1
fi

if [[ ! -d "$BAM_DIR" ]]; then
    log "ERROR" "❌ Error: BAM directory not found: $BAM_DIR"
    record_metrics "PROCESS FAILED - MISSING BAM DIRECTORY"
    exit 1
fi

# === Dry-run Simulation ===
if [ "$DRY_RUN" = true ]; then
    log "INFO" "🔍 SIMULATION RESULTS:"
    record_metrics "DRY-RUN START"
    
    total_bams=$(ls "$BAM_DIR"/*.bam 2>/dev/null | wc -l)
    log "INFO" "Grouping fields that would be used: ${GROUP_FIELDS[*]}"
    log "INFO" "BAM files that would be processed (${total_bams} total):"
    ls "$BAM_DIR"/*.bam | head -3
    [ "$total_bams" -gt 3 ] && echo "... (and $((total_bams - 3)) more)"
    
    log "INFO" "\n⚡ Example commands that would run:"
    log "INFO" "multiBamSummary bins -b sample1.bam sample2.bam -out summary.npz"
    log "INFO" "plotCorrelation -in summary.npz -c pearson -p heatmap -o heatmap.png"
    log "INFO" "plotPCA -in summary.npz -o pca.png"
    log "INFO" "Rscript - (correlation summary)"
    
    log "INFO" "\n📊 Would generate outputs in:"
    log "INFO" "  QC Results: ${QC_OUT}/deeptools/"
    log "INFO" "  Failed BAMs: ${FAIL_DIR}/"
    log "INFO" "  Filtered Metadata: ${FILTERED_MAP}"
    
    record_metrics "DRY-RUN COMPLETE: Would process ${total_bams} BAMs"
    log "INFO" "\n✅ [DRY-RUN] Simulation completed (no changes made)"
    exit 0
fi

# === Main Processing ========================================================
log "INFO" "🔧 Processing metadata and grouping BAMs..."
record_metrics "GROUPING START"

# ─── 1) Clean stale group-list artefacts ─────────────────────────────────────
log "INFO" "🧹 Cleaning up old group list files..."
rm -f "$QC_OUT/tmp_groups/"*.list 2>/dev/null || true
log "INFO" "🧹 Cleared old .list files."

# ─── 2) Build column-index map from header  (NO subshell) ────────────────────
IFS=$'\t' read -r -a COLS <"$META"
declare -A COL_IDX
for i in "${!COLS[@]}"; do
    COL_IDX["${COLS[$i]}"]=$i
done

total_samples=0
grouped_samples=0

# ─── 3) Stream the data rows -------------------------------------------------
while IFS=$'\t' read -r -a FIELDS; do
    [[ ${#FIELDS[@]} -eq 0 ]] && continue   # skip blank / trailing lines

    # Row ➜ associative array “meta”
    declare -A meta
    for k in "${!COL_IDX[@]}"; do
        meta["$k"]="${FIELDS[${COL_IDX[$k]}]}"
    done

    log "INFO" "⚡DEBUG: Sample_ID=${meta[Sample_ID]}  Sample_Type=${meta[Sample_Type]}" >&2

    # Skip controls
    [[ "${meta[Sample_Type]}" =~ ^(Input|IgG|Mock)$ ]] && continue
    ((total_samples++))

    # ── build group key (strip _repX ONLY for key) ───────────────────────────
    key_parts=()
    for field in "${GROUP_FIELDS[@]}"; do
        val="${meta[$field]}"
        [[ "$field" == "Sample_Type" ]] && val="$(echo "$val" | sed -E 's/_rep[0-9]+$//')"
        [[ -n "$val" ]] && key_parts+=( "$(echo "$val" | tr -s ' ' '_' | tr -d '/')" )
    done
    group_key="$(IFS=_; echo "${key_parts[*]}")"

    # ── build expected BAM filename (Target included, nothing hard-coded) ────
    name_parts=(
        "${meta[Sample_ID]}"
        "${meta[Instrument]}"
        "${meta[Condition]}"
        "${meta[Replicate]}"
        "${meta[Sample_Type]}"    # full value incl. _repX
        "${meta[Target]}"         # whatever target is in the row
    )
    bam_name="$(IFS=_; echo "${name_parts[*]}").bam"
    bam_path="$BAM_DIR/$bam_name"

    if [[ -f "$bam_path" ]]; then
        echo "$bam_path" >> "$QC_OUT/tmp_groups/${group_key}.list"
        ((grouped_samples++))
    else
        log "WARN" "⚠️ BAM not found: $bam_name"
    fi
done < <(tail -n +2 "$META")   # feed rows without starting a subshell

# ─── 4) De-duplicate each list (safety) ──────────────────────────────────────
for list in "$QC_OUT/tmp_groups/"*.list; do
    [[ -f "$list" ]] && sort -u "$list" -o "$list"
done

# ─── 5) Summary & metrics ────────────────────────────────────────────────────
log "INFO" "📊 Sample grouping summary:"
log "INFO" "📦 Total samples processed: $total_samples"
log "INFO" "🧾 Successfully grouped:  $grouped_samples"
record_metrics "GROUPING COMPLETE: ${grouped_samples}/${total_samples} samples grouped"

# === PBC Calculation ===
PBC_FILE="$QC_OUT/pbc_metrics.tsv"

if [[ "$FORCE_PBC" != true && -s "$PBC_FILE" ]]; then          # -s  → file exists AND is non-empty
    log "INFO" "⏩ Skipping PBC calculation – $PBC_FILE already exists"
else
    log "INFO" "📊 Calculating PBC for all BAMs..."
    echo -e "Sample\tPBC\tN1\tNd\tCutoff\tInterpretation" > "$PBC_FILE"

    for bam in "$BAM_DIR"/*.bam; do
        [[ -f "$bam" ]] || continue
        calculate_pbc "$bam" "$PBC_FILE"
    done

    log "INFO" "✅ PBC calculation complete: $PBC_FILE"
fi

###############################################################################
# === BigWig creation (bamCoverage) ==========================================
###############################################################################
log "INFO" "📈 Generating CPM-normalised BigWigs with bamCoverage …"
record_metrics "BIGWIG START"

BIGWIG_DIR="${QC_OUT}/bigwig"
mkdir -p "$BIGWIG_DIR"

generate_bw() {
    local bam="$1"
    local bw_out="$BIGWIG_DIR/$(basename "${bam%.bam}.bw")"

    if [[ -s "$bw_out" ]]; then
        log "INFO" "⏩  Skipping existing BigWig: $(basename "$bw_out")"
        return 0
    fi

    bamCoverage \
        -b "$bam" \
        -o "$bw_out" \
        --normalizeUsing CPM \
        --binSize 25 \
        --extendReads 200 \
        --effectiveGenomeSize 2913022398 \
        --numberOfProcessors 4
}

export -f generate_bw
export BIGWIG_DIR

if [[ "$DRY_RUN" = true ]]; then
    log "INFO" "[DRY-RUN] Would run bamCoverage on $(ls "$BAM_DIR"/*.bam | wc -l) BAMs"
else
    parallel --jobs "$THREADS" generate_bw ::: "$BAM_DIR"/*.bam 2>> "$MODULE_LOG"
    log "INFO" "✅ BigWig files written to $BIGWIG_DIR"
fi

record_metrics "BIGWIG COMPLETE"


###############################################################################
# === ChIP-seq Fingerprint (1 line per replicate-group) ===
###############################################################################
log "INFO" "🖋️  Running plotFingerprint (group-level curves)…"

FINGERPRINTGroup_PNG="$QC_OUT/deeptools/Group_preIDR_fingerprint.png"
FINGERPRINT_METRICSGroup="$QC_OUT/deeptools/Group_preIDR_fingerprint.tsv"

# ── pick ONE representative BAM per group ────────────────────────────────────
declare -a FP_BAMS
declare -a FP_LABELS

for list in "$QC_OUT/tmp_groups/"*.list; do
    [[ -f "$list" ]] || continue
    grp=$(basename "$list" .list)          # e.g. NovaSeq_WT_H3K27ac
    first_bam=$(head -n 1 "$list")         # take the first replicate
    [[ -n "$first_bam" ]] && FP_BAMS+=( "$first_bam" ) && FP_LABELS+=( "$grp" )
done

if [[ ${#FP_BAMS[@]} -lt 2 ]]; then
    log "WARN" "⚠️  Not enough groups for plotFingerprint (${#FP_BAMS[@]} found) – skipping."
else
    if $DRY_RUN; then
        log "INFO" "[DRY-RUN] Would run: plotFingerprint -b ${FP_BAMS[*]} --labels ${FP_LABELS[*]} …"
    else
        plotFingerprint \
          -b "${FP_BAMS[@]}" \
          --labels "${FP_LABELS[@]}" \
          --minMappingQuality 30 \
          --skipZeros \
          --numberOfSamples 250000 \
          --extendReads 200 \
          --plotTitle "ChIP-seq fingerprint before IDR" \
          --plotFile "$FINGERPRINTGroup_PNG" \
          --outQualityMetrics "$FINGERPRINT_METRICSGroup" \
          --numberOfProcessors "$THREADS" \
          2>> "$MODULE_LOG"
        log "INFO" "✅ Fingerprint plot saved to $FINGERPRINTGroup_PNG"
    fi
fi

################################################################################
# === ChIP-seq Fingerprint (global across all cleaned BAMs) ===
################################################################################
log "INFO" "🖋️ Running plotFingerprint on cleaned replicates + inputs ..."

FINGERPRINT_PNG="$QC_OUT/deeptools/preIDR_fingerprint.png"
FINGERPRINT_METRICS="$QC_OUT/deeptools/preIDR_fingerprint.tsv"

# Collect every filtered BAM in the working directory
# - Adjust the glob pattern if your BAM names differ
mapfile -t ALL_BAMS < <(ls "$BAM_DIR"/*.bam 2>/dev/null)

# Convert filenames to simple labels (basename without .bam)
declare -a LABELS
for bam in "${ALL_BAMS[@]}"; do
    LABELS+=( "$(basename "$bam" .bam | cut -d'_' -f1,4)" )
done

if [[ ${#ALL_BAMS[@]} -lt 2 ]]; then
    log "WARN" "⚠️ Not enough BAM files for plotFingerprint (${#ALL_BAMS[@]} found) – skipping."
else
    if $DRY_RUN; then
        log "INFO" "[DRY-RUN] Would run: plotFingerprint -b ${ALL_BAMS[*]} --labels ${LABELS[*]} ..."
    else
        plotFingerprint \
          -b "${ALL_BAMS[@]}" \
          --labels "${LABELS[@]}" \
          --minMappingQuality 30 \
          --skipZeros \
          --numberOfSamples 250000 \
          --extendReads 200 \
          --plotTitle "ChIP-seq fingerprint before IDR" \
          --plotFile "$FINGERPRINT_PNG" \
          --outQualityMetrics "$FINGERPRINT_METRICS" \
          --numberOfProcessors "$THREADS" \
          2>> "$MODULE_LOG"
        log "INFO" "✅ Fingerprint plot saved to $FINGERPRINT_PNG"
    fi
fi

# === Replicate Group Summary Report ===
log "INFO" "📊 Summary of replicate counts per group:"
awk -F'\t' -v OFS='\t' '
BEGIN {
    print "Instrument", "Condition", "Sample_Type", "Target", "Replicate_Count"
}
NR==1 {
    for (i=1; i<=NF; i++) {
        if ($i == "Instrument") inst_idx = i
        if ($i == "Condition") cond_idx = i
        if ($i == "Sample_Type") type_idx = i
        if ($i == "Target") targ_idx = i
    }
    next
}
{
    key = $inst_idx "_" $cond_idx "_" $type_idx "_" $targ_idx
    count[key]++
}
END {
    for (k in count) {
        split(k, arr, "_")
        printf "%-10s %-10s %-12s %-12s %d\n", arr[1], arr[2], arr[3], arr[4], count[k]
    }
}
' "$META" | tee -a "$MODULE_LOG"


# === DeepTools Analysis ===
log "INFO" "📊 Running DeepTools analyses..."
record_metrics "DEEPTOOLS ANALYSIS START"
# === Main Processing =========================================================
# -----------------------------------------------------------------------------
log "INFO" "🔧 Processing metadata and grouping BAMs..."
record_metrics "GROUPING START"

total_groups=0
processed_groups=0

for list in "$QC_OUT/tmp_groups/"*.list; do
    [[ ! -f "$list" ]] && continue
    ((total_groups++))
    key=$(basename "$list" .list)
    mapfile -t group_bams < "$list"

    if [[ "${#group_bams[@]}" -lt 2 ]]; then
        log "WARN" "⚠️ Skipping $key - fewer than 2 replicates"
        continue
    fi

    ((processed_groups++))
    log "INFO" "🔬 Processing group $processed_groups/$total_groups: $key"
    log "INFO" "🔢 Samples: ${#group_bams[@]} replicates"
    record_metrics "GROUP PROCESSING START: $key (${#group_bams[@]} samples)"

    # Define output files
    npz="$QC_OUT/deeptools/${key}_summary.npz"
    corr_png="$QC_OUT/deeptools/${key}_correlation_heatmap.png"
    corr_tsv="$QC_OUT/deeptools/${key}_correlation.tsv"
    pca_png="$QC_OUT/deeptools/${key}_PCA_plot.png"
    pca_tsv="$QC_OUT/deeptools/${key}_pca_coords.tsv"

   log "INFO" "🧮 1. Calculating multiBamSummary..."
    if ! multiBamSummary bins -b "${group_bams[@]}" -out "$npz" --numberOfProcessors "$THREADS"; then
        log "ERROR" "❌ Error: multiBamSummary failed for group $key"
        record_metrics "GROUP PROCESSING FAILED: multiBamSummary - $key"
        continue
    fi

    log "INFO" "🧠 2. Generating correlation heatmap..."
    if ! plotCorrelation -in "$npz" -c "${CORRELATION_METHOD}" -p heatmap -o "$corr_png" \
        --plotTitle "$key ${CORRELATION_METHOD^} correlation" \
        --outFileCorMatrix "$corr_tsv" \
        --colorMap "coolwarm" \
        --plotWidth 12 --plotHeight 12; then
        log "WARN" "⚠️ Warning: plotCorrelation failed for group $key"
    fi

    log "INFO" "⚙️ 3. Running PCA..."
    if ! plotPCA -in "$npz" -o "$pca_png" \
        --plotTitle "$key PCA" \
        --outFileNameData "$pca_tsv"; then
        log "WARN" "⚠️ Warning: plotPCA failed for group $key"
    fi

    record_metrics "GROUP PROCESSING COMPLETE: $key"
 done
log "INFO" "🔄 Group processing summary:"
log "INFO" " 🔍 Total groups identified: $total_groups"
log "INFO" " 🟢 Successfully processed: $processed_groups"
record_metrics "DEEPTOOLS ANALYSIS COMPLETE: ${processed_groups}/${total_groups} groups processed"

# === R Summary ===
echo -e "\n📈 Generating QC summary..."
record_metrics "R ANALYSIS START"

Rscript - <<'EOF'
tryCatch({
  qc_dir <- Sys.getenv("QC_DIR", unset = "analysis/Replicate_QC/deeptools")
  out_summary <- file.path(qc_dir, "summary_report.tsv")
  out_samples <- file.path(qc_dir, "sample_level_stats.tsv")
  low_corr_file <- file.path(qc_dir, "low_correlation_bams.txt")

  files <- list.files(qc_dir, pattern = "_correlation.tsv", full.names = TRUE)
  if (length(files) == 0) stop("No correlation.tsv files found in: ", qc_dir)

  summary_list <- list()
  sample_stats_list <- list()
  low_corr_samples <- character()

  for (f in files) {
    cat("Processing:", f, "\n")
    mat <- tryCatch({
      read.table(f, header = TRUE, row.names = 1, check.names = FALSE)
    }, error = function(e) {
      message("Error reading file: ", f, " - ", e$message)
      return(NULL)
    })
    if (is.null(mat)) next

    group <- sub("_correlation.tsv", "", basename(f))
    m <- as.matrix(mat)
    if (any(is.na(m))) {
      warning("Matrix for ", group, " contains NA values.")
      next
    }

    if (ncol(m) < 3) {
      message("Group ", group, " has fewer than 3 samples, skipping outlier detection")
      next
    }

    # Compute average correlation per replicate (excluding self)
    avg_corr <- sapply(1:nrow(m), function(i) mean(m[i, -i], na.rm = TRUE))
    names(avg_corr) <- rownames(m)

    # Mark lowest replicate if its avg corr is below threshold
    worst <- names(which.min(avg_corr))
    if (avg_corr[worst] < 0.8) {
      low_corr_samples <- c(low_corr_samples, worst)
    }

    # Summary statistics
    off_diag <- m[lower.tri(m)]
    summary_list[[group]] <- data.frame(
      Group = group,
      N_Samples = ncol(m),
      Corr_Mean = round(mean(off_diag), 3),
      Corr_Median = round(median(off_diag), 3),
      Corr_Min = round(min(off_diag), 3),
      Corr_Max = round(max(off_diag), 3),
      Corr_SD = round(sd(off_diag), 3),
      Low_Corr_Count = sum(off_diag < 0.8)
    )

    sample_stats_list[[group]] <- data.frame(
      Sample = rownames(m),
      Group = group,
      Mean_Correlation = round(avg_corr, 3),
      Min_Correlation = apply(m, 1, function(x) round(min(x[x != 1]), 3))
    )
  }

  df_summary <- do.call(rbind, summary_list)
  df_samples <- do.call(rbind, sample_stats_list)

  write.table(df_summary, out_summary, sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(df_samples, out_samples, sep = "\t", quote = FALSE, row.names = FALSE)

  low_corr_samples <- unique(gsub("'", "", low_corr_samples))
  writeLines(low_corr_samples, low_corr_file)

  cat("[SUCCESS] Written to:\n")
  cat("- ", out_summary, "\n")
  cat("- ", out_samples, "\n")
  cat("- ", low_corr_file, "\n")

}, error = function(e) {
  message("[R ERROR]: ", e$message)
  quit(status = 1)
})

EOF


R_EXIT=$?
if [ $R_EXIT -ne 0 ]; then
  log "ERROR" "❌ R analysis failed with exit code $R_EXIT"
  record_metrics "R ANALYSIS FAILED"
else
  record_metrics "R ANALYSIS COMPLETE"
fi

# === Handle Failed BAMs ===
log "INFO" "🤖 Checking for poorly correlated BAMs (r < 0.8)..."
record_metrics "FAILED BAM CHECK START"
failed_count=0

# Print out correlation matrix inspection per file
for corr_file in "$QC_OUT"/deeptools/*_correlation.tsv; do
    [[ ! -f "$corr_file" ]] && continue
    group=$(basename "$corr_file" _correlation.tsv)
    log "INFO" "🧠 Checking $group for low correlations..."

    IFS=$'\t' read -r -a header < <(grep -v '^#' "$corr_file" | head -n 1 | tr -d "'")

    tail -n +2 "$corr_file" | while IFS=$'\t' read -r -a row; do
        sample_name=$(echo "${row[0]}" | tr -d "'")
        correlations=("${row[@]:1}")  # skip sample name column

        # Filter out self-correlation
        filtered=()
        for corr in "${correlations[@]}"; do
            [[ "$corr" != "1.0000" && "$corr" =~ ^[0-9]+\.[0-9]+$ ]] && filtered+=("$corr")
        done

        if [[ "${#filtered[@]}" -eq 0 ]]; then
            log "WARN" "⚠️ No valid correlations found for $sample_name"
            continue
        fi

        # (Optional: You can also echo lowest correlation per sample if desired)
        min_corr=$(printf '%s\n' "${filtered[@]}" | sort -n | head -n 1)
        echo "🔎 $sample_name: lowest r = $min_corr"
    done
done

# Now move the worst replicate per group from the R output
LOW_CORR_LIST="$QC_OUT/deeptools/low_correlation_bams.txt"
if [[ -f "$LOW_CORR_LIST" ]]; then
    while read -r bam; do
        [[ -z "$bam" ]] && continue
        bam_clean=$(basename "$bam")

        if [[ -f "$BAM_DIR/$bam_clean" ]]; then
            echo "❌ Moving low-correlation BAM: $bam_clean"
            mkdir -p "$FAIL_DIR"
            mv "$BAM_DIR/$bam_clean" "$FAIL_DIR/"
            [[ -f "$BAM_DIR/${bam_clean}.bai" ]] && mv "$BAM_DIR/${bam_clean}.bai" "$FAIL_DIR/"
            ((failed_count++))
        else
            log "WARN" "⚠️ BAM not found: $BAM_DIR/$bam_clean"
        fi
    done < "$LOW_CORR_LIST"
else
    log "INFO" "✅ No low-correlation BAM list found - skipping"
fi

record_metrics "FAILED BAM CHECK COMPLETE: ${failed_count} BAMs flagged"

# === Update Metadata ===
log "INFO" "📝 Updating filtered metadata..."
record_metrics "METADATA UPDATE START"

if [ "$failed_count" -gt 0 ]; then
    log "INFO" "🗂️ Filtering metadata based on failed Sample_IDs..."
    
    # Extract Sample_IDs from BAM filenames (assumes Sample_ID is prefix before first "_")
    cut -d'_' -f1 <(find "$FAIL_DIR" -name "*.bam" -exec basename {} .bam \;) | sort -u > "$QC_OUT/tmp_failed_ids.txt"
    
#    grep -v -f <(find "$FAIL_DIR" -name "*.bam" -exec basename {} .bam \;) "$META" > "$FILTERED_MAP"

# Remove lines from metadata where Sample_ID matches failed ones
    awk 'NR==FNR{bad[$1]; next} !($1 in bad)' "$QC_OUT/tmp_failed_ids.txt" "$META" > "$FILTERED_MAP"
    log "INFO" "🗑️ Removed $failed_count samples from metadata based on Sample_ID prefix"
else
    cp "$META" "$FILTERED_MAP"
    log "INFO" "✅ No samples removed - keeping original metadata"
fi

record_metrics "METADATA UPDATE COMPLETE"

# === Cleanup ===
log "INFO" "🧹 Cleaning up failed BAMs..."
record_metrics "CLEANUP START"

if [ "$failed_count" -gt 0 ]; then
    find "$FAIL_DIR" -type f -name "*.bam" | while read -r failed_bam; do
        base=$(basename "$failed_bam")
        [[ -f "$BAM_DIR/$base" ]] && rm -v "$BAM_DIR/$base"
        [[ -f "$BAM_DIR/${base}.bai" ]] && rm -v "$BAM_DIR/${base}.bai"
    done
    log "INFO" "🗑️ Removed $failed_count BAMs from main directory"
else
    log "INFO" "✅ No BAMs to remove"
fi

record_metrics "CLEANUP COMPLETE"


# === Completion ===
log "INFO" "\n=============================================="
log "INFO" "✅ REPLICATE QC COMPLETED"
log "INFO" "📊 Samples processed: $total_samples"
log "INFO""📊 Groups analyzed: $processed_groups/$total_groups"
log "INFO" "⚠️  Failed BAMs identified: $failed_count"
log "INFO" "⏱️  Elapsed time: $(get_elapsed_time)"
log "INFO" "📂 Output directory: $QC_OUT/deeptools"
log "INFO" "📝 Filtered metadata: $FILTERED_MAP"
log "INFO" "📊 Performance metrics saved to: ${PERF_LOG}"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "=============================================="

record_metrics "PROCESS COMPLETE"
} | tee -a "$MODULE_LOG"
