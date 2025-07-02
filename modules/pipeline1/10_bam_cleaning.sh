#!/usr/bin/env bash
# Module: 10_bam_cleaning.sh
# Author: Nancy Anderson
# Description: Clean BAMs with MAPQ, blacklist, and mitochondrial filtering
################################################################################
# SOFTWARE REQUIREMENTS
#
# This script requires the following tools installed and accessible in PATH,
# or wrapped under the 'tools/' folder for pipeline-based resolution:
#
# Required Tools:
#   • bash        - Unix shell interpreter (v4+)
#   • samtools    - SAM/BAM operations (≥ v1.9 recommended)
#   • awk         - For BAM filtering, header editing, etc.
#   • bedtools    - For blacklist region removal (≥ v2.29)
#   • qualimap    - BAM QC and summary statistics
#   • multiqc     - For aggregating QC reports
#   • Rscript     - For spike-in QC plot (used at end)
#
# Optional:
#   • iostat      - Disk I/O monitoring (from `sysstat` package)
#   • top         - System-wide CPU stats (used for performance metrics)
#   • free        - Memory usage stats (from `procps`)
#
# Dependencies:
#   - Java is needed for Qualimap (usually ≥ Java 8)
#   - R packages may be required for spike QC plots (defined in plot script)
#
# Installation (Debian/Ubuntu-based systems):
#   sudo apt install samtools bedtools default-jre procps sysstat gawk r-base-core
#   pip install multiqc
#
# Notes:
#   - Place `qualimap` binary in: tools/qualimap/qualimap and ensure it's executable
#   - `qualimap` must be able to find Java
#
#   Example for wrapper setup:
#     chmod +x tools/qualimap/qualimap
#
################################################################################

################################################################################
# USAGE:
#   Normal run: bash modules/pipeline1/10_bam_cleaning.sh [options] <reference_folder>
#   Dry run:    bash modules/pipeline_1/10_bam_cleaning.sh --dry-run <reference_folder>
#
# Standalone usage
#   bash 10_BAM_cleaning.sh --in <BAMS_folder> --out <cleaned_bams_folder>  /Reference/hg19
# OPTIONS:
#   -d, --dry-run    Simulate execution without making changes
#   -h, --help       Show this help message
#
# PARAMETERS:
#   reference_folder   Path to directory containing reference files
#
# EXAMPLES:
#   bash 10_bam_cleaning.sh -d /path/to/reference 
#   bash 10_bam_cleaning.sh --help
#   bash modules/pipeline1/10_bam_cleaning.sh --ref Reference --prefix mm10
# bash modules/pipeline1/10_bam_cleaning.sh  --prefix mm10
################################################################################

set -uo pipefail
: "${TOOLS_DIR:=tools}"


# === IO Paths (can be overridden in standalone use) ===
    INPUT_DIR="${INPUT_DIR:-results/Filtered/Deduplicated}"
    CLEANED_DIR="${CLEANED_DIR:-results/Filtered/Cleaned}"
    FILTERED_DIR="$(dirname "$CLEANED_DIR")"
    METRICS_DIR="${FILTERED_DIR}/Metrics"
    LOG_DIR="logs"
    OUTPUT_FORMATS=("human" "csv" "tsv")  # Supported report formats
    REPORT_FORMAT="human" 
    META_TSV="${META_TSV:-metadata/mapping.tsv}"
# === Log files ===
TIMESTAMP=$(date +%F_%H-%M-%S)
MODULE_LOG="${LOG_DIR}/09_bam_cleaning_${TIMESTAMP}.log"
PERF_LOG="${LOG_DIR}/09_bam_cleaning_performance_metrics.log"

mkdir -p "$LOG_DIR" "$CLEANED_DIR" "$METRICS_DIR" "${METRICS_DIR}/BAM_QC_MultiQC"


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


# === Load wrappers if not already loaded ===
echo "[INFO] 🔍 Checking for local tool wrappers..."

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

# === Qualimap wrapper ===
log "INFO" "🔧 Loading tool wrappers..."
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


# === CLI Parsing ===
parse_args() {
    DRY_RUN=false

    # === Default values (override with CLI) ===
    REF_FOLDER="${REF_FOLDER:-Reference}"
    REF_PREFIX="${REF_PREFIX:-hg38}"

    while [[ $# -gt 0 ]]; do
        case "$1" in
            --in)
                INPUT_DIR="$2"
                shift 2
                ;;
            --out)
                CLEANED_DIR="$2"
                shift 2
                ;;
            --ref)
                REF_FOLDER="$2"
                shift 2
                ;;
            --prefix)
                REF_PREFIX="$2"
                shift 2
                ;;
            --format)
                if [[ " ${OUTPUT_FORMATS[@]} " =~ " $2 " ]]; then
                    REPORT_FORMAT="$2"
                    shift 2
                else
                    echo "Invalid format: $2. Supported: ${OUTPUT_FORMATS[*]}"
                    exit 1
                fi
                ;;
            -d|--dry-run)
                DRY_RUN=true
                shift
                ;;
            -h|--help)
                echo "Usage: $SCRIPT_NAME [--ref DIR] [--prefix NAME] [--in DIR] [--out DIR] [--format FORMAT] [-d]"
                echo "Defaults: --ref=Reference --prefix=mm10"
                exit 0
                ;;
            *)
                log "ERROR" "❌ Unknown argument: $1"
                exit 1
                ;;
        esac
    done

    # === Construct paths ===
    REF_GENOME="${REF_FOLDER}/${REF_PREFIX}/${REF_PREFIX}.fa"
    BLACKLIST="${REF_FOLDER}/${REF_PREFIX}/${REF_PREFIX}-blacklist.bed"

    # === Validate ===
    if [[ ! -f "$REF_GENOME" ]]; then
        log "ERROR" "❌ Reference genome not found: $REF_GENOME"
        exit 1
    fi

    if [[ ! -f "$BLACKLIST" ]]; then
        log "WARN" "⚠️ Blacklist BED file missing: $BLACKLIST"
    fi
}


# === Functions ===

# === Resource Tracking Variables ===
declare -A PROCESS_METRICS=()
declare -A SYSTEM_METRICS=()


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
    PROCESS_METRICS["${stage}_time"]=$(awk -v s="$start_time" -v e="$end_time" 'BEGIN {print e - s}')
    PROCESS_METRICS["${stage}_cpu"]=$(ps -p $pid -o %cpu --no-headers | awk '{print $1}')
    PROCESS_METRICS["${stage}_mem"]=$(ps -p $pid -o rss --no-headers | awk '{printf "%.2f", $1/1024}')
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
            generate_human_report
            ;;
        csv)
            generate_csv_report >> "$PERF_LOG"
            ;;
        tsv)
            generate_tsv_report >> "$PERF_LOG"
            ;;
    esac
}
generate_performance_report() {
    echo -e "\n${BLUE}=== PERFORMANCE REPORT ===${NC}"
    printf "%-20s %10s %8s %10s\n" "Stage" "Time(s)" "CPU(%)" "Mem(MB)"
    echo "-----------------------------------------------"
    local total_time=0
    for stage in mapq_filter blacklist_filter final_clean indexing; do
        local time="${PROCESS_METRICS["${stage}_time"]:-0}"
        local cpu="${PROCESS_METRICS["${stage}_cpu"]:-0}"
        local mem="${PROCESS_METRICS["${stage}_mem"]:-0}"

        printf "%-20s %10.2f %8.1f %10.1f\n" "$stage" "$time" "$cpu" "$mem"
        total_time=$(awk -v total="$total_time" -v t="$time" 'BEGIN {print total + t}')
    done

    echo "-----------------------------------------------"
     printf "%-15s %10.2f\n" "TOTAL" "$total_time"
    echo -e "${BLUE}============================================${NC}"
}

record_metrics() {
    local message="$1"
    local timestamp=$(date +%s)
    echo "[${timestamp}] ${message}" >> "$PERF_LOG"
}

get_elapsed_time() {
    local duration=$SECONDS
    echo "$((duration / 60))m $((duration % 60))s"
}

main() {
    # === Main Processing ===
exec > >(tee -a "$MODULE_LOG") 2>&1    
SECONDS=0    
    
processed=0
skipped=0
total_samples=$(ls -1 "${INPUT_DIR}/"*.dedup.bam 2>/dev/null | wc -l)
    if [[ "$total_samples" -eq 0 ]]; then
    log "WARN" "⚠️ No BAMs found in $INPUT_DIR. Exiting."
    exit 0
    fi

log "INFO" "🔍 Found ${total_samples} BAM files to process"
record_metrics "INPUT DETECTED: ${total_samples} BAMs"

for sample in "${INPUT_DIR}/"*.dedup.bam; do
    [ -e "$sample" ] || continue
    base=$(basename "$sample" .dedup.bam)
    clean_bam="${CLEANED_DIR}/${base}_clean.bam"

if [[ -f "$clean_bam" && -f "${clean_bam}.bai" ]]; then
    log "INFO" "⏩ Skipping $base — clean BAM already exists"
    ((skipped++))
    continue
fi  
    
    if "$DRY_RUN"; then
    log "INFO" "🔍 [DRY-RUN] Would process sample: $base"

    echo "  🧹 Check/add MAPQ ≥30: samtools view -h -q 30 -b $sample > ${CLEANED_DIR}/${base}_filtered.bam"
    echo "  🧫  Check for spike-in type from: $META_TSV"

    echo "  🔬 If spike-in: remove spike reads and save to _host.bam, keep spike.bam separately"
    echo "  🚫 Remove blacklist regions using: bedtools intersect -v -abam ... -b $BLACKLIST"
    echo "  🔋 Separate chrM/MT (mitochondrial) reads from host genome"
    echo "  🚿 Save clean BAM to: ${CLEANED_DIR}/${base}_clean.bam"
    echo "  🔖 Index clean BAM: samtools index ${CLEANED_DIR}/${base}_clean.bam"
    echo "  🧠 Generate QC: samtools flagstat and qualimap for clean BAM"
    echo ""

    continue
fi

    
 # === Get spike type from metadata (insert after defining 'base') ===
spike_type=$(awk -F'\t' -v sample="$base" '
    NR==1 {
        for(i=1;i<=NF;i++) {
            if($i=="Sample_ID") sid_col=i
            if($i=="Spike_Type") spike_col=i
        }
    }
    $sid_col == sample {print $spike_col; exit}' "$META_TSV")
spike_type="${spike_type:-none}"


   
    
    log "INFO" "\n⚡ Processing sample: $base"
    record_metrics "SAMPLE START: $base"
    
    
#=====================================================================
# 1. MAPQ FILTERING (EXISTING)
#=====================================================================
 # a. Count initial unfiltered reads
reads_initial=$(samtools view -c "$sample")

    log "INFO" "1. 🧹 Applying MAPQ ≥30 filter..."
    start_monitoring "mapq_filter" $$
    if ! samtools view -h -q 30 -b "$sample" > "${CLEANED_DIR}/${base}_filtered.bam"; then
        log "ERROR" "❌ Error: MAPQ filtering failed for ${base}"
        record_metrics "PROCESS FAILED: MAPQ filtering - $base"
        continue
    fi
    stop_monitoring "mapq_filter"
    
 #b.  Count reads after MAPQ filter
 reads_mapq=$(samtools view -c "${CLEANED_DIR}/${base}_filtered.bam") #=====================================================================
# 2. SPIKE FILTERING:Count spike-in reads (based on contig names starting with "spike_")
#=====================================================================  
if [[ "$spike_type" != "none" ]]; then
    log "INFO" "2. Counting and removing spike-in reads..."

    spike_reads=$(samtools view -h "${CLEANED_DIR}/${base}_filtered.bam" | awk '$1 ~ /^@/ || $3 ~ /^spike_/' | samtools view -c -)

# Save spike-only BAM (optional)
    samtools view -h "${CLEANED_DIR}/${base}_filtered.bam" | awk '$1 ~ /^@/ || $3 ~ /^spike_/' | samtools view -b -o "${CLEANED_DIR}/${base}_spike.bam"

# Remove spike reads to produce a host-only BAM (overwrite filtered.bam)
    samtools view -h "${CLEANED_DIR}/${base}_filtered.bam" | awk '$1 ~ /^@/ || $3 !~ /^spike_/' | samtools view -b -o "${CLEANED_DIR}/${base}_host.bam"

    mv "${CLEANED_DIR}/${base}_host.bam" "${CLEANED_DIR}/${base}_filtered.bam"

    log "INFO" "🧼 Stripping spike contigs from BAM header..."
    samtools view -H "${CLEANED_DIR}/${base}_filtered.bam" | grep -v '^@SQ.*spike_' > "${CLEANED_DIR}/${base}_header.sam"
    samtools reheader "${CLEANED_DIR}/${base}_header.sam" "${CLEANED_DIR}/${base}_filtered.bam" > "${CLEANED_DIR}/${base}_filtered_no_spike_header.bam"
    mv "${CLEANED_DIR}/${base}_filtered_no_spike_header.bam" "${CLEANED_DIR}/${base}_filtered.bam"
    rm "${CLEANED_DIR}/${base}_header.sam"

# Count host reads after spike removal
    host_reads=$(samtools view -c "${CLEANED_DIR}/${base}_filtered.bam")

# Compute scaling factor (e.g., per 1 million spike reads)
    if [[ "$spike_reads" -gt 0 ]]; then
        scale=$(awk -v spike="$spike_reads" 'BEGIN { printf "%.4f", 1000000 / spike }')
    else
        scale="NA"
    fi

else
    log "INFO" "2. 🧫 No spike-in declared for $base — skipping spike filtering"
    spike_reads=0
    host_reads=$(samtools view -c "${CLEANED_DIR}/${base}_filtered.bam")
    scale="NA"
fi

# === Safe normalization recording (resumable) ===
scaling_tmp="${METRICS_DIR}/${base}_scaling.tmp"

if [[ -f "$scaling_tmp" ]]; then
    log "INFO" "⏩ Skipping scaling record — already exists for $base"
else
    echo -e "${base}\t${spike_reads}\t${host_reads}\t${scale}" > "$scaling_tmp"
    log "INFO" "🔬 Spike reads: $spike_reads | Host reads: $host_reads | Scale: $scale"
    log "INFO" "📄 Scaling record saved to: $scaling_tmp"
fi
            
#=====================================================================
# 3. BLACKLIST REMOVAL (EXISTING)
#=====================================================================
    log "INFO" "3.🚫 Removing blacklist regions..."
    start_monitoring "blacklist_filter" $$
    if ! bedtools intersect -v -abam "${CLEANED_DIR}/${base}_filtered.bam" -b "$BLACKLIST" > "${CLEANED_DIR}/${base}_noblack.bam"; then
        log "ERROR" "❌ Error: Blacklist removal failed for ${base}"
        record_metrics "PROCESS FAILED: Blacklist removal - $base"
        continue
    fi
    reads_noblack=$(samtools view -c "${CLEANED_DIR}/${base}_noblack.bam")
    stop_monitoring "blacklist_filter"
    
#=====================================================================
# Step 4: Remove mito & make clean BAM
#=====================================================================
    log "INFO" "4.🔋 Separating mitochondrial reads..."
    if ! samtools view -h "${CLEANED_DIR}/${base}_noblack.bam" | \
        awk '$1 ~ /^@/ || ($3 == "chrM" || $3 == "MT")' | \
        samtools view -b -o "${CLEANED_DIR}/${base}_mito.bam"; then
        log "WARN" "⚠️ Warning: Mitochondrial separation failed for ${base}"
    fi

#=====================================================================
# Step 5: Creating clean BAM
#=====================================================================
    
    log "INFO" "5. 🚿 Creating final clean BAM..."
    start_monitoring "final_clean" $$
    if ! samtools view -h "${CLEANED_DIR}/${base}_noblack.bam" | \
        awk '$1 ~ /^@/ || ($3 != "chrM" && $3 != "MT")' | \
        samtools view -b -o "${CLEANED_DIR}/${base}_clean.bam"; then
        log "ERROR" "❌ Error: Final BAM creation failed for ${base}"
        record_metrics "PROCESS FAILED: Final BAM creation - $base"
        continue
    fi
    stop_monitoring "final_clean"
    
#=====================================================================
# Step 6: Indexing Clean BAM
#=====================================================================
   log "INFO" "6. 🔖 Indexing final BAM..."
   
   # Remove any old or invalid .bai index to prevent issues
   rm -f "${CLEANED_DIR}/${base}_clean.bam.bai"
   
   start_monitoring "indexing" $$
    if ! samtools index "${CLEANED_DIR}/${base}_clean.bam"; then
        log "ERROR" "❌ Error: Indexing failed for ${base}"
        record_metrics "PROCESS FAILED: Indexing - $base"
        continue
    fi
    stop_monitoring "indexing" $$
    
#=====================================================================
# Step 6b: Index Stats Summary (Verify spike removal)
#=====================================================================
log "INFO" "6b. 🧬 Checking BAM index stats for contig verification..."
idxstats_out="${METRICS_DIR}/${base}_idxstats.tsv"

if ! samtools idxstats "${CLEANED_DIR}/${base}_clean.bam" > "$idxstats_out"; then
    log "WARN" "⚠️ samtools idxstats failed for ${base}"
else
    if grep -q '^spike_' "$idxstats_out"; then
        spike_line=$(grep '^spike_' "$idxstats_out" | head -1)
        log "WARN" "⚠️ Spike contigs still present in index for ${base}: $spike_line"
    else
        log "INFO" "✅ No spike contigs detected in ${base}_clean.bam index"
    fi
fi
    
    
    
#=====================================================================
# Step 7: QC Metrics
#=====================================================================
    log "INFO" "7. 🧠 Generating QC metrics..."
    if ! samtools flagstat "${CLEANED_DIR}/${base}_clean.bam" > "${METRICS_DIR}/${base}_clean_flagstat.txt"; then
        log "WARN" "⚠️ Warning: Flagstat failed for ${base}"
    fi
    
    if ! qualimap bamqc \
        -bam "${CLEANED_DIR}/${base}_clean.bam" \
        -outdir "${METRICS_DIR}/${base}_qualimap_clean" \
        -outformat PDF:HTML; then
        log "WARN" "⚠️ Warning: Qualimap failed for ${base}"
    fi
    
#=====================================================================
# Step 8: Reporting
#=====================================================================
    reads_final=$(samtools view -c "${CLEANED_DIR}/${base}_clean.bam")
    retention_rate=$(awk "BEGIN {printf \"%.2f\", ($reads_final/$reads_mapq)*100}")
    
log "INFO" ""
log "INFO" "📊 ${base} processing summary:"
log "INFO" "  🧾 Initial reads:       $(printf "%'d" $reads_initial)"
log "INFO" "  🔍 After MAPQ≥30:       $(printf "%'d" $reads_mapq)"
log "INFO" "  🧬 After spike filter:  $(printf "%'d" $host_reads)"
log "INFO" "  🚫 After blacklist:     $(printf "%'d" $reads_noblack)"
log "INFO" "  ✅ Final retained:      $(printf "%'d" $reads_final) (${retention_rate}%)"
    
    ((processed++))
    record_metrics "SAMPLE COMPLETE: $base (Retention: ${retention_rate}%)"
done

# === Final MultiQC ===
log "INFO" "🔗 Aggregating QC metrics with MultiQC..."
record_metrics "MULTIQC START"
multiqc "$METRICS_DIR" -o "${METRICS_DIR}/BAM_QC_MultiQC" || log "WARN" "⚠️ Warning: MultiQC encountered issues"
record_metrics "MULTIQC COMPLETE"

generate_performance_report

if "$DRY_RUN"; then
    log "INFO" "🧪 DRY-RUN COMPLETE. No files were modified."
    return 0
fi


# === Completion ===
log "INFO" "\n=============================================="
log "INFO" "✅ CLEANING COMPLETED"
log "INFO" "🧪 Processed: ${processed}/${total_samples} samples"
log "INFO" "⏱️  Elapsed time: $(get_elapsed_time)"
log "INFO" "📂 Output directory: $CLEANED_DIR"
log "INFO" "📊 Performance metrics saved to: ${PERF_LOG}"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "=============================================="
record_metrics "PROCESS COMPLETE: ${processed} samples"

log "INFO" "📎 Rebuilding normalization_factors.tsv from per-sample records..."

{
  echo -e "Sample_ID\tSpikeReads\tHostReads\tScalingFactor"
  # Prepare map of scaling factors (no header added here)
  cat "${METRICS_DIR}/"*"_scaling.tmp"
} > "${METRICS_DIR}/normalization_factors.tsv"

# Copy to .tmp so merging works (this is what your script expects!)
#cp "${METRICS_DIR}/normalization_factors.tsv" "${METRICS_DIR}/normalization_factors.tmp"
# Merge into original metadata
awk -F'\t' '
    NR==FNR {
        if (FNR > 1 || NF == 4) {
            spike[$1]=$2;
            host[$1]=$3;
            scale[$1]=$4;
        }
        next
    }
    FNR == 1 {
    print $0 "\tSpikeReads\tHostReads\tScalingFactor";
    next
}
    {
        s = ($1 in scale) ? scale[$1] : "NA";
        sr = ($1 in spike) ? spike[$1] : "NA";
        hr = ($1 in host) ? host[$1] : "NA";
        print $0 "\t" sr "\t" hr "\t" s;
    }
' "${METRICS_DIR}/normalization_factors.tsv" metadata/mapping.tsv > metadata/mapping_scaled.tsv

# Clean up temp file
#rm -f "${METRICS_DIR}/normalization_factors.tmp"


SPIKE_QC_FLAG="${METRICS_DIR}/.spike_qc_done"

if [[ ! -f "$SPIKE_QC_FLAG" ]]; then
    log "INFO" "🔎 Performing spike-in quality check..."

    awk -F'\t' 'NR > 1 {
        sample=$1
        spike=$2 + 0
        host=$3 + 0
        scale=$4 + 0

        if (spike < 5000) {
            printf "[WARNING] ⚠️ Sample %s has low spike-in reads: %d\n", sample, spike
        } else if (spike > 50000) {
            printf "[WARNING] ⚠️ Sample %s has high spike-in reads: %d (possible contamination?)\n", sample, spike
        }

        if (scale != "NA") {
            if (scale < 10) {
                printf "[WARNING] ⚠️ Sample %s has low scale factor: %.2f\n", sample, scale
            } else if (scale > 100) {
                printf "[WARNING] ⚠️ Sample %s has high scale factor: %.2f (likely low spike)\n", sample, scale
            }
        }
    }' "${METRICS_DIR}/normalization_factors.tsv" | tee "${METRICS_DIR}/spike_qc_report.txt"

    touch "$SPIKE_QC_FLAG"
    log "INFO" "✅ Spike QC complete — flag created: $SPIKE_QC_FLAG"
    log "INFO" "📄 Spike QC report saved to: ${METRICS_DIR}/spike_qc_report.txt"
else
    log "INFO" "⏭️ Spike QC already completed — skipping"
fi
#Create spike analysis graph
log "INFO" "📊 Spike plot in production ..."
Rscript modules/pipeline1/10_plot_spike_qc_summary.R metadata/mapping_scaled.tsv results/QC_spike_plots



log "INFO" "📄 Final scaling table: ${METRICS_DIR}/normalization_factors.tsv"
log "INFO" "📄 Enriched metadata: metadata/mapping_scaled.tsv"

} 
parse_args "$@"
# This logs everything `main` prints
main 

