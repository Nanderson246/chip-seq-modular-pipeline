#!/usr/bin/env bash
###############################################################################
# bash modules/util2/setup_spikein_refs.sh –
# Author: Nancy Anderson • 2025-06-22
###############################################################################
# spikein_setup.sh – Set up spike-in genome references for ChIP-seq
#
# Author: Nancy Anderson • 2025-06-22
#
# USAGE:
#   bash modules/utils/spikein_setup.sh
#
# DESCRIPTION:
#   This script downloads, decompresses, indexes (BWA), and organizes a set of
#   spike-in genomes used in ChIP-seq normalization.
#
#   For each genome, it:
#     • Downloads the reference (FASTA or ZIP)
#     • Extracts it to:       SpikeinReference/<GENOME>/<GENOME>.fa
#     • Builds BWA index
#     • Creates symlinks in:  SpikeinReference/all/
#     • Logs all paths to:    SpikeinReference/spikeins.yaml
#
# OUTPUT:
#   • Indexed FASTA files
#   • Central symlink folder for easy access
#   • YAML file with paths to use in the pipeline
#
# LOCATION:
#   This script lives in:     modules/utils/
#   It writes to:             SpikeinReference/ (in the top-level project dir)
###############################################################################

###############################################################################
# USAGE:
#   From project root:       bash modules/utils/setup_spikein_refs.sh
#   From inside modules/:    bash ../utils/setup_spikein_refs.sh
#   From anywhere:           bash /full/path/to/modules/utils/setup_spikein_refs.sh
###############################################################################
set -euo pipefail

# === Absolute paths (must come before any usage) ===
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"

# Extract root path (e.g., if script is at project_root/modules/utils)
PACKAGE_ROOT="${SCRIPT_DIR%/modules/utils}"
BASE_DIR="${PACKAGE_ROOT}/SpikeinReference"


# ── Dependency check ──────────────────────────────────────────────────────
REQUIRED_TOOLS=(wget gunzip bwa unzip)

for tool in "${REQUIRED_TOOLS[@]}"; do
    command -v "$tool" >/dev/null || {
        echo "❌ Required tool '$tool' not found in PATH"
        exit 1
    }
done

# ── Optional logging to file ──────────────────────────────────────────────
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
LOG_DIR="${BASE_DIR}/logs/pipeline1"
LOG_FILE="${LOG_DIR}/spikein_setup_${TIMESTAMP}.log"
mkdir -p "$LOG_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1

# === SCRIPT ===
readonly VERSION="2.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly SCRIPT_BASE_NAME="${SCRIPT_NAME%.*}" 
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)

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


# Now build clean, absolute paths
ALL_LINK_DIR="${BASE_DIR}/all"
YAML_OUT="${BASE_DIR}/spikeins.yaml"

mkdir -p "$BASE_DIR" "$ALL_LINK_DIR"

# ── Optional logging to file ──────────────────────────────────────────────
LOG_DIR="${BASE_DIR}/logs/pipeline1"
LOG_FILE="${LOG_DIR}/spikein_setup_${TIMESTAMP}.log"
mkdir -p "$LOG_DIR"
exec > >(tee -a "$LOG_FILE") 2>&1


#--- BANNER ---
log "INFO" ""
log "INFO" "=============================================="
log "INFO" "🧬 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Purpose: Set up spike-in genome references for ChIP-seq"
log "INFO" "📁 Output dir: ${BASE_DIR}/"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "📦 Script version: ${VERSION} (${SCRIPT_NAME})"
log "INFO" "=============================================="


printf "# Spike-in genome reference paths\n\n" > "$YAML_OUT"

setup_spikein() {
    local name=$1
    local url=$2
    local fasta_out="$BASE_DIR/$name/$name.fa"
    local genome_dir="$BASE_DIR/$name"

    echo "📥 Setting up $name from $url..."
    mkdir -p "$genome_dir"

    # Download
    wget -q -O "$genome_dir/source" "$url" || { echo "❌ Download failed for $name"; return 1; }

    # Decompress and rename
    if [[ "$url" == *.zip ]]; then
        unzip -q "$genome_dir/source" -d "$genome_dir"
        mv "$genome_dir"/*.fa "$fasta_out"
    else
        gunzip -c "$genome_dir/source" > "$fasta_out"
    fi

    # ✅ Check FASTA was created
    if [[ ! -f "$fasta_out" ]]; then
        echo "❌ FASTA not found for $name after extraction"
        return 1
    fi

    # Index
    echo "🧬 Indexing $name..."
    bwa index "$fasta_out"

# Link FASTA and index files: Eliminated for permission level disruption when folder is moved
#    echo "🔗 Linking $name to $ALL_LINK_DIR..."
#    for ext in "" .amb .ann .bwt .pac .sa; do
#       ln -sf "$fasta_out$ext" "$ALL_LINK_DIR/$(basename "$fasta_out")$ext"
#    done

    # YAML entry
    echo "$name:" >> "$YAML_OUT"
    echo "  fasta: \"$fasta_out\"" >> "$YAML_OUT"
    echo "  index_prefix: \"$fasta_out\"" >> "$YAML_OUT"
}

# Setup genomes
setup_spikein "dm6" "https://hgdownload.soe.ucsc.edu/goldenPath/dm6/bigZips/dm6.fa.gz"
setup_spikein "sacCer3" "https://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz"
setup_spikein "ERCC" "https://assets.thermofisher.com/TFS-Assets/LSG/manuals/ERCC92.zip"
setup_spikein "phiX174" "https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz"
setup_spikein "mm10" "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"


log "INFO" "📦 Spike-in references successfully prepared at: $BASE_DIR"
log "INFO" "📝 YAML summary saved to: $YAML_OUT"
log "INFO" "✅ Reference assets ready in: ${BASE_DIR}"
log "INFO" "├── all"
log "INFO" "├── dm6"
log "INFO" "├── ERCC"
log "INFO" "├── mm10"
log "INFO" "├── phiX174"
log "INFO" "├── sacCer3"
log "INFO" "└── spikeins.yaml"


