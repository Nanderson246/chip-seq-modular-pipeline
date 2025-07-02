#!/usr/bin/env bash
###############################################################################
# bash modules/util/setup_reference.sh – download & index hg38 / mm10 assets for ChIP-seq
#
# Author: Nancy Anderson • 2025-06-22
###############################################################################

 # Install hg38 reference into the package's default Reference folder
# bash modules/utils/setup_reference.sh -g hg38
# 
# # Install mm10 reference into a custom folder
# ./setup_reference.sh -g mm10 -d /data/genomes/mm10
# Usage: setup_reference.sh -g <hg38|mm10> [-d REFERENCE_DIR]
# 
# -g   Genome build (required):  hg38  or  mm10
# 
# -d   Reference folder (optional)
# If omitted, files will be created under:
#   chip_seq_modular/Reference/<GENOME>
# 
# This script downloads and prepares:
#   • Reference FASTA
# • BWA index
# • FASTA index (.fai)
# • Picard dictionary (.dict)
# • ENCODE blacklist BED
# • TSS annotation BED (from GENCODE GTF)
###############################################################################


set -euo pipefail

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



###############################################################################
# 1.  Resolve package layout ----------------------------
###############################################################################
SCRIPT_PATH="$(readlink -f "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT_PATH")"            # .../chip_seq_modular/modules/utils
PROJECT_ROOT="$(dirname "$(dirname "$SCRIPT_DIR")")"   # .../chip_seq_modular

TOOLS_DIR="${PROJECT_ROOT}/tools"
DEFAULT_REF_DIR="${PROJECT_ROOT}/Reference"
LOG_DIR="${PROJECT_ROOT}/logs/pipeline1"
LOG_FILE="${LOG_DIR}/setup_reference${TIMESTAMP}.log"

###############################################################################
# 2.  Picard – UNCHANGED -------------------------------------------------------
###############################################################################
picard() { java -jar "${TOOLS_DIR}/picard.jar" "$@"; }
export -f picard

###############################################################################
# 3.  Command-line interface ---------------------------------------------------
###############################################################################
usage() {
    cat <<EOF
Usage: $0 -g <hg38|mm10> [-d REFERENCE_DIR]

  -g   Genome build (required):  hg38  or  mm10
  -d   Reference folder (optional, default: ${DEFAULT_REF_DIR})

Example
  $0 -g hg38            # installs into ${DEFAULT_REF_DIR}/hg38
  $0 -g mm10 -d /data/mm10_refs
EOF
    exit 1
}

GENOME=""; REF_DIR=""
while getopts ":g:d:h" opt; do
    case "$opt" in
        g) GENOME="$OPTARG" ;;
        d) REF_DIR="$OPTARG" ;;
        h|*) usage ;;
    esac
done
[[ -z $GENOME ]] && { echo "❌  -g is required"; usage; }
[[ $GENOME != "hg38" && $GENOME != "mm10" ]] && { echo "❌  -g must be hg38 or mm10"; exit 1; }
[[ -z $REF_DIR ]] && REF_DIR="$DEFAULT_REF_DIR"


###############################################################################
# 4.  URLs for each genome -----------------------------------------------------
###############################################################################
if [[ $GENOME == "hg38" ]]; then
    FASTA_URL="http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
    BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz"
    GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/gencode.v45.annotation.gtf.gz"
else
    FASTA_URL="http://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz"
    BLACKLIST_URL="https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz"
    GTF_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M31/gencode.vM31.annotation.gtf.gz"
fi

FASTA_FILE=$(basename "$FASTA_URL" .gz)
BLACKLIST_BED="${GENOME}-blacklist.bed"
TSS_BED="tss_annotations.bed"

###############################################################################
# 5.  Logging & workspace ------------------------------------------------------
###############################################################################
mkdir -p "${REF_DIR}/${GENOME}" "${LOG_DIR}"
cd       "${REF_DIR}/${GENOME}"

exec > >(tee -a "${LOG_DIR}/Reference_setup_${GENOME}_$(date +%F_%H-%M-%S).log") 2>&1

#--- BANNER ---
log "INFO" ""
log "INFO" "=============================================="
log "INFO" "🧬 Module: ${SCRIPT_BASE_NAME}"
log "INFO" "📌 Purpose: Set Genome reference, black list, and TSS annotations"
log "INFO" "📁 Output dir: ${REF_DIR}/"
log "INFO" "🕒 Start time: $(date '+%F %T')"
log "INFO" "📦 Script version: ${VERSION} (${SCRIPT_NAME})"
log "INFO" "=============================================="

start_time=$(date +%s)
echo "[INFO] Starting reference setup for $GENOME in ${REF_DIR}/${GENOME}"

###############################################################################
# 6.  Downloads & indexing -----------------------------------------------------
###############################################################################
echo "[INFO] ▶ Genome FASTA"
wget -c "$FASTA_URL"
gunzip -f "$FASTA_FILE.gz"

echo "[INFO] ▶ FASTA index (.fai)"
samtools faidx "$FASTA_FILE"

echo "[INFO] ▶ Sequence dictionary (.dict)"
picard CreateSequenceDictionary R="$FASTA_FILE" O="${FASTA_FILE%.fa}.dict"

echo "[INFO] ▶ Blacklist BED"
wget -c "$BLACKLIST_URL"
gunzip -fc "$(basename "$BLACKLIST_URL")" > "$BLACKLIST_BED"

echo "[INFO] ▶ BWA index"
bwa index "$FASTA_FILE"

###############################################################################
# 7.  TSS annotation BED -------------------------------------------------------
###############################################################################
echo "[INFO] ▶ GENCODE GTF"
wget -c "$GTF_URL"
GTF_FILE=$(basename "$GTF_URL")                           # e.g., gencode.v45.annotation.gtf.gz
GTF_UNZIPPED="${GTF_FILE%.gz}"                            # e.g., gencode.v45.annotation.gtf

echo "[INFO] ▶ Unzipping GTF to: $GTF_UNZIPPED"
gunzip -c "$GTF_FILE" > "$GTF_UNZIPPED"

echo "[INFO] ▶ Generating $TSS_BED from $GTF_UNZIPPED"
awk 'BEGIN{OFS="\t"}
     $3=="transcript" {
       match($0,/gene_id "([^"]+)"/,gid);
       match($0,/gene_name "([^"]+)"/,gname);
       full=gid[1]; id=full; sub("\\..*","",id);
       if($7=="+"){s=$4-1;e=$4}else{s=$5-1;e=$5}
       print $1,s,e,full,id,gname[1],$7
     }' "$GTF_UNZIPPED" \
  | sort -k1,1 -k2,2n > "$TSS_BED"

###############################################################################
# 8.  Done --------------------------------------------------------------------
###############################################################################
end_time=$(date +%s)
mins=$(( (end_time - start_time)/60 )); secs=$(( (end_time - start_time)%60 ))

log "INFO" "✅ Reference assets ready in: ${REF_DIR}/${GENOME}"
log "INFO" "├─ $FASTA_FILE (+.fai and .dict)"
log "INFO" "├─ $BLACKLIST_BED"
log "INFO" "├─ BWA index files"
log "INFO" "└─ $TSS_BED"
log "INFO" "⏱️  Elapsed: ${mins}m ${secs}s"
log "INFO" "=============================================="
log "INFO" "✅ MODULE COMPLETED"
log "INFO" "🕒 End time: $(date '+%F %T')"
log "INFO" "=============================================="
