#!/usr/bin/env bash
# Module: 06_igv_snapshot.sh 
# Description: Perform igv snapshot of peaks you want.
# Author: Nancy Anderson

# ###########################################################################
# Software Requirements for 08_igv_batch_snapshots.sh
#
# Required:
#   - bash              (any POSIX-compliant shell should work)
#   - awk               (for row filtering)
#   - sort, head        (for filtering top regions)
#   - bigWigSummary     (from UCSC tools, for BigWig signal filtering)
#   - IGV (command-line) (Java GUI-based genome viewer, needs X11 if GUI)
#
# Assumptions:
#   - IGV is callable from a CLI entry point (e.g. /usr/bin/igv or igv.sh)
#   - BED files are IGV-compatible (5+ columns: chr, start, end, name, score)
#   - BigWig files (.bw) exist for signal filtering (optional)
#
# Optional:
#   - X11 display (if IGV needs to launch a GUI; not required for headless batch runs on some setups)
#
# Notes:
#   - If using this in a Docker container, ensure bigWigSummary is installed or included via bind-mount.
#
# ###########################################################################
# 08_igv_batch_snapshots.sh v5 — Automated IGV snapshot generator
#
# Filters peaks by chromosome, score, length, and optional BigWig coverage.
# Supports genome checks (e.g. mm10), region caps, and logs progress.
#
# USAGE EXAMPLE:
# bash modules/pipeline2/06_igv_batch_snapshots.sh /usr/bin/igv mm10 macs3 \
#   --chr chr1 --min-score 300 --max-len 3000 --min-bw 2 --max-regions 100
# ###########################################################################
# USAGE:
#   bash modules/pipeline2/06_igv_batch_snapshots.sh /path/to/igv GENOME TAG [options]
#
# REQUIRED:
#   /path/to/igv     Path to IGV launcher (e.g., /usr/bin/igv, igv.sh)
#   GENOME           Genome ID used in IGV (e.g., hg38, mm10, hg19)
#   TAG              IDR or peak-caller tag to locate BED files (e.g., macs3, homer)
#
# OPTIONS:
#   --chr CHR             Only keep peaks from given chromosomes (comma-separated)
#   --min-score SCORE     Keep peaks with MACS/HOMER score ≥ SCORE
#   --max-len LENGTH      Keep peaks with length ≤ LENGTH (in bp)
#   --min-bw VALUE        Keep only regions with ≥ VALUE signal in any replicate BigWig
#   --max-regions N       Limit output to top N regions (default: 100)
#
# EXAMPLES:
#   bash 08_igv_batch_snapshots.sh /usr/bin/igv mm10 macs3 \
#        --chr chr1 --min-score
# ###########################################################################

set -euo pipefail

# --------------------------- 1. ARGUMENTS ------------------------------------
IGV_BIN=${1:-}    GENOME=${2:-}    IDR_TAG=${3:-}
shift 3 || true

CHR_FILTER=""     MIN_SCORE=""     MAX_LEN=""
MIN_BW=""         MAX_REGIONS="100"

while [[ $# -gt 0 ]]; do
  case $1 in
    --chr)         CHR_FILTER=$2; shift 2 ;;
    --min-score)   MIN_SCORE=$2;  shift 2 ;;
    --max-len)     MAX_LEN=$2;    shift 2 ;;
    --min-bw)      MIN_BW=$2;     shift 2 ;;
    --max-regions) MAX_REGIONS=$2;shift 2 ;;
    *) echo "[ERROR] Unknown flag: $1" >&2; exit 1 ;;
  esac
done

[[ -x $IGV_BIN ]] || { echo "[ERROR] IGV binary not executable: $IGV_BIN"; exit 1; }

# ------------------------ 2. GENOME SAFETY CHECK -----------------------------
if [[ $GENOME != "hg19" && $GENOME != "hg38" && $GENOME != "mm10" ]]; then
  echo "[WARN] Genome label '$GENOME' may not be recognized by IGV."
  echo "       Common options: hg19, hg38, mm10"
fi

# ------------------------ 3. PATH SETUP --------------------------------------
OUT_ROOT="analysis/ChIPseeker_TSS_Hommer_IDR_annotation/${IDR_TAG}"
BED_DIR="${OUT_ROOT}/Enriched_cluster/IGV_ready_bed"
SNAP_ROOT="${OUT_ROOT}/Enriched_cluster/IGV_snapshots"
BW_DIR="analysis/Replicate_QC/bigwig"
META_DIR="metadata"
mkdir -p "$SNAP_ROOT"

# ------------------------ 4. BED FILTER CONSTRUCTION -------------------------
build_awk_filter() {
  local f=""
  [[ -n $CHR_FILTER ]] && {
    IFS=',' read -ra arr <<< "$CHR_FILTER"
    pat=$(printf "|%s" "${arr[@]}")
    f+="(\$1 ~ /^(${pat:1})\$/)"
  }
  [[ -n $MIN_SCORE ]] && [[ -n $f ]] && f+=" && "; f+="(\$5 >= $MIN_SCORE)"
  [[ -n $MAX_LEN   ]] && [[ -n $f ]] && f+=" && "; f+="(\$3-\$2 <= $MAX_LEN)"
  [[ -z $f ]] && echo '{print}' || echo "{ if ($f) print }"
}
AWK_FILTER=$(build_awk_filter)

# ------------------------ 5. MAIN LOOP ---------------------------------------
for bed in "${BED_DIR}"/*.bed; do
  [[ -s $bed ]] || continue
  sample=$(basename "$bed" .bed)
  snap_dir="${SNAP_ROOT}/${sample}"
  mkdir -p "$snap_dir"
  batch_file="${snap_dir}/${sample}.bat"
  tmp_bed="$(mktemp)"

  # Filter and cap regions
  awk "$AWK_FILTER" "$bed" | sort -k5,5nr | head -n "$MAX_REGIONS" > "$tmp_bed"
  [[ ! -s $tmp_bed ]] && { echo "[SKIP] $sample — no valid peaks"; rm "$tmp_bed"; continue; }

  # Optional: BigWig intensity filtering
  declare -a bw_paths=()
  if [[ -n $MIN_BW ]]; then
    shortname=${sample%%_*}
    group=$(awk -F'\t' -v sn="$shortname" '$1==sn{print $4}' "${OUT_ROOT}/filename_mapping.tsv")
    mapfile -t bw_paths < <(awk -F'\t' -v g="$group" '$0~g{print $1}' "${META_DIR}/mapping_filtered.tsv" |
                            sed "s|^|${BW_DIR}/|;s|\$|.bw|")
    filter_bw() {
      while read -r chr s e name score rest; do
        keep=0
        for bw in "${bw_paths[@]}"; do
          [[ -s $bw ]] || continue
          val=$(bigWigSummary "$bw" "$chr" "$s" "$e" 1 2>/dev/null || echo 0)
          awk 'BEGIN{if('"$val"' >= '"$MIN_BW"') exit 0; else exit 1}'
          [[ $? -eq 0 ]] && { keep=1; break; }
        done
        (( keep )) && echo -e "$chr\t$s\t$e\t$name\t$score"
      done
    }
    <"$tmp_bed" filter_bw > "${tmp_bed}.bw" && mv "${tmp_bed}.bw" "$tmp_bed"
  fi

  [[ ! -s $tmp_bed ]] && { echo "[SKIP] $sample — no peaks passed filters"; rm "$tmp_bed"; continue; }
  region_count=$(wc -l <"$tmp_bed")

  # Generate IGV batch file
  {
    echo "new"
    echo "genome ${GENOME}"
    echo "load $bed"

    if [[ ${#bw_paths[@]} -eq 0 ]]; then
      shortname=${sample%%_*}
      group=$(awk -F'\t' -v sn="$shortname" '$1==sn{print $4}' "${OUT_ROOT}/filename_mapping.tsv")
      mapfile -t bw_paths < <(awk -F'\t' -v g="$group" '$0~g{print $1}' "${META_DIR}/mapping_filtered.tsv" |
                              sed "s|^|${BW_DIR}/|;s|\$|.bw|")
    fi
    for bw in "${bw_paths[@]}"; do [[ -s $bw ]] && echo "load $bw"; done

    echo "snapshotDirectory $snap_dir"
    awk '{printf "goto %s:%s-%s\nsnapshot %s.png\n", $1,$2,$3,$4}' "$tmp_bed"
    echo "exit"
  } > "$batch_file"

  echo "[IGV] $sample → $region_count regions"
  "$IGV_BIN" -b "$batch_file" >/dev/null 2>&1 || echo "[WARN] IGV failed for $sample"
  rm "$tmp_bed"
done

echo "[DONE] Snapshots saved in $SNAP_ROOT"

