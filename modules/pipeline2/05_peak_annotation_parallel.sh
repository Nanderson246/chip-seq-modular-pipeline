#!/usr/bin/env bash
# -----------------------------------------------------------------------------
# Name:        05_peak_annotation_dual_WITH_RIGHT_IDR_parallel.sh
# Version:     3.0.0  (2025‑06‑18)
# Author:      Nancy Anderson  
# Purpose:     Fast dual‑mode (HOMER / MACS3) peak annotation with per‑Bed
#              parallelism, cached GTF, safe logging and TSS proximity.

## === Software Requirements ===
# Required software/tools (must be in $PATH):
#   - bash (>= 4.0)
#   - coreutils (cut, grep, awk, sort, md5sum, etc.)
#   - bedtools (>= 2.29)
#   - HOMER (>= 4.11)
#       - Required: annotatePeaks.pl
#   - GNU parallel (>= 20210422)
#   - R (>= 4.0) with required CRAN and Bioconductor packages
#       - Run: modules/pipeline2/cluster_enrichment_updated_hg_mice.R
#       - See script header for exact package list (auto-installs on demand)

# Optional:
#   - tee (for log redirection and streaming output)
# ------------------------------------------------------------------------------
# =============================================================================
# Script:     cluster_enrichment.R
# Purpose:    Enrichment annotation (GO, KEGG, Reactome) from pre-annotated peak files
# Author:     Nancy Anderson
# Version:    3.0 (2025-06-XX)
#
# === Software Requirements ===
# R (>= 4.0.0) with the following CRAN and Bioconductor packages:
#
# CRAN:
#   - httr
#   - jsonlite
#   - tidyr
#   - data.table
#   - stringr
#   - dplyr
#   - readr
#   - tools
#   - ggplot2
#   - gprofiler2
#
# Bioconductor:
#   - clusterProfiler
#   - enrichplot
#   - ReactomePA
#   - reactome.db
#   - AnnotationDbi
#   - org.Hs.eg.db
#   - org.Mm.eg.db
#
# Notes:
#   - All packages are auto-installed if missing (no user input required).
#   - NO internet-based biomaRt API queries are performed.
#   - Ensembl-to-Entrez conversion is done via local `AnnotationDbi` databases.
#   - Input files must include or allow extraction of Ensembl IDs.
# =============================================================================

# ------------------------------------------------------------------------------

# Reference data:
#   - Preprocessed genome FASTA, GTF, and TSS BED files for hg38 or mm10
#     in `Reference/<genome>/` directory (e.g., Reference/hg38/gencode.v45.annotation.gtf)

# Tested environment:
#   - Linux (Ubuntu 20.04+)
#   - Not tested on macOS (due to missing GNU coreutils by default)

##############################################################################
# Usage:
#Default genome hg38
#   bash modules/pipeline2/05_peak_annotation_parallel.sh --idr macs3 --threads 4 --genome mm10
#   or
#   bash modules/pipeline2/05_peak_annotation_parallel.sh \
#        --idr hommer   
# -----------------------------------------------------------------------------
set -uo pipefail

################################################################################
#                               ── 1  SET‑UP ──                                #
################################################################################
readonly VERSION="3.0.0"
readonly SCRIPT_NAME=$(basename "$0")
readonly TIMESTAMP=$(date +%Y%m%d_%H%M%S)
THREADS=4                  # default; override with --threads
idr="homer"               # default IDR caller

# ── colours for pretty logs (TTY only) ───────────────────────────────────────
if [[ -t 1 ]]; then RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[0;33m'; BLUE='\033[0;34m'; NC='\033[0m'; else RED=''; GREEN=''; YELLOW=''; BLUE=''; NC=''; fi

log() { printf "${BLUE}[%s]${NC} %s\n" "$(date '+%F %T')" "$*"; }
err() { printf "${RED}[%s] ERROR:${NC} %s\n" "$(date '+%F %T')" "$*" >&2; }

usage() {
  cat <<EOF
Usage: $SCRIPT_NAME --idr [homer|macs3] [--threads N] --genome mm10
EOF
  exit 1
}

# ── parse CLI ────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
  case "$1" in
    --idr)     idr="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    --genome)
               GENOME_NAME="$2"; shift 2;;
    -h|--help) usage;;
    *) err "unknown option $1"; usage;;
  esac
done

[[ "$idr" =~ ^(homer|macs3)$ ]] || err "--idr must be homer or macs3"
command -v annotatePeaks.pl >/dev/null || err "HOMER not in PATH"
command -v bedtools         >/dev/null || err "bedtools not in PATH"
command -v parallel         >/dev/null || err "GNU parallel not in PATH"

################################################################################
#                         ── 2  PATHS & ENV VARS  ──                          #
################################################################################
IDR_ROOT="analysis/IDR_Results"
homer_dir="$IDR_ROOT/homer"
macs3_dir="$IDR_ROOT/macs3"
IDR_DIR=$([ "$idr" == "homer" ] && echo "$homer_dir" || echo "$macs3_dir")
# === Default CLI values ===
# === Default CLI values ===
GENOME_NAME="${GENOME_NAME:-hg38}"  # default if not passed via CLI
GENOME_DIR="Reference/${GENOME_NAME}"
GENOME="${GENOME_DIR}/${GENOME_NAME}.fa"

# Dynamically detect GTF file based on genome folder
GTF=$(find "$GENOME_DIR" -maxdepth 1 -type f -name "*.annotation.gtf" | head -n 1)
[[ -f "$GTF" ]] || err "GTF annotation file not found in $GENOME_DIR"

# Detect or set TSS annotation BED file
TSS_BED="${TSS_BED:-$GENOME_DIR/tss_annotations.bed}"
[[ -f "$TSS_BED" ]] || err "TSS BED file not found at $TSS_BED"


OUT_ROOT="analysis/ChIPseeker_TSS_Hommer_IDR_annotation/$idr"
mkdir -p "$OUT_ROOT"
LOG_DIR="logs/pipeline2"
LOG_FILE="$LOG_DIR/${idr}_peak_annotation_dual.log"
SUMMARY_LOG="$LOG_DIR/${idr}_summary_merged.log"


# ── cached GTF (parsed once) ─────────────────────────────────────────────────
ANNOT_CACHE="${GTF}.annot"
if [[ ! -s "$ANNOT_CACHE" ]]; then
  log "Building cached annotation $ANNOT_CACHE (first run only)"
  grep -v '^#' "$GTF" | cut -f1-9 | sort -k1,1 -k4,4n > "$ANNOT_CACHE"
fi

################################################################################
#                      ── 3  COLLECT INPUT BED FILES  ──                       #
################################################################################
log "Scanning $IDR_DIR for .sorted.bed files"
mapfile -t IDR_INPUT_FILES < <(find "$IDR_DIR/replicate_vs_replicate" \
                                     "$IDR_DIR/pooled_vs_replicate"     \
                                     "$IDR_DIR/merged_replicate_idr"    \
                                     -type f -name '*.sorted.bed')
[[ ${#IDR_INPUT_FILES[@]} -eq 0 ]] && err "No BED files found in $IDR_DIR"
log "Found ${#IDR_INPUT_FILES[@]} IDR peak files"

################################################################################
#                       ── 4  MAPPING TSV HEADER ──                            #
################################################################################
MAPPING_DIR="$OUT_ROOT/tmp_mapping" ; mkdir -p "$MAPPING_DIR"
MAPPING_FILE="$OUT_ROOT/filename_mapping.tsv"
printf 'shortname\toriginal_filename\tbiological_identity\tgroup\n' > "$MAPPING_FILE"
IDR_SUMMARY="$IDR_DIR/${idr}_idr_summary_full.tsv"

################################################################################
#                 ── 5  PARALLEL PEAK ANNOTATION FUNCTION ──                   #
################################################################################
annotate_one() {
  local IDR="$1"
  local base shortname bio_identity safe_bioid sample_dir annotated_file ann_log tmp_map group

  base=$(basename "$IDR" .sorted.bed)
  shortname=$(echo "$base" | md5sum | cut -d' ' -f1)
  bio_identity=$(echo "$base" | cut -d'_' -f1-4)
  safe_bioid=$(echo "$bio_identity" | tr ' /' '__')

  sample_dir="$OUT_ROOT/$safe_bioid" ; mkdir -p "$sample_dir"
  annotated_file="$sample_dir/${shortname}_${safe_bioid}.annotated_peaks.txt"

  tmp_map="$MAPPING_DIR/${shortname}.row"
  group="UNKNOWN"
  [[ -f "$IDR_SUMMARY" ]] && group=$(awk -v file="$base" 'BEGIN{FS="\t"} $0 ~ file {print $1; exit}' "$IDR_SUMMARY")

  # A) If annotation already exists, still write mapping row
  if [[ -s "$annotated_file" ]]; then
    printf '%s\t%s\t%s\t%s\n' "$shortname" "$base" "$bio_identity" "$group" > "$tmp_map"
    return 0
  fi

  # B) Run annotation
  ann_log="$sample_dir/${shortname}_${safe_bioid}.annotation_stats.log"
  annotatePeaks.pl "$IDR" "$GENOME" -gtf "$ANNOT_CACHE" \
      >  "$annotated_file" \
      2> "$ann_log"

  # C) Mapping row after annotation
  printf '%s\t%s\t%s\t%s\n' "$shortname" "$base" "$bio_identity" "$group" > "$tmp_map"
  log "Annotation for $base complete → $annotated_file"
}
export -f annotate_one
export GENOME ANNOT_CACHE OUT_ROOT IDR_SUMMARY MAPPING_DIR

parallel --verbose -j "$THREADS" --halt soon,fail=1 annotate_one ::: "${IDR_INPUT_FILES[@]}"

# concat mapping safely after parallel finished
cat "$MAPPING_DIR"/*.row >> "$MAPPING_FILE"
rm -r "$MAPPING_DIR"

################################################################################
#               ── 6  TSS PROXIMITY (also parallel, lighter) ──               #
################################################################################
TSS_OUTDIR="$OUT_ROOT/TSS" ; mkdir -p "$TSS_OUTDIR"
mapfile -t IDR_BED_FILES < <(find "$IDR_DIR/replicate_vs_replicate" \
                                   "$IDR_DIR/pooled_vs_replicate"     \
                                   "$IDR_DIR/merged_replicate_idr"    \
                                   -type f -name '*.sorted.bed')

closest_one() {
  local IDR="$1"
  [[ -s "$IDR" ]] || return 0
  local base=$(basename "$IDR" .sorted.bed)
  local shortname=$(echo "$base" | md5sum | cut -d' ' -f1)
  local bio_identity=$(echo "$base" | cut -d'_' -f1-4)
  local safe_bioid=$(echo "$bio_identity" | tr ' /' '__')
  local sample_tss_dir="$TSS_OUTDIR/$safe_bioid" ; mkdir -p "$sample_tss_dir"
  local outfile="$sample_tss_dir/${shortname}_${safe_bioid}.closest_TSS.bed"
  [[ -s "$outfile" ]] && return 0

  local tmp_input=$(mktemp)
  grep -E '^chr([0-9]{1,2}|X|Y|M)[[:space:]]' "$IDR" | sort -k1,1 -k2,2n > "$tmp_input"
  bedtools closest -a "$tmp_input" -b "$TSS_BED" -D ref > "${outfile}.tmp"
  { printf 'peak_chr\tpeak_start\tpeak_end\tTSS_chr\tTSS_start\tTSS_end\tgene_id_full\tgene_id_clean\tgene_name\tTSS_strand\tdistance_to_TSS\n';
    awk '$13 != "." && $13 != ""' "${outfile}.tmp"; } > "$outfile"
  rm "${outfile}.tmp" "$tmp_input"
}
export -f closest_one 
export TSS_OUTDIR TSS_BED

parallel --verbose -j "$THREADS" --halt soon,fail=1 closest_one ::: "${IDR_BED_FILES[@]}"

################################################################################
#                    ── 7  HEADER FIX, RELABEL & TSV COPY ──                   #
################################################################################
enriched_dir="$OUT_ROOT/Enriched"       # where the *.annotated_peaks_with_entrez.tsv live
mkdir -p "$enriched_dir"

find "$OUT_ROOT" -type f -name '*.annotated_peaks.txt' | while read -r file; do
  [[ -s "$file" ]] || continue

  header=$(head -n1 "$file")

  # 1. Fix header if it doesn't begin with 'PeakID'
  if ! grep -q '^PeakID' <<<"$header"; then
    cp "$file" "${file}.bak"
    {
      printf 'PeakID\tChr\tStart\tEnd\tStrand\tPeak_Score\tFocus_Ratio\tAnnotation\tDetailed_Annotation\tDistance_to_TSS\tNearest_PromoterID\tEnsembl_ID\tNearest_Unigene\tNearest_Refseq\tNearest_Ensembl\tGene_Name\tGene_Alias\tGene_Description\tGene_Type\n'
      tail -n +2 "${file}.bak"
    } > "$file"
    rm "${file}.bak"
  fi

  # 2. Relabel Entrez ID if needed
  entrez_col=$(echo "$header" | tr '\t' '\n' | grep -nx "Entrez ID" | cut -d: -f1 || true)
  if [[ -n "$entrez_col" ]]; then
    valid_entrez=$(awk -v col="$entrez_col" 'NR > 1 && $col ~ /^[0-9]+$/ { c++; if (c > 3) exit } END { print (c > 3 ? "yes" : "no") }' FS='\t' "$file")
    if [[ "$valid_entrez" == "no" ]]; then
      echo "[FIX] Re-labeling Entrez ID → Ensembl_ID in $file"
      sed -i '1s/Entrez ID/Ensembl_ID/' "$file"
    fi
  fi

  # 3. Create parallel .tsv version for enrichment script
  tsv_out="${file%.annotated_peaks.txt}.annotated_peaks.tsv"
  cp "$file" "$tsv_out"
  cp -f "$tsv_out" "$enriched_dir/"
done

################################################################################
# 8  OPTIONAL ENRICHMENT (R serial – one call per *group*, not per file)       #

# 8-A  first add ENTREZ IDs once per peak file (you already do this earlier)

# 8-B  now call the clustering script exactly once,
#      passing the directory and the mapping file
echo "[INFO] Previewing mapping file:"
head "$MAPPING_FILE"

n=$(awk 'END{print NR}' "$MAPPING_FILE")
[[ $n -le 1 ]] && err "Mapping file is empty or header-only!"
cut -f4 "$MAPPING_FILE" | sort | uniq -c

Rscript modules/pipeline2/cluster_enrichment_updated_hg_mice.R "$enriched_dir" "$MAPPING_FILE" "$GENOME_NAME" || err "Rscript enrichment failed"

################################################################################
#                         ── 9  GENERATE ENRICHMENT HTML REPORT ──            #
################################################################################
ENRICH_CLUSTER_DIR="$OUT_ROOT/Enriched_cluster"
GP_DIR="$ENRICH_CLUSTER_DIR/gprofiler"
REPORT_HTML="$ENRICH_CLUSTER_DIR/Functional_Enrichment_Report.html"
TSV_FILE="$ENRICH_CLUSTER_DIR/PRE_POSSIBLE_OUTCOME.tsv"
REPORT_TIMESTAMP=$(date +'%Y-%m-%d %H:%M:%S')

log "📄 Creating HTML functional enrichment report → $REPORT_HTML"

# === Add this check here ===
[[ -f "$TSV_FILE" ]] || { log "⚠️ Enrichment TSV file not found → $TSV_FILE"; exit 1; }

cat <<EOF > "$REPORT_HTML"
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>Functional Enrichment Summary</title>
  <style>
  body {
    font-family: Georgia, serif;
    font-size: 12px;
    max-width: 1000px;
    margin: auto;
    padding: 20px;
    line-height: 1.5;
    color: #333;
  }
  table {
    width: 100%;
    border-collapse: collapse;
    font-size: 11px;
    table-layout: fixed;
    word-wrap: break-word;
  }
  th, td {
    border: 1px solid #ccc;
    padding: 6px;
    text-align: left;
    vertical-align: top;
  }
  th {
    background-color: #f4f4f4;
  }
  tr:nth-child(even) {
    background-color: #fafafa;
  }
  td:nth-child(3),  /* GO_Terms */
  td:nth-child(4),  /* KEGG */
  td:nth-child(5) { /* Reactome */
    max-width: 250px;
    white-space: pre-wrap;
    word-break: break-word;
    overflow-wrap: break-word;
  }
  .image-grid {
    display: flex;
    flex-wrap: wrap;
    justify-content: center;
    gap: 20px;
    margin-top: 20px;
  }
  .image-grid img {
    max-width: 45%;
    border: 1px solid #ccc;
    padding: 4px;
    background: #fff;
    display: block;
    margin: auto;
  }
  .timestamp {
    font-size: 0.9em;
    color: #777;
    margin-bottom: 20px;
  }
</style>
</head>
<body>

<h1>🧬 Functional Enrichment Summary</h1>
<p class="timestamp">Generated: $REPORT_TIMESTAMP</p>

<h2>📊 Enrichment Plots</h2>
<div class="image-grid">
EOF

# Append enrichment plots
for img in "$GP_DIR"/*.png; do
  [[ -f "$img" ]] || continue
  rel_path=$(realpath --relative-to="$(dirname "$REPORT_HTML")" "$img")
  alt_text=$(basename "$img")
  echo "<div><img src=\"$rel_path\" alt=\"$alt_text\"></div>" >> "$REPORT_HTML"
done

cat <<EOF >> "$REPORT_HTML"
</div>

<h2>📄 PRE_POSSIBLE_OUTCOME.tsv</h2>
<table>
  <tr>
    <th>ENTREZID</th><th>Gene</th><th>GO_Terms</th>
    <th>KEGG</th><th>Reactome</th><th>Origin</th><th>Interpretation</th>
  </tr>
EOF
# Append TSV rows
tail -n +2 "$TSV_FILE" | while IFS=$'\t' read -r id gene go kegg reactome origin interp; do
  echo "<tr><td>${id}</td><td>${gene}</td><td>${go}</td><td>${kegg}</td><td>${reactome}</td><td>${origin}</td><td>${interp}</td></tr>" >> "$REPORT_HTML"
done

cat <<EOF >> "$REPORT_HTML"
</table>

</body>
</html>
EOF

log "✅ Enrichment HTML written → $REPORT_HTML"


# === Optional PDF version ===
REPORT_PDF="${REPORT_HTML%.html}.pdf"
html_size_kb=$(du -k "$REPORT_HTML" | cut -f1)
max_kb=150

if [[ "$html_size_kb" -le "$max_kb" ]]; then
  log "📄 Generating full PDF report (size: ${html_size_kb} KB)"
  if command -v chromium >/dev/null 2>&1; then
    chromium --headless --disable-gpu --no-sandbox \
             --virtual-time-budget=60000 \
             --print-to-pdf="$REPORT_PDF" "file://$(realpath "$REPORT_HTML")" \
             2>/dev/null
    log "✅ PDF report written → $REPORT_PDF"
  elif command -v google-chrome >/dev/null 2>&1; then
    google-chrome --headless --disable-gpu --no-sandbox \
                  --virtual-time-budget=60000 \
                  --print-to-pdf="$REPORT_PDF" "file://$(realpath "$REPORT_HTML")" \
                  2>/dev/null
    log "✅ PDF report written → $REPORT_PDF"
  else
    log "⚠️ No Chrome/Chromium found; skipping PDF"
  fi

else
  log "⚠️ HTML too large (${html_size_kb} KB) – generating plots-only PDF"

  PLOTS_ONLY_HTML="${REPORT_HTML%.html}_plots_only.html"
  PLOTS_ONLY_PDF="${REPORT_HTML%.html}_plots_only.pdf"

  cat <<EOF > "$PLOTS_ONLY_HTML"
<!DOCTYPE html>
<html>
<head>
  <meta charset="UTF-8">
  <title>Functional Enrichment – Plots Only</title>
  <style>
    body { font-family: Georgia, serif; padding: 20px; }
    .image-grid { display: flex; flex-wrap: wrap; justify-content: center; gap: 20px; }
    .image-grid img { max-width: 45%; border: 1px solid #ccc; padding: 4px; }
  </style>
</head>
<body>
<h1>🧬 Functional Enrichment – Plots Only</h1>
<div class="image-grid">
EOF

  for img in "$GP_DIR"/*.png; do
    [[ -f "$img" ]] || continue
    rel_path=$(realpath --relative-to="$(dirname "$PLOTS_ONLY_HTML")" "$img")
    echo "<img src=\"$rel_path\" alt=\"plot\">" >> "$PLOTS_ONLY_HTML"
  done

  cat <<EOF >> "$PLOTS_ONLY_HTML"
</div>
</body>
</html>
EOF

  if command -v chromium >/dev/null 2>&1; then
    chromium --headless --disable-gpu --no-sandbox \
             --virtual-time-budget=30000 \
             --print-to-pdf="$PLOTS_ONLY_PDF" "file://$(realpath "$PLOTS_ONLY_HTML")" \
             2>/dev/null
    log "✅ Plots-only PDF written → $PLOTS_ONLY_PDF"
  elif command -v google-chrome >/dev/null 2>&1; then
    google-chrome --headless --disable-gpu --no-sandbox \
                  --virtual-time-budget=30000 \
                  --print-to-pdf="$PLOTS_ONLY_PDF" "file://$(realpath "$PLOTS_ONLY_HTML")" \
                  2>/dev/null
    log "✅ Plots-only PDF written → $PLOTS_ONLY_PDF"
  else
    log "⚠️ No Chrome/Chromium found; skipping PDF"
  fi
fi

################################################################################
#                         ── 9  MERGE JOB LOGS ──                              #
################################################################################
{
  printf -- '==================== MASTER LOG (%s) ====================\n' "$TIMESTAMP"
  find "$OUT_ROOT" -name '*.annotation_stats.log' -exec bash -c 'printf -- "--- %s ---\n" "$1"; cat "$1"' _ {} \;
} > "$SUMMARY_LOG"

################################################################################
#                         ── 10  CLEANUP INTERMEDIATE FILES ──                #
################################################################################
log "🧹 Cleaning intermediate files created by annotation tools..."

find . -maxdepth 1 -type f -regextype posix-extended \
  -regex './[0-9]+\.[0-9]+(\.tmp|\.pos|\.clean\.pos|\.gtf\.tss)' \
  -exec rm -f {} +

log "🧹 Cleanup completed."


log "Finished peak annotation + enrichment. Master summary: $SUMMARY_LOG"

