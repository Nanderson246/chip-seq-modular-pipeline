#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  cat("Usage: Rscript 03_2_homer_to_peakFormat.R <peaks.txt> <base_output_path> <basename> <style>\n")
  quit("no", 1)
}


infile  <- args[1] #peak_file="$PEAK_DIR/peaks/${base}_vs_${Input_base}_peaks.txt"
outbase <- args[2] #"$PEAK_DIR/narrow_PEAKS/${base}_vs_${Input_base}"
prefix  <- args[3]
style   <- tolower(args[4])

# bed_out <- paste0(outbase, ".bed")
# npk_out <- paste0(outbase, ".narrowPeak")
# xml_out <- paste0(outbase, ".xml")

# Decide output extension based on HOMER style
if (style %in% c("histone", "super", "superhistone", "groseq")) {
  peak_ext <- "broadPeak"
} else if (style %in% c("factor", "dnase", "tss", "clip")) {
  peak_ext <- "narrowPeak"
} else if (style %in% c("mc")) {
  peak_ext <- "bed"
} else {
  peak_ext <- "narrowPeak"  # fallback
}


# Output file names (matching MACS3)
peak_out  <- paste0(outbase, "_peaks.", peak_ext)
xls_out  <- paste0(outbase, "_peaks.xls")
bed_out  <- paste0(outbase, "_summits.bed")
xml_out  <- paste0(outbase, "_peaks.xml")  # HOMER-only

# Read HOMER peak file (skip comments, no fread/readr)
lines <- readLines(infile)
lines <- lines[!grepl("^#", lines)]
tmpfile <- tempfile()
writeLines(lines, tmpfile)

# Read fixed-column HOMER peak data
df <- read.delim(tmpfile, header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("PeakID", "chr", "start", "end", "strand",
                               "norm_tag", "focus_ratio", "score", "total_tags",
                               "ctrl_tags", "fc_ctrl", "pval_ctrl", "fc_local", "pval_local", "clonal_fc"))
unlink(tmpfile)

# Replace zero p-values
df$pval_ctrl[df$pval_ctrl == 0] <- .Machine$double.xmin
df$pval_local[df$pval_local == 0] <- .Machine$double.xmin

# Compute q-values
df$qval <- p.adjust(df$pval_ctrl, method = "BH")

# Compute log10 p-value/q-value
safe_log10 <- function(x) ifelse(x > 0 & !is.na(x), -log10(x), 300)
df$log10_pval <- safe_log10(df$pval_ctrl)
df$log10_qval <- safe_log10(df$qval)

# Compute summit
df$summit <- as.integer((df$start + df$end) / 2) - df$start

# Generate peak names
df$name <- paste0(prefix, "_peak_", seq_len(nrow(df)))

# Create BED file
bed <- df[, c("chr", "start", "end", "name")]
write.table(bed, file = bed_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# === Write main peak output (narrowPeak / broadPeak / bed) ===
peak_df <- data.frame(
  chr        = df$chr,
  start      = as.integer(df$start),
  end        = as.integer(df$end),
  name       = df$name,
  score      = 1000,
  strand     = ifelse(df$strand %in% c("+", "-"), df$strand, "."),
  signalVal  = df$fc_ctrl,
  pValue     = df$log10_pval,
  qValue     = df$log10_qval,
  summit     = df$summit
)
peak_df <- peak_df[order(peak_df$chr, peak_df$start), ]
write.table(peak_df, file = peak_out, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# === Create XML version ===
xml_lines <- c('<?xml version="1.0"?>', '<peaks>')
for (i in seq_len(nrow(df))) {
  line <- sprintf('<peak id="%s" chr="%s" start="%d" end="%d" signal="%.2f" pval="%.2f" qval="%.2f" />',
                  df$name[i], df$chr[i], df$start[i], df$end[i],
                  df$fc_ctrl[i], df$log10_pval[i], df$log10_qval[i])
  xml_lines <- c(xml_lines, line)
}
xml_lines <- c(xml_lines, '</peaks>')
writeLines(xml_lines, xml_out)



# Create XLS-like summary file (tab-delimited)
xls_out <- paste0(outbase, "_peaks.xls")
summary_cols <- c("chr", "start", "end", "strand", "name", "norm_tag", "focus_ratio",
                  "score", "total_tags", "ctrl_tags", "fc_ctrl", "pval_ctrl",
                  "fc_local", "pval_local", "clonal_fc", "qval", "log10_pval", "log10_qval")

write.table(df[, summary_cols], file = xls_out, sep = "\t", quote = FALSE,
            row.names = FALSE, col.names = TRUE)
