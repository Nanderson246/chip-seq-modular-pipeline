#!/usr/bin/env Rscript

# === Install if needed ===
install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
bioc_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

install_if_missing("ggplot2")
install_if_missing("data.table")
install_if_missing("ggrepel")
install_if_missing("viridisLite")
install_if_missing("viridis")

# === Load packages ===
suppressMessages({
  library(ggplot2)
  library(data.table)
  library(ggrepel)
  library(viridis)
  library(viridisLite)
})

# === Parse input ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  cat("Usage: Rscript plot_spike_qc_summary.R <mapping_scaled.tsv> <output_dir>\n")
  quit("no", 1)
}

input_file <- args[1]
output_dir <- args[2]
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# === Load data ===
data <- fread(input_file)

# === Sanity check ===
required_cols <- c("Sample_ID", "SpikeReads", "HostReads", "ScalingFactor", "Condition")
missing_cols <- setdiff(required_cols, names(data))
if (length(missing_cols) > 0) {
  stop(paste("❌ Missing required columns:", paste(missing_cols, collapse = ", ")))
}

# Ensure required columns are numeric
data[, SpikeReads := as.numeric(SpikeReads)]
data[, HostReads := as.numeric(HostReads)]
data[, ScalingFactor := as.numeric(ScalingFactor)]

# ADD THIS BLOCK
if (all(is.na(data$SpikeReads) | data$SpikeReads == 0)) {
  cat("⚠️ No valid spike-in reads found — skipping plot generation.\n")
  file.create(file.path(output_dir, ".no_spike_plot"))
  quit("no", 0)
}

# Optional: Handle missing Instrument column
if (!"Instrument" %in% names(data)) {
  data[, Instrument := "Unknown"]
}

# Create grouping column based on Condition and Instrument
data[, Group := if ("Instrument" %in% colnames(data)) paste(Condition, Instrument, sep = "_") else Condition]

# Reorder Sample_ID by Group and SpikeReads (for consistent visual grouping)
data[, Sample_ID := factor(Sample_ID, levels = data[order(Group, -SpikeReads)]$Sample_ID)]
#Create legend
data[, GroupLegend := interaction(Condition, Instrument, sep = " / ")]


# === Plot 1: Spike-in reads barplot (colored by Condition) ===
p1 <- ggplot(data, aes(x = Sample_ID, y = SpikeReads, fill = GroupLegend)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = c(5000, 50000), linetype = "dashed", color = "red") +
  theme_minimal(base_size = 12) +
  coord_flip() +
  labs(title = "Spike-in Reads per Sample (Grouped by Condition/Instrument)",
       y = "Spike Reads", x = "Sample", fill = "Condition / Instrument") +
  theme(axis.text.y = element_text(size = 8))


# === Plot 2: Scaling factors (colored by Condition) ===
p2 <- ggplot(data, aes(x = Sample_ID, y = ScalingFactor, fill = GroupLegend)) +
  geom_bar(stat = "identity") +
  geom_hline(yintercept = c(10, 100), linetype = "dotted", color = "gray30") +
  coord_flip() +
  theme_minimal(base_size = 12) +
  labs(title = "Scaling Factor per Sample (Grouped by Condition/Instrument)",
       y = "Scaling Factor", x = "Sample", fill = "Condition / Instrument") +
  theme(axis.text.y = element_text(size = 8))

# === Plot 3: Scatter plot with Condition color and Instrument shape ===
p3 <- ggplot(data, aes(x = SpikeReads, y = HostReads, color = Condition, shape = Instrument, label = Sample_ID)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  scale_x_log10() + scale_y_log10() +
  theme_bw() +
  labs(title = "Spike vs Host Reads",
       x = "Spike Reads (log10)", y = "Host Reads (log10)",
       color = "Condition", shape = "Instrument")

# === Save outputs ===
ggsave(file.path(output_dir, "spike_reads_barplot.pdf"), p1, width = 8, height = 6)
ggsave(file.path(output_dir, "scaling_factor_barplot.pdf"), p2, width = 8, height = 6)
ggsave(file.path(output_dir, "spike_vs_host_scatter.pdf"), p3, width = 7, height = 6)

cat("✅ Plots saved to:", output_dir, "\n")
