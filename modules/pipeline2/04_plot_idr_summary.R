#!/usr/bin/env Rscript

# Load packages
pkgs <- c("ggplot2", "readr", "dplyr", "tidyr", "stringr")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg, repos = "https://cloud.r-project.org")
}
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript plot_idr_summary.R idr_summary.tsv [source_prefix] [unique_groups.txt] [out_dir]")
}

summary_file <- args[1]
source_prefix <- ifelse(length(args) >= 2, args[2], "IDR") # Default prefix if not provided
unique_groups <- if (length(args) >= 3) args[3] else NULL
out_dir <- if (length(args) >= 4) args[4] else "analysis/IDR_Results"  # Default fallback path

# Load unique_groups.txt if provided
if (!is.null(unique_groups)) {
  unique_groups_df <- read_tsv(unique_groups, col_names = FALSE)
}

dt <- read_tsv(args[1], col_types = cols(.default = "c")) %>%
  mutate(across(where(is.character), ~ str_replace_all(., "[\r\n]+", " "))) %>%
  filter(!is.na(Pct) & Pct != "")
# 
# dt <- dt %>%
#   mutate(Pair = str_replace_all(Pair, "\\s*\\n\\s*", " vs "))

dt <- dt %>%
  mutate(
    Percent = as.numeric(str_remove(Pct, "%")),
    Total = as.numeric(Total),
    Pass = as.numeric(Pass),
    Group = ifelse(is.na(Group) | Group =="", "Unassigned", Group)
  )%>%
  filter(!is.na(Percent))
  
dt <- dt %>%
  mutate(Group = case_when(
    Group == "Unassigned" ~ sapply(Pair, function(pair) {
      pair_words <- str_split(pair, "_", simplify = TRUE)
      
      match_found <- "Unassigned"
      for (i in seq_len(nrow(unique_groups_df))) {
        group_string <- unique_groups_df[[1]][i]
        group_words <- str_split(group_string, "_", simplify = TRUE)
        
        if (all(group_words %in% pair_words)) {
          match_found <- group_string
          break
        }
      }
      match_found
    }),
    TRUE ~ Group
  ))



# Start from the original Pair
dt <- dt %>%
  mutate(Pair_wrapped = Pair)

if (!is.null(unique_groups_df) && nrow(unique_groups_df) > 0) {
  group_tokens <- unique_groups_df[[1]] %>%
    str_split("_") %>%
    unlist() %>%
    unique()
  # Remove every token (like "G4", "BLMneg", "HiSeq_3000") from Pair_wrapped
  for (token in group_tokens) {
    dt$Pair_wrapped <- str_remove_all(dt$Pair_wrapped, fixed(token))
  }
}

# Clean up underscores and wrap
dt <- dt %>%
  mutate(
    Pair_wrapped = str_replace_all(Pair_wrapped, "_{2,}", "_"),
    Pair_wrapped = str_replace_all(Pair_wrapped, "^_|_$", ""),
    Pair_wrapped = str_wrap(Pair_wrapped, width = 40),
    Pair_short = paste0("IDR_", seq_len(n()))
  )


param_df <- str_remove_all(dt$Params, "\\[|\\]")%>%
  str_split_fixed(" ", 4) %>%
  as.data.frame(stringsAsFactors = FALSE)

colnames(param_df) <- c("mu", "sigma", "rho", "p")
param_df <- param_df %>% mutate(across(everything(), as.numeric))

dt <- bind_cols(dt, param_df)

# Add IDR PASS status (e.g. > 80% pass threshold)
dt <- dt %>%
  mutate(
    IDR_Pass = if_else(Percent >= 80, "PASS", "FAIL")
  )

# Add model quality indicator
dt <- dt %>%
  mutate(
    Quality_Indicator = case_when(
      rho >= 0.8 & p >= 0.8 & p <= 0.99 & sigma >= 0.3 & sigma <= 2 & mu >= 0.3 & mu <= 3 ~ "PASS",
      TRUE ~ "CHECK"
    )
  )


# Reorder columns to place new columns after Group
# Reorder columns to place new columns after Group
group_idx <- which(colnames(dt) == "Group")
dt <- dt %>%
  select(all_of(seq_len(group_idx)), IDR_Pass, Quality_Indicator, everything())

# Optional: move IDR_Type right after Group
idr_type_idx <- which(colnames(dt) == "IDR_Type")
if (length(idr_type_idx) > 0 && idr_type_idx > group_idx) {
  dt <- dt %>%
    relocate(IDR_Type, .after = Group)
}
dt <- dt %>%
  arrange(Group, desc(Percent)) %>%
  mutate(
    Pair_short = factor(Pair_short, levels = Pair_short),
    Pair_wrapped = factor(Pair_wrapped, levels = Pair_wrapped)
  )

write_tsv(dt, file.path(out_dir, paste0(source_prefix, "_idr_summary_full.tsv")))
write.csv(dt, file.path(out_dir, paste0(source_prefix, "_idr_summary_full.csv")))

# Plot 1: IDR Pass Rate by Replicate Pair
plot1 <- ggplot(dt, aes(x = Percent, y =Pair_short , fill = Group)) +
  geom_col() +
  coord_flip() +
  labs(title = "IDR Pass Rate by Replicate Pair", y = "% Peaks Passing", x = "Replicate Pair") +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.3))

#ggsave(paste0("analysis/IDR_Results/", source_prefix, "_idr_pass_rate_barplot_v2.png"), plot = plot1, width = 12, height = 6)
ggsave(file.path(out_dir, paste0(source_prefix, "_idr_pass_rate_barplot_v2.png")), plot = plot1, width = 12, height = 6)


plot1_2 <-ggplot(dt, aes(x = Pair_wrapped, y = Percent, fill = Group))+
  geom_col() +
  labs(title = "IDR Pass Rate by Replicate Pair", y = "% Peaks Passing", x = "Replicate Pair") +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.line = element_line(color = "black", linewidth = 0.3))

#ggsave(paste0("analysis/IDR_Results/", source_prefix, "_idr_pass_rate_barplot.png"), plot = plot1_2, width = 8, height = 8)
ggsave(file.path(out_dir, paste0(source_prefix, "_idr_pass_rate_barplot.png")), plot = plot1_2, width = 8, height = 8)
# Plot 2
melted <- dt %>%
  select(mu, sigma, rho, p) %>%
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value")

plot2 <- ggplot(melted, aes(x = Parameter, y = Value, fill = Parameter)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_jitter(width = 0.1, alpha = 0.5, size = 1.2) +
  theme_minimal(base_size = 12) +
  labs(title = "Distribution of Fitted IDR Parameters") +
  theme(panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

#ggsave(paste0("analysis/IDR_Results/", source_prefix, "_idr_params_violin.png"), plot = plot2, width = 8, height = 5)
ggsave(file.path(out_dir, paste0(source_prefix, "_idr_params_violin.png")), plot =  plot2, width = 8, height = 5)

# --- Optional: Extra plots by IDR_Type (if available) ---
if ("IDR_Type" %in% colnames(dt)) {
  idr_types <- unique(dt$IDR_Type)
  
  for (idr_t in idr_types) {
    dt_sub <- dt %>% filter(IDR_Type == idr_t)
    
    # Barplot by Pair_short
    plot_bt1 <- ggplot(dt_sub, aes(x = Percent, y = Pair_short, fill = Group)) +
      geom_col() +
      coord_flip() +
      labs(title = paste("IDR Pass Rate (", idr_t, ")", sep = ""), y = "% Peaks Passing", x = "Replicate Pair") +
      theme_minimal(base_size = 7) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_line(color = "black", linewidth = 0.3))
    
    ggsave(file.path(out_dir, paste0(source_prefix, "_idr_pass_rate_barplot_bytype_", idr_t, ".png")),
           plot = plot_bt1, width = 12, height = 6)
    
    # Barplot by wrapped name
    plot_bt2 <- ggplot(dt_sub, aes(x = Pair_wrapped, y = Percent, fill = Group)) +
      geom_col() +
      labs(title = paste("IDR Pass Rate (", idr_t, ")", sep = ""), y = "% Peaks Passing", x = "Replicate Pair") +
      theme_minimal(base_size = 7) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
            axis.line = element_line(color = "black", linewidth = 0.3))
    
    ggsave(file.path(out_dir, paste0(source_prefix, "_idr_pass_rate_barplot_wrapped_bytype_", idr_t, ".png")),
           plot = plot_bt2, width = 8, height = 8)
  }
}


