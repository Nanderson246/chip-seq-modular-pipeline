#!/usr/bin/env Rscript


## === 1. Installer functions ===
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


## === 2. Install all required packages ===
pkgs_cran <- c("httr", "jsonlite", "tidyr", "data.table", "stringr",
               "dplyr", "readr", "tools", "ggplot2", "gprofiler2")

pkgs_bioc <- c("clusterProfiler", "enrichplot",
               "ReactomePA", "reactome.db", "biomaRt", "AnnotationDbi",
               "org.Hs.eg.db", "org.Mm.eg.db")

invisible(lapply(pkgs_cran, install_if_missing))
invisible(lapply(pkgs_bioc, bioc_if_missing))

#test tomorrow
# Rscript modules/pipeline2/cluster_enrichment_updated_hg_mice.R \
# analysis/ChIPseeker_TSS_Hommer_IDR_annotation/macs3/Enriched \
# analysis/ChIPseeker_TSS_Hommer_IDR_annotation/macs3/filename_mapping.tsv \
# mm10

## === 3. Load libraries dynamically based on genome ===
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: cluster_enrichment.R <enriched_dir> <filename_mapping.tsv> [genome]")
}


in_dir    <- args[1]
map_path  <- args[2]
genome    <- ifelse(length(args) >= 3, args[3], "hg38")

gost_organism <- switch(tolower(genome),
                        "hg38" = "hsapiens",
                        "mm10" = "mmusculus",
                        stop("Unsupported genome name for gprofiler: ", genome)
)


load_orgdb <- function(genome) {
  pkg <- switch(tolower(genome),
    "hg38" = "org.Hs.eg.db",
    "mm10" = "org.Mm.eg.db",
    stop("Unsupported genome name: ", genome)
  )
  suppressMessages(library(pkg, character.only = TRUE))
  get(pkg)
}

orgdb <- load_orgdb(genome)

suppressMessages({
  library(reactome.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(AnnotationDbi)
  library(ggplot2)
  library(gprofiler2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(data.table)
  library(tools)
})


## === 4. Set plot theme ===
plot_theme <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", colour = "grey80", linewidth = 0.8),
      plot.background = element_rect(fill = "white", colour = "grey40", linewidth = 0.8),
      axis.line = element_blank(),
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank()
    )
}



## === 5. Setup options ===
options(echo = TRUE)
options(warn = 1)
options(timeout = 600)

## === 6. Handle command line arguments ===
out_dir <- file.path(dirname(in_dir), "Enriched_cluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


## === 7. Read mapping file ==

filename_df <- read.delim(map_path, stringsAsFactors = FALSE)

# Make sure expected columns exist
if (!"group" %in% names(filename_df)) filename_df$group <- filename_df$biological_identity
filename_df$Group_collapsed <- gsub("[^a-zA-Z0-9]+", "_", filename_df$group)
filename_df$Group_collapsed_intersect <- filename_df$Group_collapsed
if (!"Instrument" %in% names(filename_df)) filename_df$Instrument <- "DEFAULT"
data <- filename_df

##############################################################################
## 0. helper: Ensembl → Entrez (row-wise, keeps order)              ──────────
##############################################################################

ensembl_to_entrez_vec <- function(ens_vec) {
  if (length(ens_vec) == 0) return(character(0))
  
  # Strip version if needed (optional if already clean)
  keys <- ens_vec
  
  # Skip NA-only vectors
  if (all(is.na(keys))) return(rep(NA_character_, length(keys)))
  
  # Try to run select(), handle invalid keys gracefully
  conv <- tryCatch({
    AnnotationDbi::select(
      orgdb,
      keys    = unique(na.omit(keys)),
      keytype = "ENSEMBL",
      columns = "ENTREZID"
    )
  }, error = function(e) {
    message("[WARNING] Annotation failed: ", e$message)
    return(NULL)
  })
  
  if (is.null(conv) || nrow(conv) == 0) {
    message("[INFO] No valid mappings found for this file — filling ENTREZID with NA")
    return(rep(NA_character_, length(keys)))
  }
  
  e2e <- setNames(conv$ENTREZID, conv$ENSEMBL)
  return(unname(e2e[keys]))
}


##############################################################################
## 1. function that patches ONE file                               ──────────
##############################################################################


add_entrez_column <- function(tsv_file, overwrite = TRUE) {
  cat("\n[LOAD] Reading file:", tsv_file, "\n")
  
  # Try reading the file
  df <- tryCatch({
    read.delim(tsv_file, quote = "", comment.char = "", check.names = FALSE)
  }, error = function(e) {
    message("[ERROR] Failed to read file: ", tsv_file)
    message("        ", e$message)
    return(NULL)
  })
  
  if (is.null(df)) return(invisible(FALSE))
  
  # Check column names
  names(df) <- gsub("^ensembl_id$", "Ensembl_ID", names(df), ignore.case = TRUE)
  
  if (!"Ensembl_ID" %in% names(df)) {
    message("[SKIP] No Ensembl_ID in ", basename(tsv_file))
    return(invisible(FALSE))
  }
  
  if ("ENTREZID" %in% names(df)) {
    message("[INFO] ENTREZID already exists → ", basename(tsv_file))
    return(invisible(FALSE))
  }
  
  cat("[DEBUG] Number of Ensembl_IDs in file:", length(df$Ensembl_ID), "\n")
  cat("[DEBUG] Sample Ensembl_IDs:\n")
  print(head(df$Ensembl_ID))
  
  # Do the conversion
  df$ENTREZID <- ensembl_to_entrez_vec(df$Ensembl_ID)
  
  out_path <- if (overwrite) {
    tsv_file
  } else {
    sub("\\.tsv$", "_with_entrez.tsv", tsv_file, ignore.case = TRUE)
  }
  
  write.table(df, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[OK] Saved updated file → ", out_path)
  return(TRUE)
}



##############################################################################
## 2. run over every annotated_peaks.tsv in the Enriched directory ──────────
##############################################################################
tsv_files <- list.files(in_dir, pattern = "\\.annotated_peaks\\.tsv$",
                        full.names = TRUE, recursive = FALSE)

lapply(tsv_files, add_entrez_column, overwrite = TRUE)
##############################################################################



suppressMessages({
  library(org.Hs.eg.db)      # KEGG + GO
  has_reactome <- requireNamespace("reactome.db", quietly = TRUE)
  if (has_reactome) library(reactome.db)
})

collapse_vec <- function(x) paste(unique(na.omit(x)), collapse = ";")

add_all_pathways <- function(tsv_path) {
  
  df <- read.delim(tsv_path, quote = "", comment.char = "", check.names = FALSE)
  if (!"ENTREZID" %in% names(df)) return(invisible(FALSE))
  df$ENTREZID <- as.character(df$ENTREZID)
  ids <- unique(na.omit(df$ENTREZID))
  if (length(ids) == 0) {
    message("[SKIP] No valid ENTREZID values in file: ", basename(tsv_path))
    return(invisible(FALSE))
  }
  
  
  ## ---- KEGG (PATH) --------------------------------------------------------
   kmap <- tapply(
    AnnotationDbi::select(orgdb, ids, "PATH", "ENTREZID")$PATH,
    AnnotationDbi::select(orgdb, ids, "PATH", "ENTREZID")$ENTREZID,
    collapse_vec)
  
  ## ---- GO split into BP / CC / MF ----------------------------------------
  g <- AnnotationDbi::select(orgdb, ids, c("GO","ONTOLOGY"), "ENTREZID")
  gmap_BP <- tapply(g$GO[g$ONTOLOGY=="BP"], g$ENTREZID[g$ONTOLOGY=="BP"], collapse_vec)
  gmap_CC <- tapply(g$GO[g$ONTOLOGY=="CC"], g$ENTREZID[g$ONTOLOGY=="CC"], collapse_vec)
  gmap_MF <- tapply(g$GO[g$ONTOLOGY=="MF"], g$ENTREZID[g$ONTOLOGY=="MF"], collapse_vec)
  
  ## ---- Reactome (if package present) -------------------------------------
  if (has_reactome) {
    r <- AnnotationDbi::select(reactome.db, ids, "PATHID", "ENTREZID")
    rmap <- tapply(r$PATHID, r$ENTREZID, collapse_vec)
  } else {
    rmap <- setNames(rep(NA_character_, length(ids)), ids)
  }
  
  ## ---- attach columns ----------------------------------------------------
  df$KEGG     <- kmap[df$ENTREZID]
  df$REACTOME <- rmap[df$ENTREZID]
  df$GO_BP    <- gmap_BP[df$ENTREZID]
  df$GO_CC    <- gmap_CC[df$ENTREZID]
  df$GO_MF    <- gmap_MF[df$ENTREZID]
  
  write.table(df, file = tsv_path, sep = "\t",
              quote = FALSE, row.names = FALSE)
  message("[OK] KEGG / Reactome / GO added → ", basename(tsv_path))
  TRUE
}

## run across files
tsv_files <- list.files(in_dir, pattern="\\.annotated_peaks\\.tsv$", full.names=TRUE, recursive=TRUE)
lapply(tsv_files, add_all_pathways)




get_entrez_chunked <- function(sn) {
  enriched_dir <- file.path(dirname(in_dir), "Enriched")
  pattern <- paste0("^", sn, "_.*\\.annotated_peaks\\.tsv$")
  files <- list.files(enriched_dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    cat("[WARN] No file starting with", sn, "in", enriched_dir, "\n")
    return(NULL)
  }
  if (length(files) > 1) {
    cat("[INFO] Multiple matches for", sn, "- using first match:", basename(files[1]), "\n")
  }
  
  f <- files[1]
  df <- tryCatch(read.delim(f, comment.char = "", quote = ""), error = function(e) {
    cat("[ERROR] Failed to read:", f, "\n")
    return(NULL)
  })
  if (is.null(df)) return(NULL)
  
  # Normalize column name
  colnames(df) <- gsub("^ensembl_id$", "Ensembl_ID", colnames(df), ignore.case = TRUE)
  
  ## -------------------------------------------
  ## FIRST RETURN CASE: ENTREZID column exists
  ## -------------------------------------------
  if ("ENTREZID" %in% colnames(df)) {
    idvec <- na.omit(df$ENTREZID)
    if (length(idvec)) {
      cat("[DEBUG] Found ENTREZID directly for", sn, "- N =", length(idvec), "\n")
      return(unique(as.character(idvec)))  # ensure unique + character
    }
  }
  
  ## -------------------------------------------
  ## SECOND RETURN CASE: fallback to Ensembl
  ## -------------------------------------------
  if ("Ensembl_ID" %in% colnames(df)) {
    converted <- ensembl_to_entrez_vec(df$Ensembl_ID)
    cat("[DEBUG] Converted Ensembl_IDs to ENTREZ for", sn, "- N =", length(converted), "\n")
    return(unique(as.character(converted)))
  }
  
  cat("[WARN] No usable identifiers in", basename(f), "\n")
  return(NULL)
}

  ##############################################################################
  ## Build union + intersection + private genes by group              ─────────
  ##############################################################################
  build_gene_sets <- function(subset, group_field) {
    
    # Split data by group
    split_groups <- split(subset, subset[[group_field]])
    
    # Get all ENTREZ IDs per file (by shortname) for each group
    all_union <- lapply(split_groups, function(subgroup) {
      all_genes <- lapply(subgroup$shortname, get_entrez_chunked)
      all_genes <- Filter(Negate(is.null), all_genes)
      unique(unlist(all_genes))
    })
    
    all_intersect <- lapply(split_groups, function(subgroup) {
      all_genes <- lapply(subgroup$shortname, get_entrez_chunked)
      all_genes <- Filter(Negate(is.null), all_genes)
      if (length(all_genes) > 1) Reduce(intersect, all_genes)
      else unique(unlist(all_genes))
    })
    
    # Compute private genes: unique to each group (based on union)
    grp_names <- names(all_union)
    private <- lapply(grp_names, function(g) {
      others <- unique(unlist(all_union[setdiff(grp_names, g)]))
      setdiff(all_union[[g]], others)
    }) |> `names<-`(grp_names)
    
    # Return as a list of lists
    return(list(
      union = all_union,
      intersect = all_intersect,
      private = private
    ))
  }
  


  
  
save_all_plots <- function(gost_obj, prefix, title, folder) {
  dir.create(folder, showWarnings = FALSE, recursive = TRUE)
  
  gdf <- gost_obj$result
  top_terms <- dplyr::slice_head(gdf, n = 20)
  if (!all(c("term_name","p_value","term_size") %in% colnames(top_terms))) return()
  
  ## optional – save TSV of full enrichment table
  list_cols <- vapply(gdf, is.list, logical(1))
  if (any(list_cols)) {
    gdf[list_cols] <- lapply(gdf[list_cols], function(col)
      vapply(col, \(x) paste(as.character(unlist(x)), collapse=","), ""))
  }
  write.table(gdf,
              file = file.path(folder, paste0(prefix, "_gprofiler.tsv")),
              sep  = "\t", quote = FALSE, row.names = FALSE)
  

  # Dot plot
  p1 <- ggplot(top_terms,
               aes(x = -log10(p_value),
                   y = reorder(term_name, -p_value))) +
    geom_point(aes(size = term_size, fill = p_value),
               shape = 21, stroke = 0.8) +
    scale_fill_gradient(low = "#D73027", high = "#4575B4") +
    plot_theme() +
    labs(title = paste(title, "- Dot Plot"),
         x = "-log10(p-value)", y = "Term")
  ggsave(file.path(folder, paste0(prefix, "_dotplot.png")),
         p1, width = 12, height = 8, dpi = 300)
  
  # Bar plot
  p2 <- ggplot(top_terms,
               aes(x = reorder(term_name, -p_value),
                   y = -log10(p_value))) +
    geom_col(fill = "steelblue") +
    coord_flip() + plot_theme() +
    labs(title = paste(title, "- Bar Plot"),
         x = "Term", y = "-log10(p-value)")
  ggsave(file.path(folder, paste0(prefix, "_barplot.png")),
         p2, width = 12, height = 8, dpi = 300, bg = "white")
  
  # Manhattan plot
  gp <- gostplot(gost_obj, capped = TRUE, interactive = FALSE)
  ggsave(file.path(folder, paste0(prefix, "_manhattan.png")),
         gp, width = 14, height = 8, dpi = 300, bg = "white")
}

# === Ensure gene sets are defined ===
gs_all   <- build_gene_sets(data, group_field = "Group_collapsed")
gs_inter <- gs_all$intersect
gs_coll  <- gs_all$union



# === PERFORM g:Profiler ENRICHMENT PER GROUP (UNION SETS) ==================
# gost_organism <- switch(tolower(genome),
#                         "hg38" = "hsapiens",
#                         "mm10" = "mmusculus",
#                         stop("Unsupported genome name for gprofiler: ", genome)
# )

for (grp in names(gs_coll)) {
  genes <- gs_coll[[grp]]
  if (length(genes) < 3) next  # skip small gene sets
  
  cat("[INFO] Running enrichment for group:", grp, "\n")
  gost_res <- tryCatch({
    gprofiler2::gost(genes, organism = gost_organism, significant = TRUE)
  }, error = function(e) {
    message("[ERROR] g:Profiler failed for group ", grp, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(gost_res)) {
    plot_dir <- file.path(out_dir, "gprofiler")
    save_all_plots(gost_res, prefix = grp, title = grp, folder = plot_dir)
  }
}


##############################################################################
## === PRE_POSSIBLE_OUTCOME summary (LOCAL, no biomaRt) =====================
##############################################################################

# -- 1. Gather all ENT IDs that appear in the current analysis --------------
intersect_genes_all <- unique(unlist(gs_inter))   # built earlier
union_genes_all     <- unique(unlist(gs_coll))
all_genes           <- unique(c(intersect_genes_all, union_genes_all))
if (!length(all_genes)) quit(status = 0)

# -- 2. Read every annotated_peaks.tsv to build a local annotation table ----
tsv_files <- list.files(
  in_dir,
  pattern    = "\\.annotated_peaks\\.tsv$",
  full.names = TRUE,
  recursive  = TRUE
)

annot_list <- lapply(tsv_files, \(f) {
  df <- tryCatch(
    read.delim(f, quote = "", comment.char = "", check.names = FALSE),
    error = \(e) NULL
  )
  if (is.null(df) || !"ENTREZID" %in% names(df)) return(NULL)
  # keep only rows with an ENTREZ that we actually need
  df <- df[ df$ENTREZID %in% all_genes , ]
  if (!nrow(df)) return(NULL)
  # select the columns you care about; use NA if absent
  data.frame(
    ENTREZID  = df$ENTREZID,
    GeneName  = df$`Gene Name` %||% df$external_gene_name %||% NA,
    GO_BP     = if ("GO_BP"    %in% names(df)) df$GO_BP    else NA,
    KEGG      = if ("KEGG"     %in% names(df)) df$KEGG     else NA,
    Reactome  = if ("REACTOME" %in% names(df)) df$REACTOME else NA,
    stringsAsFactors = FALSE
  )
})
local_annot <- dplyr::bind_rows(annot_list) |> dplyr::distinct()

# -- 3. Collapse multiple values per gene -----------------------------------
collapse <- \(x) paste(unique(na.omit(x)), collapse = ";")
annotation_summary <- local_annot |>
  dplyr::group_by(ENTREZID, GeneName) |>
  dplyr::summarise(
    GO_Terms = collapse(GO_BP),
    KEGG     = collapse(KEGG),
    Reactome = collapse(Reactome),
    .groups  = "drop"
  ) |>
  dplyr::mutate(
    Origin = dplyr::case_when(
      ENTREZID %in% intersect_genes_all ~ "Intersection",
      TRUE                               ~ "Union"
    ),
    Interpretation = ifelse(
      Origin == "Intersection",
      "Likely shared / constitutive",
      "Condition-specific or noise"
    )
  )

# -- 4. Save -----------------------------------------------------------------
write.table(
  annotation_summary,
  file = file.path(out_dir, "PRE_POSSIBLE_OUTCOME.tsv"),
  sep  = "\t", quote = FALSE, row.names = FALSE
)

cat("[SUMMARY] PRE_POSSIBLE_OUTCOME.tsv written with:",
    nrow(annotation_summary), "genes\n")
    
