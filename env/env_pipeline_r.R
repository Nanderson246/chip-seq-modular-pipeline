# R Script to install Bioconductor packages for the ChIP-seq pipeline

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "AnnotationDbi",
  "clusterProfiler",
  "enrichplot",
  "ReactomePA",
  "reactome.db",
  "biomaRt",
  "org.Mm.eg.db",   # for mouse
  "org.Hs.eg.db"    # for human
))
