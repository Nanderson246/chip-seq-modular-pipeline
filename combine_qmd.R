# combine_qmd.R

# Set your directory containing .qmd files
qmd_dir <- "content"
output_file <- "wiki_combined.Rmd"

# Files in desired order (from _quarto.yml)
files_order <- c(
  "index.qmd", "overview.qmd", "installation.qmd", "reference_setup.qmd",
  "pipeline1.qmd", "pipeline2.qmd", "mapping.qmd", "docker_pipeline1.qmd",
  "docker_pipeline2.qmd", "docker_guide.qmd", "Conda.qmd", "software.qmd",
  "troubleshooting.qmd", "citations.qmd", "license.qmd"
)

# Open output file for writing
con <- file(output_file, "w")

for (i in seq_along(files_order)) {
  file_path <- file.path(qmd_dir, files_order[i])
  lines <- readLines(file_path, warn = FALSE)
  
  # Remove YAML header (only keep for first file)
  if (i > 1) {
    in_yaml <- FALSE
    lines <- lines[!grepl("^---", lines)]  # crude YAML removal
  }
  
  # Add a heading for each section
  section_name <- tools::file_path_sans_ext(basename(files_order[i]))
  writeLines(paste0("\n\n# ", gsub("_", " ", section_name), "\n"), con)
  writeLines(lines, con)
}

close(con)
cat("✅ Combined all .qmd files into", output_file, "\n")

