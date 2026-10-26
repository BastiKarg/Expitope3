# Install only into this project's private library, leaving system R untouched.
args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) normalizePath(args[[1]], mustWork = TRUE) else getwd()
library_path <- file.path(project, ".r-library")
dir.create(library_path, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library_path, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 600)
packages <- c(
  "shiny", "ggplot2", "plotly", "ggbeeswarm", "data.table", "tidyr",
  "magrittr", "dplyr", "stringr", "zip", "shinyBS", "Hmisc", "gdata",
  "shinyWidgets", "colourpicker", "htmlwidgets", "shinyjs", "shinycssloaders",
  "DT", "NGLVieweR", "Rpdb", "ggpubr", "rmarkdown", "pandoc", "BiocManager"
)
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, lib = library_path,
                   type = if (.Platform$OS.type == "windows") "binary" else "source")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) stop("BiocManager installation failed")
  BiocManager::install("Biostrings", lib = library_path, ask = FALSE, update = FALSE,
                       type = if (.Platform$OS.type == "windows") "binary" else "source")
}
required <- c(packages, "Biostrings")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Packages still unavailable: ", paste(missing, collapse = ", "))
cat("Shiny dependencies installed in", library_path, "\n")
write.csv(as.data.frame(installed.packages(lib.loc = library_path))[, c("Package", "Version")],
          file.path(library_path, "installed-versions.csv"), row.names = FALSE)
