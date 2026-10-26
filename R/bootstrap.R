# Runtime setup for the Windows UI / shared Linux worker deployment.
if (dir.exists(".r-library")) .libPaths(c(normalizePath(".r-library"), .libPaths()))
source("R/runtime.R", local = TRUE)
config <- expitope_config()
packages <- c("shiny", "ggplot2", "plotly", "ggbeeswarm", "data.table", "tidyr",
              "magrittr", "dplyr", "stringr", "zip", "shinyBS", "Hmisc", "gdata",
              "shinyWidgets", "colourpicker", "htmlwidgets", "shinyjs",
              "shinycssloaders", "DT", "NGLVieweR", "Biostrings", "Rpdb", "ggpubr")
missing_packages <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Run scripts/install_shiny_dependencies.R; missing packages: ",
                                   paste(missing_packages, collapse = ", "))
for (package in packages) suppressPackageStartupMessages(library(package, character.only = TRUE))
HLAset <- read_optional_table(file.path(config$data_dir, "hla_list_hla_only_resorted2"), "HLA", header = FALSE)
HLAset <- data.frame(HLA = unique(trimws(unlist(strsplit(HLAset$HLA, ",", fixed = TRUE)))))
ENSGset <- read_optional_table(file.path(config$data_dir, "ensg_list"), "ensg", sep = ",", header = FALSE)
tissue_groups <- read_optional_table(file.path(config$data_dir, "tcga_gtex_tissues.tsv"), c("tissue", "body"))
gtex_groups <- read_optional_table("docs/gtex_tissue_groups.tsv", c("group", "body"))
if (nrow(gtex_groups)) {
  tissue_groups <- tissue_groups[, c("tissue", "body")]
  tissue_groups <- tissue_groups[!grepl("^GTEX_", tissue_groups$tissue), ]
  tissue_groups <- rbind(tissue_groups, unique(data.frame(tissue = gtex_groups$group, body = gtex_groups$body)))
}
ensg_genename <- read_optional_table(file.path(config$data_dir, "ensg_genename.tsv"), c("ENSG", "gene"))
expression_paths <- expression_files(config$expression_dir)
if (!nrow(ensg_genename)) ensg_genename <- expression_gene_metadata(expression_paths)
if (!nrow(ENSGset)) ENSGset <- data.frame(ensg = unique(ensg_genename$ENSG))
# Supplied test.tsv contains test1..test4 placeholders rather than a TCRdock schema.
tcr_dock <- data.frame(organism = "human", mhc_class = 1L, mhc = "", peptide = "",
                       va = "", ja = "", cdr3a = "", vb = "", jb = "", cdr3b = "")
for (folder in c(config$queue_dir, config$results_dir, "recent_finish"))
  dir.create(folder, recursive = TRUE, showWarnings = FALSE)
for (folder in c("todo", "todo_mutate", "todo_alascan", "todo_mutate_rosetta", "upload"))
  dir.create(file.path(config$structure_dir, folder), recursive = TRUE, showWarnings = FALSE)
source("R/legacy_helpers.R", local = TRUE)
source("R/session.R", local = TRUE)
source("R/expression.R", local = TRUE)
source("R/distribution.R", local = TRUE)
source("R/ui.R", local = TRUE)
source("R/server.R", local = TRUE)
required_data <- c("tcga_gtex_tissues.tsv", "tcga_abbreviations")
unavailable <- required_data[!file.exists(file.path(config$data_dir, required_data))]
if (!length(expression_paths)) unavailable <- c(unavailable, "expression tables")
if (!nrow(ENSGset)) unavailable <- c(unavailable, "gene identifiers")
if (length(unavailable)) {
  ui <- tagList(tags$div(role = "status", style = "padding:12px;background:#fff3cd;color:#664d03",
    strong("Some datasets are unavailable. "), paste("Missing:", paste(unavailable, collapse = ", ")),
    " — TCGA/GTEx views require the original assets; no placeholder expression values are used."), ui)
}
legacy_server <- server
server <- function(input, output, session) {
  handlers <- configure_session(input, output, session, config, HLAset)
  session_server <- legacy_server
  scope <- new.env(parent = environment(legacy_server))
  scope$submit_prediction <- handlers$submit_prediction
  environment(session_server) <- scope
  session_server(input, output, session)
}
shinyApp(ui, server)
