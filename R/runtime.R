# Shared-filesystem contract: Windows Shiny and Linux workers see the same files.
beeswarm_data <- function(df, tissues, limits) {
  # Keep source-column order so point order and swarm placement stay unchanged.
  selected <- intersect(names(df)[seq.int(29L, ncol(df))], tissues)
  narrow <- df[, c(names(df)[c(1, 2, 4, 28)], selected), drop = FALSE]
  data <- tidyr::gather(narrow, tissue, expression, tidyselect::all_of(selected))
  data <- data[which(limits[1] < data$expression & data$expression < limits[2]), , drop = FALSE]
  names(data)[1:4] <- c("epitope", "distance", "bind", "gene_id")
  data$info <- paste0(data$epitope, "\ndistance: ", data$distance,
                      "\nbind: ", data$bind, "\n", data$gene_id, "\n",
                      data$tissue, ": ", data$expression)
  data
}

queue_listing <- function(folders) {
  files <- unlist(lapply(folders, list.files), use.names = FALSE)
  if (!length(files)) return("No jobs in queue")
  paste(paste(files, "— queued (worker status unknown)"), collapse = "\n")
}

read_structure_scores <- function(path) {
  if (!file.exists(path)) stop("Missing structure score file: ", basename(path))
  lines <- readLines(path, warn = FALSE)
  lines <- sub("[\t ]+$", "", lines)
  lines <- lines[nzchar(lines)]
  if (!length(lines)) stop("Structure score file is empty")
  data.table::fread(text = paste(lines, collapse = "\n"), data.table = FALSE,
                   check.names = FALSE)
}

pdb_results_summary <- function(structure_dir) {
  files <- list_pdb_files(structure_dir)
  if (!length(files)) return("No PDB results available")
  paste(length(files), "PDB files available. Select a file below to view or download it.")
}

expitope_config <- function(root = getwd()) {
  root <- normalizePath(root, winslash = "/", mustWork = TRUE)
  setting <- function(name, default) Sys.getenv(name, unset = default)
  list(root = root,
       data_dir = setting("EXPITOPE_DATA_DIR", root),
       expression_dir = setting("EXPITOPE_EXPRESSION_DIR", if (dir.exists(file.path(root, "tcga"))) file.path(root, "tcga") else file.path(root, "tcga_read_counts")),
       structure_dir = setting("EXPITOPE_STRUCTURE_DIR", file.path(root, "3d/testing_2")),
       queue_dir = setting("EXPITOPE_QUEUE_DIR", file.path(root, "running")),
       results_dir = setting("EXPITOPE_RESULTS_DIR", file.path(root, "results")))
}

read_optional_table <- function(path, columns, sep = "\t", header = TRUE) {
  if (!file.exists(path)) {
    out <- setNames(rep(list(character()), length(columns)), columns)
    return(as.data.frame(out, stringsAsFactors = FALSE))
  }
  out <- read.table(path, sep = sep, header = header, quote = "", comment.char = "",
                    check.names = FALSE, stringsAsFactors = FALSE)
  if (!header && ncol(out) == length(columns)) names(out) <- columns
  if (!all(columns %in% names(out)))
    stop("Unexpected columns in ", path, "; required: ", paste(columns, collapse = ", "))
  out
}

list_pdb_files <- function(structure_dir) {
  folders <- file.path(structure_dir, c("best_pdbs", "mutated"))
  unique(unlist(lapply(folders, list.files, pattern = "\\.pdb$", ignore.case = TRUE), use.names = FALSE))
}

resolve_pdb_file <- function(name, structure_dir) {
  if (length(name) != 1 || is.na(name) || basename(name) != name)
    stop("Select a PDB filename from the list")
  candidates <- file.path(structure_dir, c("best_pdbs", "mutated", "upload"), name)
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("PDB file not found: ", name)
  found[[1]]
}

pdb_residue_selection <- function(path, chain = "B") {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[substr(lines, 1, 6) == "ATOM  " & substr(lines, 22, 22) == chain]
  residues <- unique(trimws(substr(lines, 23, 27)))
  residues <- residues[nzchar(residues)]
  if (!length(residues)) return("none")
  # Explicit residues preserve gaps and insertion codes; avoid a misleading range.
  paste(paste0(residues, ":", chain), collapse = " or ")
}

read_gene_expression <- function(path, gene) {
  if (!file.exists(path)) stop("Missing expression table: ", path)
  tab <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (grepl("_tpm\\.tsv$", path)) {
    names(tab)[1] <- "id"
    tab$id <- sub("[._].*$", "", tab$id)
  }
  if (!ncol(tab) || !"id" %in% names(tab)) stop("Expression table requires an id column: ", path)
  rows <- tab[tab$id == gene, , drop = FALSE]
  if (nrow(rows) != 1) stop("Expected one row for ", gene, " in ", path, "; found ", nrow(rows))
  samples <- setdiff(names(rows), "id")
  data.frame(ENSG = gene, tissue = sub("_tpm$", "", tools::file_path_sans_ext(basename(path))),
             case_id = samples, expression = as.numeric(unlist(rows[1, samples, drop = FALSE], use.names = FALSE)), check.names = FALSE)
}

expression_files <- function(folder) {
  tpm <- list.files(folder, pattern = "_tpm\\.tsv$", full.names = TRUE)
  files <- if (length(tpm)) c(tpm, list.files(folder, pattern = "^GTEX.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)) else
    list.files(folder, pattern = "\\.tsv$", full.names = TRUE)
  files <- unique(files)
  files <- files[!grepl("\\.samples\\.tsv$", files)]
  setNames(files, sub("_tpm$", "", tools::file_path_sans_ext(basename(files))))
}

expression_gene_metadata <- function(files) {
  if (!length(files)) return(data.frame(ENSG = character(), gene = character()))
  tpm <- files[grepl("_tpm\\.tsv$", files)]
  path <- unname(if (length(tpm)) tpm[[1]] else files[[1]])
  ids <- as.character(data.table::fread(path, select = 1L, data.table = FALSE)[[1]])
  if (!grepl("_tpm\\.tsv$", path)) return(data.frame(ENSG = ids, gene = ids))
  data.frame(ENSG = sub("[._].*$", "", ids),
             gene = sub("^[^_]+_([^_]+)_.*$", "\\1", ids), stringsAsFactors = FALSE)
}

write_prediction_job <- function(peptide, hla, mismatches, project, cleavage, model,
                                 mode = "normal_seq", custom = "", queue_dir) {
  peptide <- trimws(peptide); custom <- toupper(trimws(custom)); hla <- trimws(hla)
  if (length(peptide) != 1 || !grepl("^[ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy]{8,}$", peptide))
    stop("Enter a peptide of at least 8 standard amino acids")
  if (!grepl("^HLA-[ABC]\\*[0-9]{2,3}:[0-9]{2,3}(:[0-9]{2,3})*$", hla))
    stop("Select one HLA class-I allele")
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]*$", project) || project %in% c(".", ".."))
    stop("Project names may contain letters, numbers, underscores, dots and hyphens")
  integer_in <- function(x, low, high) length(x) == 1 && !is.na(x) && x == as.integer(x) && x >= low && x <= high
  if (!integer_in(mismatches, 0, 4) || !integer_in(cleavage, 1, 5) || !integer_in(model, 1, 3))
    stop("Invalid mismatch count or model selection")
  if (!mode %in% c("normal_seq", "indel", "multigap") && !grepl("^(8|9|10|11|12)$", mode))
    stop("Choose normal, indel, multigap, or an 8–12 amino-acid k-mer")
  if (nzchar(custom) && (!grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", custom) || nchar(custom) < nchar(peptide)))
    stop("The custom protein must contain standard amino acids and be at least as long as the peptide")
  dir.create(queue_dir, recursive = TRUE, showWarnings = FALSE)
  fields <- c(peptide, hla, mismatches, project, cleavage, model, mode)
  if (nzchar(custom)) fields <- c(fields, custom)
  # Stage outside the queue: the worker scans every file inside it.
  temporary <- tempfile("expitope-", tmpdir = dirname(queue_dir))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(paste(fields, collapse = " "), temporary, useBytes = TRUE)
  target <- file.path(queue_dir, paste0(project, "-", basename(temporary), ".job"))
  if (!file.rename(temporary, target)) stop("Could not atomically submit the job")
  target
}
