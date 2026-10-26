source("R/distribution.R")
matches <- as.data.frame(matrix("", nrow = 3, ncol = 28))
matches[[1]] <- c("AAA", "BBB", "CCC")
matches[[2]] <- c(0, 1, 2)
matches[[4]] <- c("strong", "weak", "non")
matches[[28]] <- c("g1", "g1", "missing")
ref <- data.frame(id = c("g1", "g2"), tissue = c(4, 0))
x <- expression_distribution(ref, matches, "tissue")
stopifnot(sum(x$data$count) == 4, all(is.finite(x$data$center)), nrow(x$data) <= 99)
zero <- expression_distribution(transform(ref, tissue = 0), matches, "tissue")
stopifnot(sum(zero$data$count) == 4, nrow(zero$data) == 1)
stopifnot(length(plotly::plotly_build(distribution_plot(x))$x$data) > 0)
cat("PASS duplicate matches, missing reference genes, zero range and Plotly rendering\n")
path <- "results/publication_clinical_trials_test/test_LLGPLLVLQA_EX_PAX_v4.tsv"
if (file.exists(path)) {
  matches <- read.csv(path, sep = "\t")
  refpath <- "Datenbanken_Ordner_shiny/PAX_v4.tsv"
  headers <- names(data.table::fread(refpath, nrows = 0))
  column <- match(names(matches)[29], make.names(headers, unique = TRUE))
  ref <- data.table::fread(refpath, select = c(1L, column), data.table = FALSE)
  timing <- system.time(x <- expression_distribution(ref, matches, headers[column]))
  cat("Demo:", nrow(ref), "reference rows;", sum(x$data$count), "joined rows;",
      nrow(x$data), "bins; seconds:", timing[["elapsed"]], "\n")
}
