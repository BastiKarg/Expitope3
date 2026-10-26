source("R/runtime.R")
library(tidyr)
old_data <- function(df, tissues, limits) {
  narrow <- as.data.frame(cbind(df[, c(1, 2, 4, 28)], df[, 29:ncol(df)]))
  unused <- gather(df, tissue, expression, 2:ncol(df))
  data <- gather(narrow, tissue, expression, 5:ncol(narrow))
  data <- data[which(limits[1] < data$expression & data$expression < limits[2]), ]
  names(data)[1:4] <- c("epitope", "distance", "bind", "gene_id")
  data$info <- paste0(data$epitope, "\ndistance: ", data$distance, "\nbind: ", data$bind,
                      "\n", data$gene_id, "\n", data$tissue, ": ", data$expression)
  data[data$tissue %in% tissues, , drop = FALSE]
}
path <- "results/publication_clinical_trials_test/test_LLGPLLVLQA_EX_CCLE.tsv"
stopifnot(file.exists(path))
df <- read.csv(path, sep = "\t")
tissues <- names(df)[29:30]
for (selection in list(tissues, rev(tissues), names(df)[29:ncol(df)])) {
  before <- old_data(df, selection, c(25, 1000))
  after <- beeswarm_data(df, selection, c(25, 1000))
  rownames(before) <- rownames(after) <- NULL
  stopifnot(identical(before, after))
}
elapsed <- function(fun) median(replicate(3, system.time(fun(df, tissues, c(25, 1000)))[["elapsed"]]))
cat("Rows:", nrow(df), "tissues:", ncol(df)-28, "selected:", length(tissues), "\n")
cat("Preparation median seconds: old", elapsed(old_data), "new", elapsed(beeswarm_data), "\n")
parse("R/server.R") |> invisible()
cat("PASS: identical values, types, row order and hover/click keys\n")
