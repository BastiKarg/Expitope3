# Integration check against locally imported data; research files are not bundled.
source("R/runtime.R")
source("R/expression.R")
files <- expression_files(expitope_config()$expression_dir)
mapping <- read.delim("docs/gtex_tissue_groups.tsv")
stopifnot(all(unique(mapping$group) %in% names(files)))
stopifnot(!any(grepl("samples", names(files))))
gene <- "ENSG00000000003"
data <- rbind(read_gene_expression(files[["GTEX_liver"]], gene),
              read_gene_expression(files[["TCGA-LIHC"]], gene))
stopifnot(all(is.finite(data$expression)), all(data$expression >= 0),
          length(unique(data$tissue)) == 2)
groups <- data.frame(tissue=c("GTEX_liver", "TCGA-LIHC"), body="Liver")
plot <- expression_plot(data, "boxplot", "log2", "tissue", groups, gene)
stopifnot(inherits(plotly::plotly_build(plot), "plotly"))
cat(length(unique(mapping$group)), "GTEx groups discovered; liver TCGA/GTEx plot passed with", nrow(data), "samples\n")
