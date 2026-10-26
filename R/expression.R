# All display variants share the same numeric data and transformation.
expression_plot <- function(data, kind, scale, ordering, tissue_groups, title) {
  stopifnot(all(c("tissue", "expression") %in% names(data)))
  if (!nrow(data) || any(!is.finite(data$expression)) || any(data$expression < 0))
    stop("Expression values must be finite and non-negative")
  medians <- tapply(data$expression, data$tissue, median)
  data$tissue <- factor(data$tissue, levels = names(sort(medians, decreasing = TRUE)))
  data$value <- if (scale == "log2") log2(data$expression + 1) else data$expression
  if (ordering == "tissue") {
    data$group <- tissue_groups$body[match(as.character(data$tissue), tissue_groups$tissue)]
    data$group[is.na(data$group)] <- "Unmapped tissue"
  } else if (ordering == "tcga_gtex") {
    data$group <- sub("[-_].*$", "", as.character(data$tissue))
  }
  layer <- switch(kind, boxplot = ggplot2::geom_boxplot(),
                  violin = ggplot2::geom_violin(),
                  beeswarm = ggbeeswarm::geom_beeswarm(), stop("Unknown plot type"))
  p <- ggplot2::ggplot(data, ggplot2::aes(x = tissue, y = value, colour = tissue)) + layer +
    ggplot2::labs(title = title, x = "Study / tissue",
                  y = if (scale == "log2") "log2(expression + 1)" else "Expression") +
    ggplot2::theme_bw() + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                                       legend.position = "none")
  if ("group" %in% names(data)) p <- p + ggplot2::facet_wrap(~group, scales = "free_x")
  plotly::ggplotly(p)
}
