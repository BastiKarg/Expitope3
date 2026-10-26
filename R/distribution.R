# Summarise before sending data to Plotly; retain the historical joined-row counts.
expression_distribution <- function(reference, matches, tissue) {
  ref <- data.table::as.data.table(reference)
  hits <- unique(data.table::as.data.table(matches[, c(1, 2, 4, 28), drop = FALSE]))
  data.table::setnames(hits, c("epitope", "distance", "bind", "gene_id"))
  ref <- data.table::copy(ref)
  data.table::setnames(ref, c("gene_id", "expression"))
  joined <- merge(ref, hits, by = "gene_id", all = TRUE, allow.cartesian = TRUE, sort = FALSE)
  joined[is.na(expression), expression := 0]
  if (!nrow(joined) || any(!is.finite(joined$expression) | joined$expression < 0))
    stop("Reference expression must be finite and non-negative")
  joined[, value := log2(expression + 1)]
  breaks <- seq(0, max(1, max(joined$value)), length.out = 100)
  joined[, bin := pmin(99L, findInterval(value, breaks, all.inside = TRUE))]
  # Original colour precedence: non overrides weak, which overrides strong.
  genes <- unique(hits[, .(gene_id, bind)])
  strong <- genes[bind == "strong", gene_id]
  weak <- genes[bind == "weak", gene_id]
  non <- genes[bind == "non", gene_id]
  joined[, category := data.table::fcase(gene_id %in% strong, 1L,
      gene_id %in% weak, 2L, gene_id %in% non, 3L, default = 0L)]
  joined[, combined := paste(epitope, collapse = ", "), by = .(gene_id, bind)]
  counts <- joined[, .(count = .N, category = max(category)), by = bin]
  hover <- joined[!is.na(bind), .(info = paste(head(unique(paste0(gene_id,
      ": ", expression, "; mm: ", distance, "; bind: ", bind, "; seq: ", combined)), 54),
      collapse = "<br>")), by = bin]
  counts <- merge(counts, hover, by = "bin", all.x = TRUE)
  counts[is.na(info), info := "No epitope matches"]
  counts[, highlight := c("without NE", "strong", "weak", "non")[category + 1L]]
  counts[, center := (breaks[bin] + breaks[bin + 1L])/2]
  list(data = counts, width = diff(breaks)[1], tissue = tissue)
}

distribution_plot <- function(summary) {
  plotly::plot_ly(summary$data) |>
    plotly::add_bars(x = ~center, y = ~count, color = ~highlight,
      width = summary$width, text = ~info,
      hovertemplate = "log2(expression + 1): %{x}<br>Count: %{y}<br>%{text}<extra></extra>") |>
    plotly::layout(title = paste("Gene Expression Distribution in", summary$tissue),
      xaxis = list(title = "log2(gene expression + 1)"),
      yaxis = list(title = "Number of matched/reference rows", type = "log"), barmode = "stack")
}
