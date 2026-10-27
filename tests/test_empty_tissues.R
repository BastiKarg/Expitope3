app <- source('R/bootstrap.R', local=TRUE)$value
shiny::testServer(app$serverFuncSource(), {
  for (empty in list(NULL, character())) {
    session$setInputs(toggle=FALSE, tcga_tissues=empty, tissue_select_plot_server=empty, tissue_select_plot=empty)
    for (id in c('tcga_boxplot', 'tcga_boxplot_with_new_features', 'beeswarm_plot_gg_server', 'beeswarm_plot_gg')) {
      err <- tryCatch(output[[id]], error=identity)
      stopifnot(inherits(err, 'shiny.silent.error'), identical(conditionMessage(err), 'Please select at least one tissue.'))
    }
  }
})
cat('Empty tissue hints passed for all four plot outputs\n')

