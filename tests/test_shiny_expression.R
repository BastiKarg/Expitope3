app <- source('app.R')$value
shiny::testServer(server, {
  session$setInputs(toggle=FALSE, rows='2', tcga_ensg='ENSG00000000003', tcga_tissues='TCGA-ACC', tcga_plot='boxplot', tcga_order='decreasing', tcga_scale='no scale')
  for (kind in c('boxplot','violin','beeswarm')) for (scale in c('no scale','log2')) for (ordering in c('decreasing','tissue','tcga_gtex')) {
    session$setInputs(tcga_plot=kind, tcga_scale=scale, tcga_order=ordering)
    value <- output$tcga_boxplot_with_new_features
    stopifnot(nchar(value)>100)
    cat('PASS',kind,scale,ordering,'\n')
  }
})
