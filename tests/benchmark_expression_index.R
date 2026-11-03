source('R/runtime.R')
legacy <- read_gene_expression
environment(legacy) <- list2env(list(indexed_expression_row=function(...) NULL), parent=environment(read_gene_expression))
for (tissue in c('GTEX_liver', 'GTEX_brain', 'TCGA-LIHC_tpm')) {
 path <- file.path('tcga_read_counts', paste0(tissue,'.tsv'))
 gene <- 'ENSG00000000003'
 old_time <- system.time(expected <- legacy(path,gene))['elapsed']
 first_time <- system.time(actual <- read_gene_expression(path,gene))['elapsed']
 stopifnot(identical(expected,actual))
 warm_time <- system.time(for (i in 1:20) read_gene_expression(path,gene))['elapsed']/20
 cat(tissue,'old=',old_time,'first indexed=',first_time,'warm=',warm_time,'seconds\n')
}
