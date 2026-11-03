source("R/runtime.R")
local({
  folder <- tempfile("expression-index-")
  dir.create(folder)
  on.exit(unlink(folder, recursive = TRUE))
  for (newline in c("\n", "\r\n")) {
    path <- file.path(folder, "TCGA-test_tpm.tsv")
    header <- paste0("id\tsample1\tsample2", newline)
    rows <- paste0(c("ENSG1.2_GENE_coding\t0\t1.5", "ENSG2.1_OTHER_coding\t2\t3"), newline)
    writeBin(charToRaw(paste0(header, paste0(rows, collapse=""))), path)
    info <- file.info(path)
    writeLines(c(paste(info$size, sprintf("%.6f", as.numeric(info$mtime)), sep="\t"),
                 "id\toffset\tlength",
                 paste("ENSG1", nchar(header, type="bytes"), nchar(rows[1], type="bytes"), sep="\t"),
                 paste("ENSG2", nchar(paste0(header,rows[1]), type="bytes"), nchar(rows[2], type="bytes"), sep="\t")), paste0(path,".gidx"))
    stopifnot(identical(read_gene_expression(path,"ENSG1")$expression, c(0,1.5)),
              identical(read_gene_expression(path,"ENSG2")$expression, c(2,3)))
    stopifnot(inherits(tryCatch(read_gene_expression(path,"absent"),error=identity),"error"))
    # A changed source invalidates the cached index and uses the ordinary reader.
    writeBin(charToRaw(paste0(header,"ENSG1.2_GENE_coding\t10\t20",newline)), path)
    stopifnot(is.null(indexed_expression_row(path,"ENSG1")),
              identical(read_gene_expression(path,"ENSG1")$expression,c(10,20)))
  }
})
cat("PASS indexed LF/CRLF rows, gene normalization, missing genes and stale-index fallback\n")
