args <- commandArgs(trailingOnly = TRUE)
port <- if (length(args)) as.integer(args[[1]]) else 3838L
if (dir.exists(".r-library")) .libPaths(c(normalizePath(".r-library"), .libPaths()))
shiny::runApp(".", host = "127.0.0.1", port = port, launch.browser = FALSE)
