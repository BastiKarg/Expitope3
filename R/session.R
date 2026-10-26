# Shared handlers for the standard and tumor submission forms.
configure_session <- function(input, output, session, config, HLAset) {
  submit_prediction <- function(suffix = "") {
    value <- function(id, default = "") {
      x <- input[[paste0(id, suffix)]]
      if (is.null(x) || !length(x) || identical(x, "")) default else x
    }
    tryCatch({
      peptide <- value("input1")
      project <- if (suffix == "") value("input4", "default") else
        paste0("tumor-", format(Sys.time(), "%Y%m%d-%H%M%S"), "-", sample.int(999999, 1))
      cleavage <- match(value("input_pep", "epitope (default)"), c("epitope (default)",
        "constitutive proteasomal (gbm)", "immunoproteasomal (gbm)",
        "constitutive proteasomal (nn)", "immunoproteasomal (nn)"))
      model <- match(value("model", "atlas_iedb"), c("atlas_iedb", "atlas", "iedb"))
      write_prediction_job(peptide, value("hla2", "HLA-A*02:01"),
        as.numeric(value("radio_mismatches", 2)), project, cleavage, model,
        if (suffix == "") value("kmer", "normal_seq") else "normal_seq",
        value("input6"), config$queue_dir)
      message <- paste("Job queued for project", project,
                       "— processing requires the Linux prediction worker.")
      output[[paste0("results", suffix)]] <- renderText(message)
    }, error = function(e) {
      output[[paste0("results", suffix)]] <- renderText(conditionMessage(e))
      showNotification(conditionMessage(e), type = "error")
    })
  }
  list(submit_prediction = submit_prediction)
}
