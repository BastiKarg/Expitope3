# Link Plotly keys to original rows, independent of trace/colour grouping.
link_plot_selection <- function(source, session) {
  rows <- shiny::reactiveVal(data.frame())
  selected <- shiny::reactiveVal(integer())
  shiny::observeEvent(rows(), selected(integer()), ignoreNULL = FALSE)
  indices <- function(event) {
    if (is.null(event$key)) return(integer())
    match(unique(as.character(unlist(event$key))), as.character(seq_len(nrow(rows()))), nomatch = 0L) |>
      (function(x) x[x > 0L])()
  }
  shiny::observeEvent(plotly::event_data("plotly_selected", source = source, session = session, priority = "event"), {
    selected(indices(plotly::event_data("plotly_selected", source = source, session = session, priority = "event")))
  }, ignoreNULL = FALSE, ignoreInit = TRUE)
  shiny::observeEvent(plotly::event_data("plotly_deselect", source = source, session = session, priority = "event"), {
    selected(integer())
  })
  shiny::observeEvent(plotly::event_data("plotly_click", source = source, session = session, priority = "event"), {
    ids <- indices(plotly::event_data("plotly_click", source = source, session = session, priority = "event"))
    selected(c(setdiff(selected(), ids), setdiff(ids, selected())))
  })
  list(rows = rows, table = shiny::reactive(rows()[selected(), , drop = FALSE]))
}

# Submission handler for the prediction form.
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
                       "â€” processing requires the Linux prediction worker.")
      output[[paste0("results", suffix)]] <- renderText(message)
    }, error = function(e) {
      output[[paste0("results", suffix)]] <- renderText(conditionMessage(e))
      showNotification(conditionMessage(e), type = "error")
    })
  }
  list(submit_prediction = submit_prediction)
}
