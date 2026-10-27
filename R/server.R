server <- function(input, output, session) {
    updateSelectizeInput(session, "hla2", choices = HLAset$HLA, server = TRUE)
    addTooltip(session, id = "input1", title = "Enter an epitope sequence (8-12aa)\nIf lower and uppercase is used, uppercase aa are fixed positions",
        placement = "right", trigger = "hover")
    addTooltip(session, id = "input4", title = "(opt) save the results under a project folder, defaults to 'default'",
        placement = "right", trigger = "hover")
    addTooltip(session, id = "input6", title = "Use a custom sequence for the epitope identification instead of the seuqence database.\n\n                           Useful for mutated sequences",
        placement = "right", trigger = "hover")
    observeEvent(input$run_button, submit_prediction(""))
    observe({
        if (input$toggle == TRUE) {
            shinyjs::show(id = "Color-1")
            shinyjs::show(id = "Color-2")
            shinyjs::show(id = "Color-3")
            shinyjs::show(id = "Color-4")
        }
        else {
            shinyjs::hide(id = "Color-1")
            shinyjs::hide(id = "Color-2")
            shinyjs::hide(id = "Color-3")
            shinyjs::hide(id = "Color-4")
        }
    })
    directory <- config$results_dir
    directory_running <- config$queue_dir
    updateSelectInput(session, "file_selection", choices = basename(list.dirs(directory, recursive = FALSE)))
    updateSelectInput(session, "file_selection_server", choices = basename(list.dirs(directory, recursive = FALSE)))
    observeEvent(input$file_selection, {
        directory <- file.path(directory, input$file_selection)
        updateSelectInput(session, "subdir_selection", choices = list.files(directory, full.names = FALSE))
    })
    observeEvent(input$file_selection_server, {
        directory <- file.path(directory, input$file_selection_server)
        updateSelectInput(session, "subdir_selection_server", choices = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), selected = grep(pattern = "PAX", x = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), value = T)[1])
    })
    output$running_directory_list <- renderText({
        files <- files_reactive_running()
        if (files == "") {
            paste("No runs in queue, waiting for input")
        }
        else {
            paste(files, collapse = ", ")
        }
    })
    files_reactive <- reactivePoll(intervalMillis = 1000, session = session, checkFunc = function() {
        list.files(directory)
    }, valueFunc = function() {
        files <- list.files(directory)
        updateSelectInput(session, "file_selection", choices = files)
        files
    })
    files_reactive_server <- reactivePoll(intervalMillis = 1000, session = session, checkFunc = function() {
        list.files(directory)
    }, valueFunc = function() {
        files <- list.files(directory)
        updateSelectInput(session, "file_selection_server", choices = files)
        files
    })
    files_reactive_running <- reactivePoll(1000, session,
        checkFunc = function() queue_listing(config$queue_dir),
        valueFunc = function() queue_listing(config$queue_dir))
    files_reactive_finished <- reactivePoll(intervalMillis = 1000, session = session, checkFunc = function() {
        files <- list.files("recent_finish")
        if (length(files) > 0) {
            result <- as.data.frame(fread(paste0("recent_finish", "/", files[1]), sep = " "))
            paste0("The run with epitope ", result[1], " on ", result[2], " finished!<br/>Find the results under the directory ",
                result[4])
        }
        else {
            paste0("")
        }
    }, valueFunc = function() {
        files <- list.files("recent_finish")
        if (length(list.files("recent_finish")) > 0) {
            result <- as.data.frame(fread(paste0("recent_finish", "/", files[1]), sep = " "))
            paste0("The run with epitope ", result[1], " on ", result[2], " finished!<br/>Find the results under the directory ",
                result[4])
        }
        else {
            paste0("")
        }
    })
    output$finished_directory_list <- renderText({
        files <- files_reactive_finished()
        if (files == "") {
            paste("Waiting for jobs to finish")
        }
        else {
            paste0(files)
        }
    })
    observeEvent(session, {
        updateSelectInput(session, "file_selection", choices = files_reactive())
    })
    observeEvent(session, {
        updateSelectInput(session, "file_selection_server", choices = files_reactive_server())
    })
    new_tables <- TRUE
    observeEvent(input$view_button_done, {
        print(paste0(directory, "/", input$file_selection, "/", input$subdir_selection))
        if (file.exists(paste0(directory, "/", input$file_selection, "/", input$subdir_selection))) {
            df_done <- read.csv(paste0(directory, "/", input$file_selection, "/", input$subdir_selection),
                sep = "\t")
            if (!grepl(x = paste0(directory, "/", input$file_selection, "/", input$subdir_selection),
                pattern = "_ICR_")) {
                print("hat nicht ICR")
                df_done <- unique(df_done)
            }
            if (new_tables == FALSE) {
                print("this is the old plot")
                output$table_server <- DT::renderDataTable({
                  targets <- c()
                  if (grepl(x = input$subdir_selection, pattern = "search|_EX_")) {
                    targets = c(0, 13, 14, 15, 17, 18, 28)
                    if (grepl(x = input$subdir_selection, pattern = "scores")) {
                      targets = c(0, 13, 14, 15, 17, 18)
                    }
                  }
                  if (grepl(x = input$subdir_selection, pattern = "lowCompl")) {
                    targets = c(0)
                  }
                  if (grepl(x = input$subdir_selection, pattern = "epitope_search.tsv")) {
                    targets = c(0)
                  }
                  DT::datatable(df_done, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                    order = list(2, "asc"), dom = "Bfrtip", colReorder = list(realtime = FALSE), buttons = list(list(extend = "colvis",
                      columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                      show = ":hidden")), columnDefs = list(list(visible = FALSE, targets = targets))))
                }, server = TRUE, processing = TRUE, DT.warn.size = FALSE, reactivity = F)
            }
            else {
                print("this is the new plot")
                if (!grepl(x = paste0(directory, "/", input$file_selection, "/", input$subdir_selection),
                  pattern = "_ICR_")) {
                  output$table_server <- DT::renderDataTable({
                    targets <- c()
                    if (grepl(x = input$subdir_selection, pattern = "search|EX")) {
                      targets = c(0, 13, 14, 15, 17, 18, 28)
                      if (grepl(x = input$subdir_selection, pattern = "scores")) {
                        targets = c(0, 13, 14, 15, 17, 18)
                      }
                    }
                    if (grepl(x = input$subdir_selection, pattern = "lowCompl")) {
                      targets = c(0)
                    }
                    if (grepl(x = input$subdir_selection, pattern = "epitope_search.tsv")) {
                      targets = c(0, 13)
                    }
                    if (!grepl(x = paste0(directory, "/", input$file_selection, "/", input$subdir_selection),
                      pattern = "_ICR_")) {
                      print("hat nicht ICR loop 2")
                      df_done$Binding <- factor(df_done$Binding, levels = c("strong", "weak", "non"))
                      df_done$Uniprot_db <- factor(df_done$Uniprot_db, levels = c("reviewed", "unreviewed",
                        "not included"))
                      df_done$Expression <- factor(df_done$Expression, levels = c("Expression", "no Expression"))
                    }
                    DT::datatable(df_done, filter = "top", extensions = "Buttons", options = list(paging = TRUE,
                      pageLength = 25, scrollX = TRUE, searching = TRUE, ordering = TRUE, dom = "Bfrtip",
                      buttons = list(list(extend = "colvis", columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup",
                        text = "Show All", show = ":hidden")), columnDefs = list(list(visible = FALSE,
                        targets = targets))))
                  })
                }
                else {
                  cnames <- colnames(df_done)
                  for (i in 1:dim(df_done)[1]) {
                    df_done[i, 1] <- gsub("\\.$", "", gsub("\\]", "", gsub("\\[", "", df_done[i, 1])))
                  }
                  df_done <- cbind(mismatch = 0:(dim(df_done)[1] - 1), df_done)
                  for (i in 2:dim(df_done)) {
                    df_done[i, ] <- substr(sprintf("%.10f", df_done[i, ]), 0, 10)
                    df_done[i, ] <- as.numeric(df_done[i, ])
                  }
                  df_done <- unique(df_done)
                  output$table_server <- DT::renderDataTable({
                    DT::datatable(df_done, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                      order = list(1, "desc"), dom = "Bfrtip", colReorder = list(realtime = FALSE), buttons = list(list(extend = "colvis",
                        columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                        show = ":hidden")), columnDefs = list(list(visible = FALSE, targets = c(0), colnames = cnames))))
                  })
                }
            }
        }
    })
    observeEvent(input$view_score_button, {
        score_file <- grep(x = list.files(paste0(directory, "/", input$file_selection), full.names = T),
            pattern = "epitope_search_scores_mlp.tsv", value = TRUE)
        if (length(score_file) > 0) {
            df_done <- read.csv(paste0(score_file[1]), sep = "\t")
            cnames <- colnames(df_done)
            df_done <- unique(df_done)
            print(dim(df_done))
            output$table_server <- DT::renderDataTable({
                if (new_tables == FALSE) {
                  DT::datatable(df_done, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                    order = list(2, "asc"), dom = "Bfrtip", colReorder = list(realtime = FALSE), buttons = list(list(extend = "colvis",
                      columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                      show = ":hidden")), columnDefs = list(list(visible = FALSE, targets = c(0, 13,
                      14, 15, 17, 18)))))
                }
                else {
                  df_done$Binding <- factor(df_done$Binding, levels = c("strong", "weak", "non"))
                  df_done$Uniprot_db <- factor(df_done$Uniprot_db, levels = c("reviewed", "unreviewed",
                    "not included"))
                  df_done$Expression <- factor(df_done$Expression, levels = c("Expression", "no Expression"))
                  DT::datatable(df_done, filter = "top", extensions = "Buttons", options = list(paging = TRUE,
                    pageLength = 25, scrollX = TRUE, searching = TRUE, ordering = TRUE, dom = "Bfrtip",
                    buttons = list(list(extend = "colvis", columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup",
                      text = "Show All", show = ":hidden")), columnDefs = list(list(visible = FALSE,
                      targets = c(0, 13, 14, 15, 17, 18)))))
                }
            }, server = FALSE, processing = TRUE, DT.warn.size = FALSE, reactivity = F)
        }
    })
    observeEvent(input$view_icr_button, {
        score_file <- grep(x = list.files(paste0(directory, "/", input$file_selection), full.names = T),
            pattern = "ICR_PAX_v5", value = TRUE)
        if (length(score_file) > 0) {
            df_done <- read.csv(paste0(score_file[1]), sep = "\t")
            cnames <- colnames(df_done)
            for (i in 1:dim(df_done)[1]) {
                df_done[i, 1] <- gsub("\\.$", "", gsub("\\]", "", gsub("\\[", "", df_done[i, 1])))
            }
            df_done[, 1] <- as.numeric(df_done[, 1])
            df_done <- cbind(mismatch = 0:(dim(df_done)[1] - 1), df_done)
            for (i in 2:dim(df_done)[1]) {
                df_done[i, ] <- substr(sprintf("%.10f", df_done[i, ]), 0, 10)
            }
            df_done[, 1] <- as.numeric(df_done[, 1])
            output$table_server <- DT::renderDataTable({
                DT::datatable(df_done, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                  order = list(1, "desc"), dom = "Bfrtip", colReorder = list(realtime = FALSE), buttons = list(list(extend = "colvis",
                    columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                    show = ":hidden")), columnDefs = list(list(visible = FALSE, targets = c(0), colnames = cnames))))
            })
        }
    })
    observeEvent(input$view_icr_exp_button, {
        score_file <- grep(x = list.files(paste0(directory, "/", input$file_selection), full.names = T),
            pattern = "ICR_PAX_v4", value = TRUE)
        if (length(score_file) > 0) {
            df_done <- read.csv(paste0(score_file[1]), sep = "\t")
            cnames <- colnames(df_done)
            for (i in 1:dim(df_done)[1]) {
                df_done[i, 1] <- gsub("\\.$", "", gsub("\\]", "", gsub("\\[", "", df_done[i, 1])))
            }
            df_done <- cbind(mismatch = 0:(dim(df_done)[1] - 1), df_done)
            score_file_expression <- gsub("_ICR_", "_EX_", score_file)
            df_done_expression <- read.csv(paste0(score_file_expression[1]), sep = "\t")
            df_expression <- df_done_expression[c(2, which(colnames(df_done_expression) %in% cnames))]
            colnames(df_expression)[1] <- "mismatch"
            df_expression_table <- df_expression %>% group_by(mismatch) %>% summarise(across(everything(),
                sum))
            df_expression_table <- as.data.frame(df_expression_table)
            if (dim(df_done)[1] > dim(df_expression_table)[1]) {
                for (m in df_done$mismatch) {
                  if (!(m %in% df_expression_table$mismatch)) {
                    df_expression_table <- rbind(df_expression_table, c(m, rep(0, dim(df_expression_table)[2] -
                      1)))
                  }
                }
            }
            df_expression_table <- left_join(df_done[c(1, 2)], df_expression_table, by = "mismatch")
            df_to_show <- df_expression_table
            output$table_server <- DT::renderDataTable({
                DT::datatable(df_to_show, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                  order = list(1, "desc"), dom = "Bfrtip", colReorder = list(realtime = FALSE), buttons = list(list(extend = "colvis",
                    columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                    show = ":hidden")), columnDefs = list(list(visible = FALSE, targets = c(0), colnames = cnames))))
            })
        }
    })
    observeEvent(input$view_icr_button_testing, {
        score_file <- grep(x = list.files(paste0(directory, "/", input$file_selection), full.names = T),
            pattern = "ICR_PAX_v4", value = TRUE)
        if (length(score_file) > 0) {
            df_done <- read.csv(paste0(score_file[1]), sep = "\t")
            cnames <- colnames(df_done)
            for (i in 1:dim(df_done)[1]) {
                df_done[i, 1] <- gsub("\\.$", "", gsub("\\]", "", gsub("\\[", "", df_done[i, 1])))
            }
            df_done <- cbind(mismatch = 0:(dim(df_done)[1] - 1), df_done)
            score_file_expression <- gsub("_ICR_", "_EX_", score_file)
            df_done_expression <- read.csv(paste0(score_file_expression[1]), sep = "\t")
            df_expression <- df_done_expression[c(2, which(colnames(df_done_expression) %in% cnames))]
            colnames(df_expression)[1] <- "mismatch"
            df_expression_table <- df_expression %>% group_by(mismatch) %>% summarise(across(everything(),
                sum))
            df_expression_table <- as.data.frame(df_expression_table)
            if (dim(df_done)[1] > dim(df_expression_table)[1]) {
                for (m in df_done$mismatch) {
                  if (!(m %in% df_expression_table$mismatch)) {
                    df_expression_table <- rbind(df_expression_table, c(m, rep(0, dim(df_expression_table)[2] -
                      1)))
                  }
                }
            }
            df_expression_table <- left_join(df_done[c(1, 2)], df_expression_table, by = "mismatch")
            df_to_show <- df_done
            output$table_server <- DT::renderDataTable({
                DT::datatable(df_to_show, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                  order = list(1, "desc"), dom = "Bfrtip", colReorder = list(realtime = FALSE), buttons = list(list(extend = "colvis",
                    columns = c(0:dim(df_done)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                    show = ":hidden"), list(extend = "collection", text = "My Action", action = DT::JS("function() {\n                                                                                                                    var node = this[0].node;\n                                                                                                                    var value = $(node).attr('data-value') || 0;\n                                                                                                                    var value = df_expression_table;\n                                                                                                                    $(node).attr('data-value', value);\n                                                                                                                      Shiny.setInputValue('dt_to_show', df_expression_table, {priority: 'event'});\n                                                                                                                      }"))),
                  columnDefs = list(list(visible = FALSE, targets = c(0), colnames = cnames))))
            })
        }
    })
    observeEvent(input$view_button_done_tum, {
        print(paste0(directory, "/", input$file_selection_tum, "/", input$subdir_selection_tum))
        df_done_tum <- read.csv(paste0(directory, "/", input$file_selection_tum, "/", input$subdir_selection_tum),
            sep = "\t")
        df_done_tum <- unique(df_done_tum)
        output$table_server_tum <- DT::renderDataTable({
            DT::datatable(df_done_tum, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
                order = list(0, "desc"), dom = "Bfrtip", colReorder = TRUE, buttons = list(list(extend = "colvis",
                  columns = c(0:dim(df_done_tum)[2] - 1)), list(extend = "colvisGroup", text = "Show All",
                  show = ":hidden")), columnDefs = list(list(visible = FALSE, targets = c(0, 13, 14,
                  15, 17, 18, 28)))))
        })
    })
    files_reactive_server <- reactivePoll(intervalMillis = 1000, session = session, checkFunc = function() {
        list.files(directory)
    }, valueFunc = function() {
        files <- list.files(directory)
        updateSelectInput(session, "file_selection_server", choices = files)
        files
    })
    updateSelectInput(session, "file_selection_server", choices = list.dirs(directory, recursive = FALSE))
    observeEvent(input$file_selection_server, {
        directory <- file.path(directory, input$file_selection_server)
        updateSelectInput(session, "subdir_selection_server", choices = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), selected = grep(pattern = "PAX", x = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), value = T)[1])
        updateSelectInput(session, "subdir_selection_server_de1", choices = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), selected = grep(pattern = "gtex_v8", x = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), value = T)[1])
        updateSelectInput(session, "subdir_selection_server_de2", choices = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), selected = grep(pattern = "TCGA_q05", x = list.files(directory, full.names = FALSE,
            pattern = "_EX_"), value = T)[1])
        updateSelectInput(session, "icr_select_server", choices = list.files(directory, full.names = FALSE,
            pattern = "_ICR_"), selected = grep(pattern = "PAX", x = list.files(directory, full.names = FALSE,
            pattern = "_ICR_"), value = T)[1])
        updateSelectInput(session, "mlp_select_server", choices = list.files(directory, full.names = FALSE,
            pattern = "_EX_|epitope_search_scores_mlp"), selected = list.files(directory, full.names = FALSE,
            pattern = "epitope_search_scores_mlp")[1])
    })
    observeEvent(input$subdir_selection_server_de1, {
        req(input$subdir_selection_server_de1)
        if (file.exists(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server_de1))) {
            df <- read.csv(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server_de1),
                sep = "\t")
            startselection <- colnames(df)[29]
            if ("breast" %in% colnames(df)[-c(1:28)]) {
                startselection <- "breast"
            }
            updateSelectInput(session, "tissue_select_df1_server", choices = colnames(df)[-c(1:28)],
                selected = startselection)
        }
    })
    observeEvent(input$subdir_selection_server_de2, {
        req(input$subdir_selection_server_de2)
        if (file.exists(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server_de2))) {
            df <- read.csv(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server_de2),
                sep = "\t")
            startselection <- colnames(df)[30]
            if ("tcga.brca" %in% colnames(df)[-c(1:28)]) {
                startselection <- "tcga.brca"
            }
            updateSelectInput(session, "tissue_select_df2_server", choices = colnames(df)[-c(1:28)],
                selected = startselection)
        }
    })
    observeEvent(input$subdir_selection_server, {
        req(input$file_selection_server, input$subdir_selection_server)
        if (file.exists(file.path(directory, input$file_selection_server, input$subdir_selection_server))) {
            df <- data_expression()
            updateSelectInput(session, "tissue_select_server", choices = colnames(df)[-c(1:28)])
            updateSelectInput(session, "tissue_select_plot_server", choices = colnames(df)[-c(1:28)])
        }
    })
    data_expression <- reactive({
        req(input$file_selection_server, input$subdir_selection_server)
        df <- as.data.frame(read.csv(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server),
            sep = "\t"))
    })
    reference_tissue <- reactive({
        req(input$subdir_selection_server, input$tissue_select_server)
        name <- sub("\\.tsv$", "", strsplit(input$subdir_selection_server, "_EX_", fixed = TRUE)[[1]][2])
        path <- file.path(config$data_dir, "Datenbanken_Ordner_shiny", paste0(name, ".tsv"))
        validate(need(file.exists(path), "The matching reference expression database is missing"))
        headers <- names(data.table::fread(path, nrows = 0))
        column <- match(input$tissue_select_server, make.names(headers, unique = TRUE))
        validate(need(!is.na(column) && column > 1, "Tissue not found in the reference database"))
        data.table::fread(path, select = c(1L, column), data.table = FALSE)
    })
    output$distPlot_server_new <- renderPlotly({
        req(input$tissue_select_server)
        distribution_plot(expression_distribution(reference_tissue(), data_expression(), input$tissue_select_server))
    })
    bee_selection <- link_plot_selection("bee_source", session)
    volcano_selection <- link_plot_selection("vulcano_source", session)
    output$beeswarm_plot_gg_server <- renderPlotly({
        validate(need(length(input$tissue_select_plot_server) > 0, "Please select at least one tissue."))
        req(input$subdir_selection_server)
        if (length(input$tissue_select_plot_server) != 0) {
            df2 <- data_expression()
            validate(need(dim(df2)[1] > 0, "Expression file seems to be empty! Try another file. "))
            data_gather2 <- beeswarm_data(df2, input$tissue_select_plot_server, my_range_server())
            bee_selection$rows(setNames(data_gather2[, c("epitope", "distance", "bind", "gene_id", "tissue", "expression")],
                c("Epitope", "Mismatch", "Bind", "Gene", "Tissue", "Expression")))
            key <- as.character(seq_len(nrow(data_gather2)))
            p <- ggplot(data_gather2, aes(x = tissue, y = expression, color = tissue, shape = bind, infotext = info,
                key = key)) + geom_beeswarm(cex = 1.2) + theme(axis.text.x = element_text(angle = 60,
                hjust = 1)) + ggtitle(paste0("Expression in ", strsplit(x = input$subdir_selection_server,
                split = ".", fixed = T)[[1]][1], " per tissue and values between ", my_range_server()[1],
                " and ", my_range_server()[2]))
            p_bee <- ggplotly(p, tooltip = "infotext", source = "bee_source")
            ggp700(p_bee, "info") %>% layout(dragmode = "lasso") %>% event_register("plotly_click") %>% event_register("plotly_selected") %>% event_register("plotly_deselect")
        }
    })
    output$selected_data_table <- DT::renderDataTable({
        DT::datatable(bee_selection$table(), extensions = "Buttons", options = list(paging = TRUE, pageLength = 25,
            scrollX = TRUE, searching = TRUE, ordering = TRUE, dom = "Bfrtip", buttons = c("copy", "csv",
                "excel", "pdf")))
    })
    output$vulcano_plot_gg_server <- renderPlotly({
        req(input$tissue_select_df1_server, input$tissue_select_df2_server)
        cols_interest <- c(1, 2, 4, 28)
        df_de1 <- read.csv(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server_de1),
            sep = "\t")
        de1 <- as.data.frame(cbind(df_de1[, cols_interest], df_de1[, input$tissue_select_df1_server]))
        colnames(de1) <- c("epitope", "distance", "bind", "gene_id", input$tissue_select_df1_server)
        de1 <- unique(de1)
        df_de2 <- read.csv(paste0(directory, "/", input$file_selection_server, "/", input$subdir_selection_server_de2),
            sep = "\t")
        de2 <- as.data.frame(cbind(df_de2[, cols_interest], df_de2[, input$tissue_select_df2_server]))
        colnames(de2) <- c("epitope", "distance", "bind", "gene_id", input$tissue_select_df2_server)
        de2 <- unique(de2)
        de_data <- full_join(de1, de2)
        de_data[is.na(de_data)] <- 0
        de_matrix <- as.matrix.data.frame(de_data[, 5:dim(de_data)[2]])
        de_matrix <- matrix(c(trunc(as.numeric(de_matrix)), round(as.numeric(de_matrix)), ceiling(as.numeric(de_matrix))),
            ncol = ncol(de_matrix) * 3)
        colnames(de_matrix) <- as.vector(outer(colnames(de_data[5:dim(de_data)[2]]), c("_low", "_mid",
            "_high"), paste0, sep = ""))
        rownames(de_matrix) <- de_data$gene_id
        ttestRat <- function(df, grp1, grp2) {
            x = df[grp1]
            y = df[grp2]
            x = as.numeric(x)
            y = as.numeric(y)
            results = t.test(x, y)
            results$p.value
        }
        pvalue = apply(de_matrix, 1, ttestRat, grp1 = c(1, 3, 5), grp2 = c(2, 4, 6))
        rat = log2(de_matrix)
        control = apply(rat[, c(1, 3, 5)], 1, mean)
        test = apply(rat[, c(2, 4, 6)], 1, mean)
        foldchange <- control - test
        results = cbind(de_data[, 5:dim(de_data)[2]], foldchange, pvalue)
        results = as.data.frame(results)
        results$gene_id <- de_data$gene_id
        results$epitope <- de_data$epitope
        results$distance <- de_data$distance
        results$bind <- de_data$bind
        results$info <- paste0(results$epitope, "\ndistance", results$distance, "\nbind: ", results$bind,
            "\n", results$gene_id, "\n", input$tissue_select_df1_server, ": ", results[, input$tissue_select_df1_server],
            "\n", input$tissue_select_df2_server, ": ", results[, input$tissue_select_df2_server], "\n",
            results$bind_prediction)
        volcano_selection$rows(data.frame(Epitope = results$epitope, Mismatch = results$distance,
            Bind = results$bind, Pvalue = results$pvalue, Gene = results$gene_id,
            Tissue1 = paste(input$tissue_select_df1_server, results[[input$tissue_select_df1_server]], sep = ": "),
            Tissue2 = paste(input$tissue_select_df2_server, results[[input$tissue_select_df2_server]], sep = ": ")))
        results$key <- as.character(seq_len(nrow(results)))
        p <- ggplot(data = results, aes(x = foldchange, y = -1 * log10(pvalue), color = bind, infotext = info,
            key = key)) + geom_point(size = 2) + ggtitle(paste0("Differential Expression of ", input$tissue_select_df1_server,
            " vs ", input$tissue_select_df2_server)) + ylab("-log10(p-value)") + xlab("log2 fold change")
        p_vul <- ggplotly(p, tooltip = "infotext", source = "vulcano_source")
        ggp700(p_vul, "info") %>% layout(dragmode = "lasso") %>% event_register("plotly_click") %>% event_register("plotly_selected") %>% event_register("plotly_deselect") %>% onRender(js)
    })
    output$selected_data_table_vulcano <- DT::renderDataTable({
        DT::datatable(volcano_selection$table(), extensions = "Buttons", options = list(paging = TRUE,
            pageLength = 25, scrollX = TRUE, searching = TRUE, ordering = TRUE, dom = "Bfrtip", buttons = c("copy",
                "csv", "excel", "pdf")))
    })
    selected <- reactive({
        event_data("plotly_selected", source = "master")
    })
    output$blosum_server_plot2 <- renderPlotly({
        req(input$file_selection_server, input$mlp_select_server)
        req(results_blos())
        blosum_epitope <- input$blosum_epitope
        if (blosum_epitope == "") {
            blosum_epitope <- strsplit(input$mlp_select_server, "_")[[1]][1]
        }
        results <- results_blos()$data
        size_epitopes <- as.character(dim(results)[1])
        x_axis <- input$x_axis
        if (input$blosum_scale == "no scale") {
            p <- ggplot(results, aes(x = .data[[x_axis]], y = BLOSUM, color = Binding, infotext = info)) +
                geom_point() + ggtitle(paste0("BLOSUM", input$blosum_radio, " score of ", size_epitopes,
                " epitopes to the search epitope ", blosum_epitope, " versus ", x_axis)) + ylab(paste0("BLOSUM",
                input$blosum_radio)) + xlab(x_axis) + theme_bw()
            print(selected()$y)
            obj <- tryCatch(results[which(results$BLOSUM %in% selected()$y & round(results[[input$x_axis]],
                8) %in% round(as.numeric(selected()$x), 8), ), ], error = function(e) {
                NULL
            })
            if (nrow(obj) != 0) {
                p <- p + geom_point(data = obj, color = "yellow", size = 4)
            }
            ggp700_source(p, "info", "master") %>% onRender(js)
        }
        else {
            p <- ggplot(results, aes(x = log(.data[[x_axis]]), y = BLOSUM, color = Binding, infotext = info)) +
                geom_point() + ggtitle(paste0("BLOSUM", input$blosum_radio, " score of ", size_epitopes,
                " epitopes to the search epitope ", blosum_epitope, " versus ", x_axis)) + ylab(paste0("BLOSUM",
                input$blosum_radio)) + xlab(paste0("log2( ", x_axis, " )")) + theme_bw()
            obj <- tryCatch(results[which(results$BLOSUM %in% selected()$y & round(results[[input$x_axis]],
                8) %in% round(as.numeric(selected()$x), 8), ), ], error = function(e) {
                NULL
            })
            print(dim(obj))
            if (nrow(obj) != 0) {
                p <- p + geom_point(data = obj, color = "yellow", size = 4)
            }
            ggp700_source(p, "info", "master") %>% onRender(js)
        }
    })
    results_blos <- reactive({
        req(input$file_selection_server, input$mlp_select_server)
        df2 <- as.data.frame(read.csv(paste0(directory, "/", input$file_selection_server, "/", input$mlp_select_server),
            sep = "\t"))
        cols_interest <- c(1, 2, 4, 12, 20, 21)
        print("vor dem if")
        print(grepl(x = input$mlp_select_server, pattern = "_EX_"))
        if (grepl(x = input$mlp_select_server, pattern = "_EX_")) {
            cols_interest <- c(1, 2, 4, 12, 20, 21, 29:dim(df2)[2])
        }
        print("Nach dem IF")
        blosum_df <- unique(df2[, cols_interest])
        colnames(blosum_df)[1:6] <- c("Epitope", "Mismatch", "Binding", "Affinity", "Processing_score_flanks",
            "Presentation_score_flanks")
        updateSelectInput(session, "x_axis", choices = colnames(blosum_df)[4:length(blosum_df)], selected = "Presentation_score_flanks")
        BLOSUM <- as.matrix(fread(paste0("blosum/BLOSUM", input$blosum_radio, ".tsv"), header = T), rownames = 1)
        blosum_epitope <- input$blosum_epitope
        if (blosum_epitope == "") {
            blosum_epitope <- strsplit(input$mlp_select_server, "_")[[1]][1]
        }
        blosum_df$BLOSUM <- unlist(lapply(blosum_df$Epitope, function(x) pwalign::pairwiseAlignment(blosum_epitope,
            x, substitutionMatrix = BLOSUM, gapOpening = -2, gapExtension = -8, scoreOnly = TRUE)))
        results <- blosum_df
        results$info <- paste0(results$Epitope, "\nmismatch ", results$Mismatch, "\nbind: ", results$Binding,
            "\nBLOSUM: ", results$BLOSUM, "\nPresentation_score_flanks: ", results$Presentation_score_flanks)
        results$Binding <- factor(results$Binding, levels = c("strong", "weak", "non"))
        results$tcr_ranking <- 1.1 * ifelse(results$Binding == "strong", yes = 2, no = ifelse((results$Binding ==
            "weak"), yes = 1.5, no = 1)) + 16 * (results$BLOSUM/max(results$BLOSUM)) + 5 * results$Presentation_score_flanks
        results$tcr_ranking <- 1.1 * ifelse(results$Binding == "strong", yes = 2, no = ifelse((results$Binding ==
            "weak"), yes = 1.5, no = 1)) + 3 * (results$BLOSUM/max(results$BLOSUM)) + 2 * results$Presentation_score_flanks
        list(data = results)
    })
    output$download_blos <- downloadHandler(filename = function() {
        file2 <- paste0(gsub(x = input$mlp_select_server, pattern = "epitope_search_scores_mlp.tsv",
            replacement = "blosum"), ".xlsx")
    }, content = function(file) {
        file2 <- as.data.frame(tryCatch(results_blos()$data[which(results_blos()$data$BLOSUM %in% selected()$y &
            round(results_blos()$data[[input$x_axis]], 8) %in% round(as.numeric(selected()$x), 8), ),
            ], error = function(e) {
            NULL
        }))
        print(dim(file2))
        writexl::write_xlsx(x = file2, file, col_names = T, format_headers = T)
    }, contentType = "application/zip")
    output$se <- DT::renderDataTable(tryCatch(results_blos()$data[which(results_blos()$data$BLOSUM %in%
        selected()$y & round(results_blos()$data[[input$x_axis]], 8) %in% round(as.numeric(selected()$x),
        8), ), ], error = function(e) {
        NULL
    }), filter = "top", extensions = "Buttons", options = list(autoWidth = TRUE, paging = TRUE, pageLength = 15,
        dom = "Bfrtip", buttons = c("copy", "csv", "excel", "pdf")))
    output$downloadData <- downloadHandler(filename = function() {
        file2 <- paste0(input$subdir_selection, ".xlsx")
    }, content = function(file) {
        file2 <- as.data.frame(fread(paste0(directory, "/", input$file_selection, "/", input$subdir_selection),
            sep = "\t"))
        print(dim(file2))
        file2 <- unique(file2)
        writexl::write_xlsx(x = file2, file, col_names = T, format_headers = T)
    }, contentType = "application/zip")
    output$download_zip <- downloadHandler(filename = function() {
        zip_file <- paste0(input$file_selection, ".zip")
    }, content = function(file) {
        directory <- file.path(directory, input$file_selection)
        zip_file <- paste0(input$file_selection, ".zip")
        zip(file, files = list.files(directory, full.names = TRUE, recursive = TRUE), mode = "cherry-pick")
        print("zip ersellt")
    }, contentType = "application/zip")
    observeEvent(input$download_button, {
        directory <- file.path(directory, input$file_selection)
        zip_file <- paste0(input$file_selection, ".zip")
        zip(zip_file, files = list.files(directory, full.names = TRUE, recursive = TRUE))
        output$status_message <- renderText({
            sprintf("Zip file created: %s", zip_file)
        })
        print("fertig gezippt")
        print("hier da hast du")
        output$download_zip <- downloadHandler(filename = function() {
            zip_file
        }, content = function(file) {
            file.copy(zip_file, file)
        })
        shinyjs::show("download_zip")
    })
    observeEvent(input$view_button, {
        df <- read.csv(input$file$datapath, sep = "\t")
        output$table <- DT::renderDataTable({
            DT::datatable(df, extensions = "Buttons", options = list(pageLength = 20, order = list(0,
                "desc"), dom = "Bfrtip", buttons = list(list(extend = "colvis", columns = c(0:dim(df)[2] -
                1)), list(extend = "colvisGroup", text = "Show All", show = ":hidden")), columnDefs = list(list(visible = FALSE,
                targets = c(0, 13, 14, 15, 17, 18, 28)))))
        })
    })
    observeEvent(input$file_ex_de1, {
        req(input$file_ex_de1)
        df_de1 <- read.table(input$file_ex_de1$datapath, header = TRUE, sep = "\t")
        updateSelectInput(session, "tissue_select_df1", choices = colnames(df_de1)[c(29:dim(df_de1)[2])])
    })
    observeEvent(input$file_ex_de2, {
        req(input$file_ex_de2)
        df_de2 <- read.table(input$file_ex_de2$datapath, header = TRUE, sep = "\t")
        updateSelectInput(session, "tissue_select_df2", choices = colnames(df_de2)[c(29:dim(df_de2)[2])])
    })
    observeEvent(input$plot_button, {
        req(input$file_ex)
        df <- read.table(input$file_ex$datapath, header = TRUE, sep = "\t")
        updateSelectInput(session, "tissue_select", choices = colnames(df)[-c(1:28)])
        updateSelectInput(session, "tissue_select_plot", choices = colnames(df)[-c(1:28)])
    })
    output$dotplot_tissue <- renderPlotly({
        if (input$tissue_select != "") {
            df2 <- as.data.frame(read.table(input$file_ex$datapath, header = TRUE, sep = "\t"))
            cols_interest <- c(1, 2, 4, 28)
            tissue_df <- as.data.frame(cbind(df2[, cols_interest], df2[, input$tissue_select]))
            tissue_df[, 5] <- as.numeric(tissue_df[, 5])
            tissue_df <- tissue_df[order(tissue_df[, 5]), ]
            colnames(tissue_df) <- c("epitope", "distance", "bind", "gene_id", "expression")
            tissue_df$ids <- 1:dim(tissue_df)[1]
            p2 <- ggplot(tissue_df, aes(x = ids, y = expression, Protein = gene_id, color = bind)) +
                geom_point() + scale_y_continuous(breaks = seq(min(tissue_df$expression), max(tissue_df$expression),
                (max(tissue_df$expression) - min(tissue_df$expression))/5))
            ggplotly(p2)
        }
    })
    output$distPlot <- renderPlotly({
        if (input$tissue_select != "") {
            df <- as.data.frame(read.table(input$file_ex$datapath, header = TRUE, sep = "\t"))
            if (input$tissue_select != "") {
                tissue_df <- df[, input$tissue_select]
                x <- tissue_df
                bins <- seq(min(x), max(x), length.out = input$bins + 1)
                df <- df[order(df[, input$tissue_select], decreasing = T), ]
                for (i in 1:(length(bins) - 1)) {
                  df$info[which(df[, input$tissue_select] >= bins[i] & df[, input$tissue_select] <= bins[i +
                    1])] <- str_c((paste0(df[which(df[, input$tissue_select] >= bins[i] & df[, input$tissue_select] <=
                    bins[i + 1]), ]$ENST, ": ", df[which(df[, input$tissue_select] >= bins[i] & df[,
                    input$tissue_select] <= bins[i + 1]), input$tissue_select])), collapse = "\n")
                }
                p <- ggplot(data = df, aes(.data[[input$tissue_select]], infotext = .data[["info"]])) +
                  geom_histogram(bins = input$bins, breaks = bins, fill = "steelblue3", colour = "grey30") +
                  xlab("Gene Expression Distribution ") + ylab("Number of proteins with that expression") +
                  theme_minimal()
                ggplotly(p, tooltip = "infotext")
            }
        }
    })
    data_gather2 <- reactive({
        req(input$file_ex)
        data2 <- as.data.frame(read.csv(input$file$datapath, sep = "\t"))
        data2 %>% gather(tissue, expression, 2:dim(data2)[2])
    })
    my_range <- reactive({
        cbind(input$min_tpm[1], input$min_tpm[2])
    })
    my_range_server <- reactive({
        cbind(input$min_tpm_server[1], input$min_tpm_server[2])
    })
    output$beeswarm_plot_gg <- renderPlotly({
        validate(need(length(input$tissue_select_plot) > 0, "Please select at least one tissue."))
        req(input$file_ex)
        if (length(input$tissue_select_plot) != 0) {
            df2 <- as.data.frame(read.table(input$file_ex$datapath, header = TRUE, sep = "\t"))
            cols_interest <- c(1, 2, 4, 28)
            tissue_df <- as.data.frame(cbind(df2[, cols_interest], df2[, 29:dim(df2)[2]]))
            data_gather <- df2 %>% gather(tissue, expression, 2:dim(df2)[2])
            data_gather <- tissue_df %>% gather(tissue, expression, 5:dim(tissue_df)[2])
            data_gather2 <- data_gather
            data_gather2 <- data_gather[which((my_range()[1] < data_gather$expression) & (data_gather$expression <
                my_range()[2])), ]
            colnames(data_gather2)[1:4] <- c("epitope", "distance", "bind", "gene_id")
            data_gather2$info <- paste0(data_gather2$epitope, "\ndistance: ", data_gather2$distance,
                "\nbind: ", data_gather2$bind, "\n", data_gather2$gene_id, "\n", data_gather2$tissue,
                ": ", data_gather2$expression)
            data_gather2 <- data_gather2[which(data_gather2$tissue %in% input$tissue_select_plot), ]
            p <- ggplot(data_gather2, aes(x = tissue, y = expression, color = tissue, shape = bind, infotext = info)) +
                geom_beeswarm(cex = 1.2) + theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
                ggtitle(paste0("Expression in ", strsplit(x = input$file_ex$name, split = ".", fixed = T)[[1]][1],
                  " per tissue and values between ", my_range()[1], " and ", my_range()[2]))
            ggplotly(p, tooltip = "infotext")
        }
    })
    output$vulcano_plot_gg <- renderPlotly({
        req(input$tissue_select_df1, input$tissue_select_df2)
        cols_interest <- c(1, 2, 4, 28)
        df_de1 <- read.table(input$file_ex_de1$datapath, header = TRUE, sep = "\t")
        de1 <- as.data.frame(cbind(df_de1[, cols_interest], df_de1[, input$tissue_select_df1]))
        colnames(de1) <- c("epitope", "distance", "bind", "gene_id", input$tissue_select_df1)
        de1 <- unique(de1)
        df_de2 <- read.table(input$file_ex_de2$datapath, header = TRUE, sep = "\t")
        de2 <- as.data.frame(cbind(df_de2[, cols_interest], df_de2[, input$tissue_select_df2]))
        colnames(de2) <- c("epitope", "distance", "bind", "gene_id", input$tissue_select_df2)
        de2 <- unique(de2)
        de_data <- full_join(de1, de2)
        de_data[is.na(de_data)] <- 0
        de_matrix <- as.matrix.data.frame(de_data[, 5:dim(de_data)[2]])
        de_matrix <- matrix(c(trunc(as.numeric(de_matrix)), round(as.numeric(de_matrix)), ceiling(as.numeric(de_matrix))),
            ncol = ncol(de_matrix) * 3)
        colnames(de_matrix) <- as.vector(outer(colnames(de_data[5:dim(de_data)[2]]), c("_low", "_mid",
            "_high"), paste0, sep = ""))
        rownames(de_matrix) <- de_data$gene_id
        ttestRat <- function(df, grp1, grp2) {
            x = df[grp1]
            y = df[grp2]
            x = as.numeric(x)
            y = as.numeric(y)
            results = t.test(x, y)
            results$p.value
        }
        pvalue = apply(de_matrix, 1, ttestRat, grp1 = c(1, 3, 5), grp2 = c(2, 4, 6))
        rat = log2(de_matrix)
        control = apply(rat[, c(1, 3, 5)], 1, mean)
        test = apply(rat[, c(2, 4, 6)], 1, mean)
        foldchange <- control - test
        results = cbind(de_data[, 5:dim(de_data)[2]], foldchange, pvalue)
        results = as.data.frame(results)
        results$gene_id <- rownames(results)
        results$bind <- de_data$bind
        results$info <- paste0(results$gene_id, "\nbind: ", results$bind, "\nExpression ", input$tissue_select_df1,
            ": ", results[, input$tissue_select_df1], "\nExpression ", input$tissue_select_df2, ": ",
            results[, input$tissue_select_df2], "\n", results$bind_prediction)
        p <- ggplot(data = results, aes(x = foldchange, y = -1 * log10(pvalue), color = bind, infotext = info)) +
            geom_point(size = 2) + ggtitle(paste0("Differential Expression of ", input$tissue_select_df1,
            " vs ", input$tissue_select_df2)) + ylab("-log10(p-value)") + xlab("log2 fold change")
        ggplotly(p, tooltip = "infotext", height = 800)
    })
    output$icr_server_plot <- renderPlotly({
        req(input$icr_select_server)
        data_icr <- as.data.frame(fread(paste0(directory, "/", input$file_selection_server, "/", input$icr_select_server),
            sep = "\t"))
        for (i in 1:dim(data_icr)[1]) {
            data_icr$ICR[i] <- as.numeric(gsub(pattern = "\\]", replacement = "", x = (gsub(pattern = "\\[",
                replacement = "", x = data_icr$ICR[i]))))
        }
        data_icr$mismatches <- (0:(dim(data_icr)[1] - 1))
        d <- gather(data = data_icr, tissue, CR, 3:dim(data_icr)[2] - 1)
        d$info <- paste0("\nICR: ", d$ICR)
        p <- ggplot(d, aes(x = tissue, y = mismatches, infotext = info)) + geom_tile(aes(fill = CR),
            alpha = 0) + scale_fill_distiller(palette = "YlGnBu", direction = -1) + theme(axis.text.x = element_text(angle = 60,
            hjust = 1)) + labs(title = "On Target off Tumour and Cross Reactivity per tissue", x = "tissues",
            y = "mismatches", )
        ggplotly(p, info = "infotext", height = 800)
    })
    updateSelectizeInput(session, "tcga_ensg", choices = ENSGset$ensg, server = TRUE)
    updateSelectizeInput(session, "tcga_tissues", choices = names(expression_paths), server = TRUE)
    session_store <- reactiveValues()
    expression_data <- reactive({
        validate(need(length(input$tcga_tissues) > 0, "Please select at least one tissue."))
        req(input$tcga_ensg)
        validate(need(all(input$tcga_tissues %in% names(expression_paths)), "Choose an available tissue"))
        do.call(rbind, lapply(input$tcga_tissues, function(tissue)
            read_gene_expression(unname(expression_paths[[tissue]]), input$tcga_ensg)))
    })
    gene_title <- reactive({
        req(input$tcga_ensg)
        name <- ensg_genename$gene[match(input$tcga_ensg, ensg_genename$ENSG)]
        if (is.na(name) || !nzchar(name)) name <- input$tcga_ensg
        paste("Expression of", name, paste0("(", input$tcga_ensg, ")"))
    })
    output$tcga_boxplot <- renderPlotly({
        validate(need(length(input$tcga_tissues) > 0, "Please select at least one tissue."))
        expression_plot(expression_data(), input$tcga_plot, input$tcga_scale,
                        "decreasing", tissue_groups, gene_title())
    })
    selected_expression_plot <- reactive({
        validate(need(length(input$tcga_tissues) > 0, "Please select at least one tissue."))
        expression_plot(expression_data(), input$tcga_plot, input$tcga_scale,
                        input$tcga_order, tissue_groups, gene_title())
    })
    output$tcga_boxplot_with_new_features <- renderPlotly(selected_expression_plot())
    output$download_plotly_widget <- downloadHandler(filename = function() {
        genename <- ensg_genename$gene[which(ensg_genename$ENSG == input$tcga_ensg)]
        if (genename != "") {
            paste("Expitope3_TCGA_", genename, "_", input$tcga_order, ".html", sep = "")
        }
        else {
            paste("Expitope3_TCGA_", input$tcga_ensg, "_", input$tcga_order, ".html", sep = "")
        }
    }, content = function(file) {
        library(pandoc)
        sys_pandoc <- rmarkdown:::find_pandoc()
        saveWidget(as_widget(selected_expression_plot()), file, selfcontained = TRUE)
    })
    output$download_tcga <- downloadHandler(filename = function() paste0(input$tcga_ensg, "_expression.tsv"),
        content = function(file) {
            req(input$tcga_ensg, input$tcga_tissues)
            tables <- lapply(input$tcga_tissues, function(tissue) read_gene_expression(unname(expression_paths[[tissue]]), input$tcga_ensg))
            data.table::fwrite(do.call(rbind, tables), file, sep = "\t")
        }, contentType = "text/tab-separated-values")
    output$tcga_abb <- DT::renderDataTable({
        validate(need(file.exists(file.path(config$data_dir, "tcga_abbreviations")), "Missing tcga_abbreviations dataset"))
        TCGA_abb <- read.csv(file.path(config$data_dir, "tcga_abbreviations"), header = TRUE, sep = "\t")
        DT::datatable(TCGA_abb, extensions = c("ColReorder", "Buttons"), options = list(pageLength = 20,
            columnDefs = list(list(visible = FALSE, targets = c(0)))))
    }, server = FALSE, processing = TRUE, DT.warn.size = FALSE, reactivity = F)
    score_files <- c(`3d_scores` = "score.sc", tcrdock_scores = "lowest_pae_per_run",
        `3d_rosetta_scores` = "score_mut_ros.sc", `3d_ala_scores` = "score_ala.sc",
        `3d_sim_scores` = "score_mut_sim.sc")
    for (output_id in names(score_files)) local({
        id <- output_id
        path <- file.path(config$structure_dir, score_files[[id]])
        output[[id]] <- DT::renderDT({
            validate(need(file.exists(path), paste("Missing score file:", basename(path))))
            scores <- read_structure_scores(path)
            DT::datatable(scores[nrow(scores):1, , drop = FALSE], rownames = FALSE,
                filter = "top", extensions = "Buttons", options = list(scrollX = TRUE,
                    pageLength = 10, dom = "Bfrtip", buttons = list("colvis")))
        }, server = FALSE)
    })
    observeEvent(input$submit_tcr, {
        showNotification("Structure job submission is unavailable until the Linux TCRdock/Rosetta worker is restored and validated.", type = "warning", duration = 10)
    })
    output$running_directory_list_tcrdock <- renderText({
        files <- files_reactive_running_tcrdock()
        if (files == "") {
            paste("No runs in queue, waiting for input")
        }
        else {
            files
        }
    })
    output$finished_directory_list_tcrdock <- renderText({
        files <- files_reactive_finished_tcrdock()
        if (files == "") {
            paste("Waiting for jobs to finish")
        }
        else {
            paste0(files)
        }
    })
    files_reactive_running_tcrdock <- reactivePoll(1000, session,
        checkFunc = function() queue_listing(file.path(config$structure_dir, "todo")),
        valueFunc = function() queue_listing(file.path(config$structure_dir, "todo")))
    output$running_directory_mutate <- renderText({
        files <- files_reactive_running_mutate()
        if (files == "") {
            paste("No runs in queue, waiting for input")
        }
        else {
            files
        }
    })
    files_reactive_running_mutate <- reactivePoll(1000, session,
        checkFunc = function() queue_listing(file.path(config$structure_dir, c("todo_mutate", "todo_alascan", "todo_mutate_rosetta"))),
        valueFunc = function() queue_listing(file.path(config$structure_dir, c("todo_mutate", "todo_alascan", "todo_mutate_rosetta"))))
    output$finished_directory_mutate <- renderText({
        files <- files_reactive_finished_mutate()
        if (files == "") {
            paste("No runs in queue, waiting for input")
        }
        else {
            files
        }
    })
    files_reactive_finished_mutate <- reactive({
        files <- list.files(file.path(config$structure_dir, "mutated"), pattern = "\\.pdb$", ignore.case = TRUE)
        paste(length(files), "mutated PDB files available")
    })
    observeEvent(input$pdb_upload, {
        req(input$pdb_upload)
        name <- input$pdb_upload$name
        validate(need(grepl("^[A-Za-z0-9_.-]+\\.pdb$", name, ignore.case = TRUE), "Use a simple .pdb filename"))
        validate(need(any(grepl("^(ATOM  |HETATM)", readLines(input$pdb_upload$datapath, warn = FALSE))), "No PDB atom records found"))
        target <- file.path(config$structure_dir, "upload", name)
        validate(need(!file.exists(target), "A file with this name already exists; rename your upload"))
        validate(need(file.copy(input$pdb_upload$datapath, target), "Could not save the upload"))
        updateSelectInput(session, "pdb_uploaded", choices = c(list.files(paste0(config$structure_dir,
            "/upload"), full.names = FALSE, pattern = "pdb")), selected = list.files(paste0(config$structure_dir,
            "/upload"), full.names = FALSE, pattern = "pdb")[1])
    })
    output$downloadPDB <- downloadHandler(filename = function() {
        file2 <- input$pdb
    }, content = function(file) {
        if (!file.copy(resolve_pdb_file(input$pdb, config$structure_dir), file, overwrite = TRUE))
            stop("Could not prepare PDB download")
    }, contentType = "pdb")
    observeEvent(input$mutate_button, {
        showNotification("Structure job submission is unavailable until the Linux TCRdock/Rosetta worker is restored and validated.", type = "warning", duration = 10)
    })
    observeEvent(input$mutate_rosetta_button, {
        showNotification("Structure job submission is unavailable until the Linux TCRdock/Rosetta worker is restored and validated.", type = "warning", duration = 10)
    })
    observeEvent(input$ala_button, {
        showNotification("Structure job submission is unavailable until the Linux TCRdock/Rosetta worker is restored and validated.", type = "warning", duration = 10)
    })
    files_reactive_finished_tcrdock <- reactivePoll(1000, session,
        checkFunc = function() list_pdb_files(config$structure_dir),
        valueFunc = function() {
            files <- list_pdb_files(config$structure_dir)
            pdb_results_summary(config$structure_dir)
        })
    updateSelectInput(session, "pdb", choices = list_pdb_files(config$structure_dir))
    updateSelectInput(session, "pdb_uploaded", choices = c(list.files(paste0(config$structure_dir, "/upload"),
        full.names = FALSE, pattern = "pdb")), selected = c(list.files(paste0(config$structure_dir, "/upload"),
        full.names = FALSE, pattern = "pdb"))[1])
    get_epitope_positions <- function(pdbfile = input$pdb) pdb_residue_selection(resolve_pdb_file(pdbfile,
        config$structure_dir), chain = "B")
    tospin <- reactive(input$spin)
    output$structure <- renderNGLVieweR({
        req(input$pdb)
        pdb_to_read <- resolve_pdb_file(input$pdb, config$structure_dir)
        NGLVieweR(pdb_to_read) %>%
            stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
            addRepresentation("cartoon", param = list(name = "cartoon", colorScheme = "chainid")) %>%
            addRepresentation("ball+stick", param = list(name = "peptide_chain_B",
                colorScheme = "element", sele = get_epitope_positions()))

    })
    sele <- reactiveValues()
    observe({
        sele$aa <- str_extract(input$structure_selection, "(?<=[\\[])(.*?)(?=\\])")
        sele$aa_bond <- str_extract(input$structure_selection, "(?<=[\\]])(.*?)(?=[:space:])")
        sele$resiChain <- str_extract(input$structure_selection, "(?<=[]])(.*?)(?=[.])")
        sele$resi <- str_extract(input$structure_selection, "(?<=[]])(.*?)(?=[:])")
        sele$fileName <- str_extract(input$structure_selection, "(?<=[(])(.*?)(?=[.])")
    })
    output$selection = renderPrint({
        print(input$structure_selection)
        print(sele$aa)
        print(sele$aa_bond)
        print(sele$resiChain)
        print(sele$resi)
        print(sele$fileName)
        print(input$structure_selAround)
    })
    observeEvent(input$structure_selAround, {
        NGLVieweR_proxy("structure") %>% removeSelection("selAround")
        NGLVieweR_proxy("structure") %>% addSelection("ball+stick", param = list(name = "selAround",
            sele = input$structure_selAround, colorValue = "grey"))
    })
    observeEvent(sele$resiChain, {
        NGLVieweR_proxy("structure") %>% removeSelection("label")
        NGLVieweR_proxy("structure") %>% removeSelection("contact")
        NGLVieweR_proxy("structure") %>% addSelection("label", param = list(name = "label", sele = sele$resiChain,
            labelType = "format", labelFormat = "[%(resname)s]%(resno)s", labelGrouping = "residue",
            color = "white", xOffset = 1, fixedSize = TRUE, radiusType = 1, radiusSize = 1.5)) %>% addSelection("contact",
            param = list(name = "contact", sele = "*", filterSele = list(sele$resiChain, "*"), labelVisible = TRUE,
                labelFixedSize = FALSE, labelUnit = "angstrom", labelSize = 2))
    })
    observeEvent(input$animate, {
        if (input$animate == "Rock") {
            NGLVieweR_proxy("structure") %>% updateRock(TRUE)
        }
        else if (input$animate == "Spin") {
            NGLVieweR_proxy("structure") %>% updateSpin(TRUE)
        }
        else {
            NGLVieweR_proxy("structure") %>% updateRock(FALSE) %>% updateSpin(FALSE)
        }
    })
    observeEvent(input$add, {
        NGLVieweR_proxy("structure") %>% addSelection(isolate(input$type), param = list(name = "sel1",
            sele = isolate(input$selection), colorValue = isolate(input$color)))
    })
    observeEvent(input$remove, {
        NGLVieweR_proxy("structure") %>% removeSelection("sel1")
    })
    data_tcrinput <- reactive({
        req(input$rows)
        count <- suppressWarnings(as.integer(input$rows))
        validate(need(length(count) == 1 && !is.na(count) && count >= 1 && count <= 100,
                      "Enter a row count between 1 and 100"))
        do.call(rbind, base::replicate(count, tcr_dock[1, , drop = FALSE], simplify = FALSE))
    })
    dfon <- reactiveValues(data = NULL)
    observe({
        dfon$data <- data_tcrinput()
    })
    output$table3 <- renderDataTable({
        datatable(data_tcrinput(), editable = "cell", options = list(dom = "t"))
    })
    observeEvent(input$table3_cell_edit, {
        info <- input$table3_cell_edit
        i <- info$row
        j <- info$col
        v <- info$value
        dfon$data[i, j] <- v
    })
}
