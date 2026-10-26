upsideDown <- function(data) {
    new <- data[order(rev(1:dim(data)[1])), ]
    new
}

ggp <- function(obj, info) {
    g <- ggplotly(obj, tooltip = info)
    d <- g$x$data
    for (i in 1:length(d)) {
        w <- d[[i]]$text
        if (length(w))
            g$x$data[[i]]$text <- gsub("info: ", "", w)
    }
    g
}

ggp800 <- function(obj, info) {
    g <- ggplotly(obj, tooltip = info, height = 1000)
    d <- g$x$data
    for (i in 1:length(d)) {
        w <- d[[i]]$text
        if (length(w))
            g$x$data[[i]]$text <- gsub("info: ", "", w)
    }
    g
}

ggp700 <- function(obj, info) {
    g <- ggplotly(obj, tooltip = info, height = 700)
    d <- g$x$data
    for (i in 1:length(d)) {
        w <- d[[i]]$text
        if (length(w))
            g$x$data[[i]]$text <- gsub("info: ", "", w)
    }
    g
}

ggp700_source <- function(obj, info, source) {
    g <- ggplotly(obj, tooltip = info, height = 700, source = source)
    d <- g$x$data
    for (i in 1:length(d)) {
        w <- d[[i]]$text
        if (length(w))
            g$x$data[[i]]$text <- gsub("info: ", "", w)
    }
    g
}

ggp2 <- function(obj, info) {
    g <- ggplotly(obj, tooltip = info)
    d <- g$x$data
    for (i in 1:length(d)) {
        w <- d[[i]]$text
        g$x$data[[i]]$text[length(g$x$data[[i]]$text) + 1] <- c("")
        if (length(w))
            g$x$data[[i]]$text <- gsub("info: ", "", w)
    }
    g
}

js <- c("function(el,x){", "  $('[id^=Color]').on('change', function(){", "    var color = this.value;",
    "    var id = this.id;", "    var index = parseInt(id.split('-')[1]) - 1;", "    var data = el.data;",
    "    var marker = data[index].marker;", "    marker.color = color;", "    marker.size = 15;", "    Plotly.restyle(el, {marker: marker}, [index]);",
    "  });", "}")

radioTooltip <- function(id, choice, title, placement = "bottom", trigger = "hover", options = NULL) {
    options = shinyBS:::buildTooltipOrPopoverOptionsList(title, placement, trigger, options)
    options = paste0("{'", paste(names(options), options, sep = "': '", collapse = "', '"), "'}")
    bsTag <- shiny::tags$script(shiny::HTML(paste0("\n    $(document).ready(function() {\n      setTimeout(function() {\n        $('input', $('#",
        id, "')).each(function(){\n          if(this.getAttribute('value') == '", choice, "') {\n            opts = $.extend(",
        options, ", {html: true});\n            $(this.parentElement).tooltip('destroy');\n            $(this.parentElement).tooltip(opts);\n          }\n        })\n      }, 500)\n    });\n  ")))
    htmltools::attachDependencies(bsTag, shinyBS:::shinyBSDep)
}

options(DT.warn.size = FALSE)
