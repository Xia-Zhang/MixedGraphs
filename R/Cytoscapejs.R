Cytoscapejs <- function(cy, width = NULL, height = NULL, elementId = NULL) {

  # forward options using x
  x <- cy
  # create widget
  htmlwidgets::createWidget(
    name = 'Cytoscapejs',
    x,
    width = width,
    height = height,
    package = 'MixedGraphs',
    elementId = elementId
  )
}

#' Shiny bindings for Cytoscapejs
#'
#' Output and render functions for using Cytoscapejs within Shiny
#' applications and interactive Rmd documents.
#'
#' @param outputId output variable to read from
#' @param width,height Must be a valid CSS unit (like \code{'100\%'},
#'   \code{'400px'}, \code{'auto'}) or a number, which will be coerced to a
#'   string and have \code{'px'} appended.
#' @param expr An expression that generates a Cytoscapejs
#' @param env The environment in which to evaluate \code{expr}.
#' @param quoted Is \code{expr} a quoted expression (with \code{quote()})? This
#'   is useful if you want to save an expression in a variable.
#'
#' @name Cytoscapejs-shiny
#'
CytoscapejsOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'Cytoscapejs', width, height, package = 'MixedGraphs')
}

#' @rdname Cytoscapejs-shiny
renderCytoscapejs <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, CytoscapejsOutput, env, quoted = TRUE)
}
