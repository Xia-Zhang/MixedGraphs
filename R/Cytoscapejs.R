Cytoscapejs <- function(cy, width = NULL, height = NULL, elementId = NULL) {
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

CytoscapejsOutput <- function(outputId, width = '100%', height = '400px'){
  htmlwidgets::shinyWidgetOutput(outputId, 'Cytoscapejs', width, height, package = 'MixedGraphs')
}

renderCytoscapejs <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) } # force quoted
  htmlwidgets::shinyRenderWidget(expr, CytoscapejsOutput, env, quoted = TRUE)
}
