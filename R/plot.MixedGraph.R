
darken_color <- function(color) {
    rgb(t(col2rgb(color)/ 1.2), maxColorValue=255)
}

produce_colors <- function(K) {
    colors <- rainbow(K)
    substr(colors, start = 1, stop = 7)
}

#' Plot function for S3 class "MixedGraph"
#'
#' @description Plot the graph of MixedGraph object using igraph, cytoscape and cytpscape.js.
#'
#' @method plot MixedGraph 
#'
#' @param x is a MixedGraph object.
#' @param method is the related package or lib used in the function. When "igraph", the function will use the R package "igraph". When "cytoscape", the function will use the R package RCy3, the user should start the Cytoscape software before call the function. And when "cytoscape.js", the function will use the lib cytoscape.js.
#' @param weighted is a boolean value, which indicates if we would plot the width of edge according to the weight of edge.
#' @param stability is the stability threshold, when the stability score of the edge is more than it means the coefficient between the two vertices can be trusted. 
#' @param out.file is the file name to save the plot of MixedGraph object. The default value is NULL, and the graph will be plotted to the screen.
#' @param ... other arguments for different methods.
#' \itemize{
#' \item{igraph}{the generic arguments for plot.igraph is available.}
#' \item{cytoscape}{layout: the names should be in RCy3::getLayoutNames(CytoscapeWindow), the default layout is "attributes-layout". And you can also modify the layout throungh Cytoscape software.}
#' \item{cytoscape.js}{To be extend}
#' }
#'
#' @examples
#' X1 <- matrix(rnorm(12), nrow = 4)
#' X2 <- matrix(rnorm(12), nrow = 4)
#' X <- list(X1, X2)
#' crf_structure <- matrix(c(1, 1, 0, 1), nrow = 2)
#' brail_control <- list(B = 5, tau = 0.6)
#' G <- MixedGraph(X, crf_structure, brail_control = brail_control)
#' plot(G, method = "igraph", weighted = TRUE)
#'
#' @importFrom graphics plot
#' @import igraph
#' @import grDevices
#' @import htmlwidgets
#' @import RCy3
#' @export


plot.MixedGraph <- function(x, method = c("igraph", "cytoscape", "cytoscape.js"), weighted = FALSE, stability = 0.0, out.file = NULL, ...) {
    method <- tolower(method)
    method <- match.arg(method)
    K <- length(x$data)
    size_list <- cumsum(c(1, sapply(x$data, ncol)))
    p <- ncol(x$network)

    # produce color list
    colors <- produce_colors(K)
    graph_color <- as.vector(sapply(seq_along(x$data), function(i){
        rep(colors[i], ncol(x$data[[i]]))
    }))
    graph_color <- unlist(graph_color)

    # construct the directed network
    network <- x$network
    indexes <- (x$stability >= stability)
    network[!indexes] <- 0
    directed_network <- network
    sapply(1:K, function(i) {
        indexes <- size_list[i] : (size_list[i + 1] - 1)
        directed_network[indexes, indexes] <<- 0
    })

    # construct the undirected network
    undirected_network <- network
    undirected_network[directed_network!=0] <- 0.0
    if (x$rule == "AND") {
        undirected_network[!(undirected_network & base::t(undirected_network))] <- 0
    }
    else {
        undirected_network[!(undirected_network | base::t(undirected_network))] <- 0
    }
    undirected_network <- (undirected_network + base::t(undirected_network)) / 2

    # set the ids and get the label names from the input data
    ids <- 1 : length(graph_color)
    labelnames <- NA
    if (!is.null(colnames(x$data[[1]]))) {
        labelnames <- as.vector(sapply(x$data, function(x) {
            colnames(x)
        }))
    }

    if (method == "igraph") {
        # set the arrow size of directed edges
        argv_list <- list(...)
        directed_graph <- graph.adjacency(directed_network, weighted = TRUE, mode = "directed")
        V(directed_graph)$name <- ids
        directed_arrow_size <- 0.3
        if ("arrow.size" %in% names(argv_list)) {         
            directed_arrow_size <- argv_list["arrow.size"]
            argv_list["arrow.size"] <- NULL  
        }

        # set the color of undirected edges
        undirected_graph <- graph.adjacency(undirected_network, weighted = TRUE, mode = "directed")
        V(undirected_graph)$name <- ids
        edge_start <- ends(undirected_graph, es=E(undirected_graph), names=F)[,1]
        undirected_edge_color <- graph_color[edge_start]

        # set the edge width according the edge weight
        if(weighted){
            E(directed_graph)$width <- 0.5 + abs(E(directed_graph)$weight)
            E(undirected_graph)$width <- 0.5 + abs(E(undirected_graph)$weight)
        }
        argv_list <- c(argv_list, list(vertex.color = graph_color, vertex.label = labelnames, 
                                       vertex.size = 30 / as.integer(vcount(directed_graph)/100 + 1) ))

        # if user didn't set the layout, set the default layout
        if (!"layout" %in% names(argv_list)) {
            weight_network <- matrix(0, p, p)
            weight_network[directed_network != 0] <- 1
            sapply(1:K, function(i) {
                indexes <- size_list[i] : (size_list[i + 1] - 1)
                weight_network[indexes, indexes] <<- 1})
            weight_graph <- graph.adjacency(weight_network, weighted = TRUE)
            layout <- layout.fruchterman.reingold(weight_graph, weights=E(weight_graph)$weight)
            argv_list <- c(argv_list, list(layout = layout))
        }

        if(is.null(out.file) == FALSE){
            pdf(out.file)
            do.call(plot.igraph, c(list(x = directed_graph, edge.arrow.size = directed_arrow_size), argv_list))
            do.call(plot.igraph, c(list(x = undirected_graph, add = T, edge.arrow.size = 0, edge.color = undirected_edge_color), argv_list))
            dev.off()
            cat(paste("Output file: ", out.file, "\n", sep=""))
        }
        else {
            do.call(plot.igraph, c(list(x = directed_graph, edge.arrow.size = directed_arrow_size), argv_list))
            do.call(plot.igraph, c(list(x = undirected_graph, add = T, edge.arrow.size = 0, edge.color = undirected_edge_color), argv_list))
        }
    }
    else if (method == "cytoscape") {
        undirected_network[lower.tri(undirected_network)] <- 0
        mixed_network <- directed_network + undirected_network
        mixed_graph <- graph.adjacency(mixed_network, weighted = TRUE, mode = "directed")
        g <- igraph.to.graphNEL(mixed_graph)
        
        # set node group
        g <- initNodeAttribute (graph = g, attribute.name = 'group', 
                                attribute.type = 'integer', 
                                default.value = 0)
        count <- 0
        index <- 1
        for (node in nodes(g)) {
            count <- count + 1
            if (count >= size_list[index]) index <- index + 1
            graph::nodeData(g, node, 'group') <- index - 1
        }

        # set edge type attribute
        g <- initEdgeAttribute (graph = g,  attribute.name = 'edgeType',
                                attribute.type ='char',
                                default.value = "undefined")
        g <- initEdgeAttribute (graph = g,  attribute.name = 'weight',
                                attribute.type ='numeric',
                                default.value = 0)
        undirected_indexes <- which(undirected_network != 0, arr.ind = T)
        for (i in 1 : nrow(undirected_indexes)) {
            node <- formatC(undirected_indexes[i,"row"])
            group <- graph::nodeData(g, node, 'group')
            graph::edgeData(g, from = node, to = formatC(undirected_indexes[i,"col"]), 'edgeType') <- paste("Type", group, sep = "")
        }

        # display cytoscape windows
        cw <- CytoscapeWindow ('test', graph = g, overwriteWindow = TRUE)
        displayGraph(cw)
        argv_list <- list(...)
        if ("layout" %in% names(argv_list)) {
            if (!argv_list["layout"] %in% getLayoutNames(cw)) {
                stop("The input layout not in getLayoutNames(CytoscapeWindow)!")
            }
            layout <- argv_list["layout"]
        }
        else {
            layout <- "attributes-layout"
        }
        layoutNetwork(cw, layout.name = layout)
        
        # set rules
        setDefaultNodeShape (cw, 'ELLIPSE')
        setDefaultNodeSize  (cw, 35)
        setDefaultNodeFontSize (cw, 10)
        setNodeColorRule(cw, 'group', c(1 : K), colors, mode = 'lookup')
        setEdgeTargetArrowRule(cw, 'edgeType', paste('Type', c(1 : K), sep = ""), rep('None', K), default = 'Arrow')
        setEdgeTargetArrowColorRule(cw, 'edgeType', paste('Type', c(1 : K), sep = ""), colors, mode = 'lookup')
        setEdgeColorRule(cw, 'edgeType', paste('Type', c(1 : K), sep = ""), colors, mode  = 'lookup', default.color = '#000000')

        if(is.null(out.file) == FALSE) {
            if (nchar(out.file) < 4) stop("Out.file should match one of the type c('png', 'pdf', 'svg')")
            image_type <- substr(out.file, nchar(out.file) - 2, nchar(out.file))
            file_name <- substr(out.file, 1, nchar(out.file) - 4)
            match.arg(image_type, c("png", "pdf", "svg"))
            redraw (cw)
            saveImage(cw, file_name, image.type = image_type)
        }
    }
    else if (method == "cytoscape.js") {
        # set node entries
        nodes <- data.frame(id = ids, name = labelnames, color = graph_color)
        node_entries <- apply(nodes, 1, function(x) {
            x[1] = trimws(x[1])
            x[2] = trimws(x[2])
            list(data = as.list(x))
        })

        # set edge entries
        directed_indexes <- which(directed_network != 0, arr.ind=T)
        undirected_network[lower.tri(undirected_network)] <- 0
        undirected_indexes <- which(undirected_network != 0, arr.ind=T)
        if (weighted) {
            edges <- data.frame(source = directed_indexes[,"row"], target = directed_indexes[,"col"], 
                                weight = directed_network[directed_indexes], directed = rep(TRUE, nrow(directed_indexes)))
            edges <- rbind(edges, data.frame(source = undirected_indexes[,"row"], target = undirected_indexes[,"col"], 
                                weight = undirected_network[undirected_indexes], directed = rep(FALSE, nrow(undirected_indexes))))
        }
        else {
            edges <- data.frame(source = directed_indexes[,"row"], target = directed_indexes[,"col"], 
                                directed = rep(TRUE, nrow(directed_indexes)))
            edges <- rbind(edges, data.frame(source = undirected_indexes[,"row"], target = undirected_indexes[,"col"], 
                                directed = rep(FALSE, nrow(undirected_indexes))))
        }
        edges_entries <- apply(edges, 1, function(x) {
            list(data = as.list(x))
        })
        cy <- list(nodes = node_entries, edges = edges_entries)

        # save file
        if(is.null(out.file) == FALSE) {
            htmlwidgets::saveWidget(Cytoscapejs(cy), file = out.file, selfcontained = TRUE)
            cat(paste("Output file: ", out.file, "\n", sep=""))
        }
        else {
            Cytoscapejs(cy)
        }
    }
}