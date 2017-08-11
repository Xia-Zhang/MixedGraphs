
darken_color <- function(color) {
    rgb(t(col2rgb(color)/ 1.2), maxColorValue=255)
}

produce_colors <- function(K, node = TRUE) {
    colors <- rainbow(K)
    print(colors <- rainbow(K))
    
    node_colors <- c("#8dd3c7", "#fb8072", "#ffffb3", "#bebada", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
    edge_colors <- node_colors
    colors <- node_colors
    if (!node) {
        colors <- edge_colors
    }
    if (K > length(node_colors)) {
        extra_color <- rainbow(K - 12, alpha=.8)
        if (node) {
            colors <- c(node_colors, extra_color)
        }
        else {
            extra_color <- vapply(extra_color, darken_color)
            colors <- c(edge_colors, extra_color)
        }
    }
    colors[1 : K]
}

#' Plot the graph from MixedGraph object.
#' 
#' @method plot MixedGraph 
#'
#' @param x is a MixedGraph object.
#' @param method is the package or the return type used in the function. When is "igraph", use the package "igraph". When is "cytoscape", record the network of graph modeling language (GML) format, which can be imported in cytoscape.Or we can use the R package RCytoscape. I haven't decided now.
#' @param weighted is a boolean value, which indicate if we would plot the width of edge according to the weight of edge.
#' @param stability is the stability threshold, more than it indicate the coefficient between the two vertices can be trusted. 
#' @param out.file is the file name to save the plot of MixedGraph object. The default value is NULL, and the graph will be plotted to the screen.
#' @param ... other generic arguments for plot method
#'
#' @examples
#' X1 <- matrix(rnorm(12), nrow = 4)
#' X2 <- matrix(rnorm(12), nrow = 4)
#' X <- list(X1, X2)
#' crf_structure <- matrix(rep(1, 4), nrow = 2)
#' brail_control <- list(B = 2, tau = 0.6)
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

    colors <- produce_colors(K)
    graph_color <- as.vector(sapply(seq_along(x$data), function(i){
        rep(colors[i], ncol(x$data[[i]]))
    }))
    graph_color <- unlist(graph_color)

    network <- x$network
    indexes <- (x$stability >= stability)
    network[!indexes] <- 0
    directed_network <- network
    sapply(1:K, function(i) {
        indexes <- size_list[i] : (size_list[i + 1] - 1)
        directed_network[indexes, indexes] <<- 0
    })

    undirected_network <- network
    undirected_network[directed_network!=0] <- 0.0
    if (x$rule == "AND") {
        undirected_network[!(undirected_network & base::t(undirected_network))] <- 0
    }
    else {
        undirected_network[!(undirected_network | base::t(undirected_network))] <- 0
    }
    undirected_network <- (undirected_network + base::t(undirected_network)) / 2

    ids <- 1 : length(graph_color)
    labelnames <- NA
    if (!is.null(colnames(x$data[[1]]))) {
        labelnames <- as.vector(sapply(x$data, function(x) {
            colnames(x)
        }))
    }

    if (method == "igraph") {
        argv_list <- list(...)
        directed_graph <- graph.adjacency(directed_network, weighted = TRUE, mode = "directed")
        V(directed_graph)$name <- ids
        directed_arrow_size <- 0.3
        if ("arrow.size" %in% names(argv_list)) {         
            directed_arrow_size <- argv_list["arrow.size"]
            argv_list["arrow.size"] <- NULL  
        }

        undirected_graph <- graph.adjacency(undirected_network, weighted = TRUE, mode = "directed")
        V(undirected_graph)$name <- ids
        # undirected edge color using start node color, and darker the color
        edge_start <- ends(undirected_graph, es=E(undirected_graph), names=F)[,1]
        undirected_edge_color <- darken_color(graph_color[edge_start])

        if(weighted){
            E(directed_graph)$width <- 0.5 + abs(E(directed_graph)$weight)
            E(undirected_graph)$width <- 0.5 + abs(E(undirected_graph)$weight)
        }
        argv_list <- c(argv_list, list(vertex.color = graph_color, vertex.label = labelnames, 
                                       vertex.size = 20 / as.integer(vcount(directed_graph)/100 + 1) ))

        if (!"layout" %in% names(argv_list)) {
            weight_network <- matrix(1, p, p)
            sapply(1:K, function(i) {
                indexes <- size_list[i] : (size_list[i + 1] - 1)
                weight_network[indexes, indexes] <<- 10})
            weight_graph <- graph.adjacency(weight_network, weighted = TRUE)
            layout <- layout.fruchterman.reingold(weight_graph, weights=E(weight_graph)$weight)
            argv_list <- c(argv_list, list(layout = layout))
        }

        if(is.null(out.file) == FALSE){
            pdf(out.file)
            do.call(plot.igraph, c(list(x = directed_graph, edge.arrow.size = directed_arrow_size), argv_list))
            do.call(plot.igraph, c(list(x = undirected_graph, add = T, edge.arrow.size = 0), argv_list))
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

        # set edge type attributes
        g <- initEdgeAttribute (graph = g,  attribute.name = 'edgeType',
                                attribute.type ='char',
                                default.value = "undefined")
        g <- initEdgeAttribute (graph = g,  attribute.name = 'weight',
                                attribute.type ='numeric',
                                default.value = 0)
        undirected_indexes <- which(undirected_network != 0, arr.ind = T)
        print(undirected_indexes)
        for (i in 1 : nrow(undirected_indexes)) {
            node <- formatC(undirected_indexes[i,"row"])
            group <- graph::nodeData(g, node, 'group')
            print(group)
            graph::edgeData(g, from = node, to = formatC(undirected_indexes[i,"col"]), 'edgeType') <- paste("Type", group, sep = "")
        }

        # display cytoscape windows
        cw <- CytoscapeWindow ('test', graph = g, overwriteWindow = TRUE)
        displayGraph(cw)
        layoutNetwork(cw, layout.name = "attributes-layout")
        
        # set rules
        setDefaultNodeShape (cw, 'ELLIPSE')
        setDefaultNodeSize  (cw, 35)
        setDefaultNodeFontSize (cw, 10)
        setNodeColorRule(cw, 'group', c(1 : K), produce_colors(K), mode = 'lookup')
        setEdgeTargetArrowRule(cw, 'edgeType', paste('Type', c(1 : K), sep = ""), rep('None', K), default = 'Arrow')
        setEdgeTargetArrowColorRule(cw, 'edgeType', paste('Type', c(1 : K), sep = ""), produce_colors(K, FALSE), mode = 'lookup')
        setEdgeColorRule(cw, 'edgeType', paste('Type', c(1 : K), sep = ""), produce_colors(K, FALSE), mode  = 'lookup', default.color = '#000000')
    }
    else if (method == "cytoscape.js") {
        nodes <- data.frame(id = ids, name = labelnames, color = graph_color)
        node_entries <- apply(nodes, 1, function(x) {
            x[1] = trimws(x[1])
            x[2] = trimws(x[2])
            list(data = as.list(x))
        })

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
        Cytoscapejs(cy)
    }
}