guess_family <- function(X) {
    if (is.list(X)) {
        result <- sapply(X, function(x){
            guess_family(x)
        })
        return(result)
    }
    if (all(is.element(X, c(0, 1)))) {
        return("binomial")
    }
    else if (all(X == as.integer(X))) {
        return("poisson")
    }
    else {
        return("gaussian")
    }
}

produce_colors <- function(K) {
    colors <- c("#8dd3c7", "#fb8072", "#ffffb3", "#bebada", "#80b1d3", "#fdb462", "#b3de69", "#fccde5", "#d9d9d9", "#bc80bd", "#ccebc5", "#ffed6f")
    if (K > 12)
        colors <- c(colors, rainbow(K - 12, alpha=.8))
    colors
}

#' MixedGraph is used to fit mixed graph model.
#'
#' @param X is a list containing k matrices, each of the matrices has n rows and pk columns. They are the 'blocks'. Each block should have all columns of the same type.
#' @param crf_structure is a K-by-K matrix giving the structure of the CRF
#'  \itemize{
#'  \item{element (i, j) is 1}{if element (i, j) is 1 if there are directed edges going from the i-th block to the j-th block;}
#'  \item{element (i, j) is 0}{if there are no edges going from the i-th block to the j-th block;}
#'  \item{element (i, i) is 1}{if there are undirected edges within the i-th block;}
#'  \item{element (i, i) is 0}{if there are not undirected edges within the i-th block}
#' }
#' @param family is a K vector giving the data types for each blocks. The type could be "gaussian", "binomial" or "poisson". If NULL, we can use a heuristic to guess the right value
#' @param rule says whether to use the "AND" or "OR" rule
#' @param brail_control is a list of optional arguments to be passed to BRAIL internally
#'
#' @return the coefficients vector
#'  \item{data}{ is a list containing k matrices from input}
#'  \item{network}{ is a p * p (sparse) matrix giving the edge weights returned by BRAIL}
#'  \item{family}{ is a K vector giving the data types for each blocks from input, if the input family is NULL, there will be the values guessed in the MixedGraph function}
#'  \item{crf_structure}{ is a K-by-K matrix giving the structure of the CRF from input}
#'  \item{stability}{ is a p * p (sparse) matrix giving the stability scores}
#'  \item{rule}{"AND" or "OR" rule for edge selection within blocks}
#'
#' @examples
#' X1 <- matrix(rnorm(12), nrow = 4)
#' X2 <- matrix(rnorm(12), nrow = 4)
#' X <- list(X1, X2)
#' crf_structure <- matrix(rep(1, 4), nrow = 2)
#' brail_control <- list(B = 5, tau = 0.6)
#' MixedGraph(X, crf_structure, brail_control = brail_control)
#'
#' @importFrom foreach %dopar% %:%
#' @export


MixedGraph <-function(X, crf_structure, family = NULL, rule = c("AND", "OR"), brail_control = NULL) {
    K <- length(X)
    p <- sum(sapply(X, ncol))
    if (is.null(family) == FALSE && length(family) != K) {
        stop("Input family error!")
    }
    if (is.null(family)) {
        family <- guess_family(X)
    }
    rule <- toupper(rule)
    rule <- match.arg(rule)

    ncols <- sapply(X, ncol)
    size_list <- cumsum(c(1, ncols))
    indexes_mask <- matrix(unlist(sapply(1:K, function(i){
        rep(rep(crf_structure[,i], ncols), ncols[i])
    })), nrow = p)
    diag(indexes_mask) <- 0

    k <- 1
    i <- 1
    graph_list <- foreach::foreach(k = 1:K, .packages='MixedGraphs') %:% 
        foreach::foreach(i = 1:ncol(X[[k]])) %dopar% {
            crf_k <- crf_structure[,k]
            block_indexes <- which(crf_k == 1)
            if (length(block_indexes)==0) {
                return(list(coef = rep(0, p), score = rep(0, p)))
            }
            brail_X <- X[block_indexes]
            index_k <- match(k, block_indexes) # whether block has undirected edges
            if (is.na(index_k) == FALSE) {
                brail_X[[index_k]] <- X[[k]][,-i]
            }
            brail_y <- X[[k]][,i]
            brail_family <- family[k]
            brail_argv <- list(X = brail_X, y = brail_y, family = brail_family, doPar = FALSE)
            brail_res <- do.call(BRAIL, c(brail_argv, brail_control))
            
            coef <- rep(0, p)
            coef[which(indexes_mask[,size_list[k] + i - 1] == 1)] <- unlist(brail_res$coefficients)

            score <- rep(0, p)
            score[which(indexes_mask[,size_list[k] + i - 1] == 1)] <- unlist(brail_res$score)

            list(coef = coef, score = score)
        }

    result <- list()
    result$data <- X
    result$network <- matrix(unlist(sapply(graph_list, function(x) {
        tmp_list <- sapply(x, function(subx) {
            unlist(subx$coef)
        })
    })), p, p)
    result$family <- family
    result$crf_structure <- crf_structure
    result$stability <- matrix(unlist(sapply(graph_list, function(x) {
        tmp_list <- sapply(x, function(subx) {
            unlist(subx[[2]])
        })
    })), p, p)
    result$rule <- rule
    class(result) <- "MixedGraph"
    result
}


#' Plot the graph from Mixedgraph object.
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
        undirected_edge_color <- rgb(t(col2rgb(graph_color[edge_start])/ 1.2), maxColorValue=255)

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
        # save the graph in cytoscape
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