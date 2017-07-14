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
#'
#' @examples
#' X1 <- matrix(rnorm(12), nrow = 4)
#' X2 <- matrix(rnorm(12), nrow = 4)
#' X <- list(X1, X2)
#' crf_structure <- matrix(rep(1, 4), nrow = 2)
#' brail_control <- list(B = 10)
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
    rule <- match.arg(rule)
    k <- 1
    i <- 1
    graph_list <- foreach::foreach(k=1:K, .packages='MixedGraphs') %:% 
        foreach::foreach(i = 1:ncol(X[[k]])) %dopar% {
            crf_k <- crf_structure[,k]
            block_indexes <- which(crf_k == 1)
            brail_X <- X[block_indexes]
            index_k <- match(k, block_indexes)
            if (is.na(index_k) == FALSE) {
                brail_X[[index_k]] <- X[[k]][,-i]
            }
            brail_y <- X[[k]][,i]
            brail_family <- family[k]
            brail_argv <- list(X = brail_X, y = brail_y, family = brail_family, doPar = FALSE)
            brail_res <- do.call(BRAIL, c(brail_argv, brail_control))
            index <- 0
            coef <- lapply(crf_k, function(tmp_k) {
                if (tmp_k == 0) {
                    rep(0, length(X[crf_k]))
                }
                else {
                    index <<- index + 1
                    if (block_indexes[index] == k) {
                        append(brail_res$coefficients[[index]], 0, after = i - 1)
                    }
                    else {
                        brail_res$coefficients[[index]]
                    }
                }
            })
            index <- 0
            score <- lapply(crf_k, function(tmp_k) {
                if (tmp_k == 0) {
                    rep(0, length(X[crf_k]))
                }
                else {
                    index <<- index + 1
                    if (block_indexes[index] == k) {
                        append(brail_res$score[[index]], 0, after = i - 1)
                    }
                    else {
                        brail_res$score[[index]]
                    }
                }
            })
            list(unlist(coef), unlist(score))
        }
    result <- list()
    result$data <- X
    tmp_graph <- sapply(graph_list, function(x) {
        tmp_list <- sapply(x, function(subx) {
            unlist(subx[[1]])
        })
    })
    result$network <- matrix(unlist(tmp_graph), nrow = p)
    result$family <- family
    result$crf_structure <- crf_structure
    tmp_stability <- sapply(graph_list, function(x) {
        tmp_list <- sapply(x, function(subx) {
            unlist(subx[[2]])
        })
    })
    result$stability <- matrix(unlist(tmp_stability), nrow = p)
    class(result) <- "MixedGraph"
    result
}


#' Plot the graph from Mixedgraph object.
#' 
#' @param x is a MixedGraph object.
#' @param method is the package or the return type used in the function. When is "igraph", use the package "igraph". When is "cytoscape", record the network of graph modeling language (GML) format, which can be imported in cytoscape.Or we can use the R package RCytoscape. I haven't decided now.
#' @param weighted is a boolean value, which indicate if we would plot the width of edge according to the weight of edge.
#' @param thresh is the threshold which decide if the two vertices is connected according to the coefficient value. The default value is 1e-6.
#' @param stability is the stability threshold, more than it indicate the coefficient between the two vertices can be trusted. 
#' @param save.fn is the file name to save the plot of MixedGraph object. The default value is NULL, and the graph will be plotted to the screen.
#' @param ... other generic arguments for plot method
#'
#' @examples
#' X1 <- matrix(rnorm(12), nrow = 4)
#' X2 <- matrix(rnorm(12), nrow = 4)
#' X <- list(X1, X2)
#' crf_structure <- matrix(rep(1, 4), nrow = 2)
#' brail_control <- list(B = 10)
#' G <- MixedGraph(X, crf_structure, brail_control = brail_control)
#' plot(G, method = "igraph", weighted = TRUE)
#' 
#' @importFrom grDevices dev.off pdf rainbow
#' @importFrom graphics plot
#' @export


plot.MixedGraph <- function(x, method = c("igraph", "cytoscape", "cytoscape.js"), weighted = FALSE, thresh = 1e-6, stability = 1e-6, save.fn = NULL, ...) {
    method <- tolower(method)
    method <- match.arg(method)
    colors <- rainbow(length(x$data), alpha=.7)
    if (method == "igraph") {
        network <- x$network
        indexes <- (abs(network) > thresh & x$stability > stability)
        # network[indexes] <- 1
        network[!indexes] <- 0
        size_list <- cumsum(c(1, sapply(x$data, ncol)))
        graph_color <- as.vector(sapply(seq_along(colors), function(i){
            rep(colors[i], ncol(x$data[[i]]))
        }))
        graph <- igraph::graph.adjacency(network, weighted = TRUE, mode = "directed")
        # graph_list <- lapply(seq_along(x$data), function(i){
        #     sub_network <- network[size_list[i]: size_list[i+1] - 1, size_list[i]: size_list[i+1] - 1]
        #     graph_sub <- igraph::graph.adjacency(sub_network, "directed")
        #     V(graph_sub)$color <- colors[i]
        #     if(is.null(x$data) == FALSE)
        #         V(graph_sub)$label <- colnames(x$data[[i]])

        #     graph_sub
        # })
        # %u%
        #rgb(10, 100, 100, maxColorValue=255)
        # graph <- graph.union(graph_list)
        # if(is.null(x$data) == FALSE)
        #     V(graph_sub)$label <- colnames(x$data[[i]])
        if(weighted)
            igraph::E(graph)$width <- 0.5 + abs(igraph::E(graph)$weight)
        igraph::V(graph)$color <- unlist(graph_color)
        igraph::V(graph)$size <- 40
        igraph::E(graph)$arrow.size <- .5
    }
    else if (method == "cytoscape") {
        # save the graph in cytoscape
    }
    else if (method == "cytoscape.js") {

    }
    if(is.null(save.fn) == FALSE){
        pdf(save.fn)
        plot(graph)
        dev.off()
        cat(paste("Output file: ", save.fn, "\n",sep=""))
    }
    else {
        plot(graph)
    }
}