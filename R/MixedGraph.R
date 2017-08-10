
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
