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
    graph_list <- foreach::foreach(k = 1:K, .packages='MixedGraphs') %:% 
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
            coef <- sapply(crf_k, function(tmp_k) {
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
            score <- sapply(crf_k, function(tmp_k) {
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
            list(as.vector(coef), as.vector(score))
        }
    result <- list()
    result$data <- X
    tmp_graph <- sapply(graph_list, function(x) {
        tmp_list <- sapply(x, function(subx) {
            subx[[1]]
        })
    })
    result$graph <- matrix(tmp_graph, nrow = p)
    result$family <- family
    result$crf_structure <- crf_structure
    tmp_stability <- sapply(graph_list, function(x) {
        tmp_list <- sapply(x, function(subx) {
            subx[[2]]
        })
    })
    result$stability <- matrix(tmp_stability, nrow = p)
    class(result) <- "MixedGraph"
    result
}