
get_c <- function(X_k, beta_k, family) {
    n <- nrow(X_k)
    W <- matrix(0, nrow = n, ncol = n)
    if (family == "gaussian") {
        return(1.0)
    }
    else if(family == "logistic") {
        for (i in 1:n) {
            tmp <- 1 /(1 + exp(-(X_k[i,] %*% beta_k)))
            W[i, i] <- tmp * (1 - tmp) / n
        }
    }
    else if(family == "poisson") {
        for (i in 1:n) {
            W[i, i] <- exp(X_k[i,] %*% beta_k) / n
        }
    }
    return(max(eigen(t(X_k) %*% W %*% X_k)$values) / max(eigen(t(X_k) %*% X_k)$values))
}

check_stop_criteria <- function(pre_beta, beta) {
    result <- mapply(function(a, b) {all(sign(a) == sign(b))}, pre_beta, beta)
    all(result == TRUE)
}

#' BRAIL is Block-Randomized Adaptive Iterative Lasso
#'
#' @param X is a list containing k matrices, each of the matrices has n rows. They are the 'blocks'.
#' @param y is the response vector of length n.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, "binomial", "gaussian" and  "poisson" are available.
#' @param tau is the stability cut-off, which is used for choosing the support features in each block.
#' @param B is the number of bootstraps in each iteration.
#' @param lasso.control is a list of glmLasso related arguments, like KLB, thresh and max.iter.
#' @param ridge.control is a list of ridgeLasso related arguments, like lambda, max.iter and thresh.
#'
#' @return a list containing n vectors. Each element of the list will be a vector of length p_k, where p_k is the number of columns in the k-th block.
#'
#' @examples
#' X1 <- matrix(rnorm(500), ncol = 10)
#' X2 <- matrix(rnorm(500), ncol = 10)
#' y <- rbinom(50, 1, 0.6)
#' X <- list(X1, X2)
#' BRAIL(X, y, family = "binomial", tau = 0.8, B = 200, lasso.control= list(KLB = 10, max.iter = 1e6, thresh = 0.005), ridge.control = list(lambda = 1, max.iter = 1e4, thresh = 1e-5))

BRAIL <- function(X, y, family = "gaussian", tau = 0.8, B = 200, lasso.control = list(), ridge.control = list()) {
    beta <- lapply(X, function(x){rep(0, ncol(x))})
    K <- length(X)
    n <- length(y)
    family <- tolower(family)
    if (family == "binomial")
        family <- "logistic"

    while (TRUE) {
        pre_beta <- beta
        for (k in 1 : K) {
            pk <- ncol(X[[k]])
            beta_norm0 <- sum(abs(sign(beta[[k]])))
            c <- get_c(X[[k]], beta[[k]], family)
            tmp_lambda <- c / n * norm(as.matrix(beta[[k]]), "2") * sqrt(log(pk) * beta_norm0)
            lambda <- sapply(beta[[k]], function(x){if(x!=0) {tmp_lambda} else {2*tmp_lambda}})
            betak_samples <- 
            foreach::foreach(b = 1:B, .combine = cbind, .inorder = FALSE, .packages='MixedGraphs') %dopar% {
                indexes <- sample(nrow(X[[k]]), size = nrow(X[[k]]), replace = TRUE)
                sample_X <- lapply(X, function(x){x[indexes,]})
                sample_y <- y[indexes]
                multi_tmp <- mapply(function(x, y) {x %*% y}, sample_X, beta)
                o <- rowSums(as.matrix(multi_tmp[, -k]))
                w <- lambda * runif(pk, 0.5, 1.5)
                lasso_argv <- list(X = sample_X[[k]], y = sample_y, o = o, family = family, lambda = w)
                do.call(glmLasso, c(lasso_argv, lasso.control))$Coef[-1]
            }
            support_indexes <- which(rowSums(sign(betak_samples))/B >= tau)
            betak_support <- beta[[k]][support_indexes]
            multi_tmp <-mapply(function(x, y) {x %*% y}, X, beta)
            o <- rowSums(as.matrix(multi_tmp[, -k]))
            ridge_argv <- list(X = X[[k]][,support_indexes], y = y, o = o, family = family)
            sub_beta <- do.call(glmRidge, c(ridge_argv, ridge.control))$Coef[-1]
            beta[[k]] <- rep(0, length(beta[[k]]))
            mapply(function(index, value){beta[[k]][index] <<- value}, support_indexes, sub_beta)
        }
        if (check_stop_criteria(pre_beta, beta)) {
            break
        }
    }
    return(beta)
}