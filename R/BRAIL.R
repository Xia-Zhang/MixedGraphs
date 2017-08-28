
get_c <- function(X, beta, family) {
    n <- nrow(X)
    W <- matrix(0, nrow = n, ncol = n)
    if (family == "gaussian") {
        return(1.0)
    }
    else if(family == "binomial") {
        for (i in 1:n) {
            tmp <- 1 /(1 + exp(-(X[i,] %*% beta)))
            W[i, i] <- tmp * (1 - tmp) / n
        }
    }
    else if(family == "poisson") {
        for (i in 1:n) {
            W[i, i] <- exp(X[i,] %*% beta) / n
        }
    }
    return(max(eigen(t(X) %*% W %*% X)$values) / max(eigen(t(X) %*% X)$values))
}

check_stop_criteria <- function(prev_beta, beta) {
    result <- mapply(function(a, b) {all(sign(a) == sign(b))}, prev_beta, beta)
    all(result)
}

#' Block-Randomized Adaptive Iterative Lasso
#'
#' @description Block-Randomized Adaptive Iterative Lasso
#'
#' @param X is a list containing k matrices, each of the matrices has n rows. They are the 'blocks'.
#' @param y is the response vector of length n.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, the character string can be "binomial", "gaussian" or "poisson".
#' @param tau is the stability cut-off, which is used for choosing the support features in each block.
#' @param B is the number of bootstraps in each iteration.
#' @param doPar is logical value to indicate if the process should be run parallelly by foreach, if FALSE, the doPar will be disabled. The default value is TRUE.
#' @param lasso.control is a list of glmLasso related arguments, like support_stability, thresh and max.iter.
#' @param ridge.control is a list of ridgeLasso related arguments, like lambda, max.iter and thresh.
#'
#' @return a list containing n vectors. Each element of the list will be a vector of length p_k, where p_k is the number of columns in the k-th block.
#'
#' @examples
#' n <- 5
#' p <- 10
#' X <- lapply(1:2, function(x) {matrix(rnorm(n * p), n, p)})
#' y <- rnorm(n)
#' BRAIL(X, y, family = "gaussian", tau = 0.8, B = 20, doPar = TRUE, 
#' lasso.control= list(support_stability = 10, max.iter = 1e6),
#' ridge.control = list(max.iter = 1e4, thresh = 1e-5))
#'
#' @importFrom foreach %dopar% %do%
#' @export


BRAIL <- function(X, y, family = c("gaussian", "binomial", "poisson"), tau = 0.8, B = 200, doPar = TRUE, lasso.control = list(), ridge.control = list()) {
    K <- length(X)
    n <- length(y)
    beta <- lapply(X, function(x){
        non_zero <- min(as.integer(0.2 * n), ncol(x))
        betak <- numeric(ncol(x))
        if (non_zero > 0) betak[1 : non_zero] <- 1
        betak
    })
    scores <- vector("list", length = length(X))
    family <- tolower(family)
    family <- match.arg(family)
    if (B <= 0) stop("Invailid B!")

    while (TRUE) {
        prev_beta <- beta
        for (k in 1 : K) {
            pk <- ncol(X[[k]])
            beta_norm0 <- sum(abs(sign(beta[[k]])))
            c <- get_c(Reduce(cbind, X), Reduce(c, beta), family)
            eta <- c / n * norm(as.matrix(beta[[k]]), "2") * sqrt(log(pk) * beta_norm0)
            lambda <- ifelse(beta[[k]], eta, 2*eta)

            # Estimate support indexes, using foreach to parallelize the process
            `%myfun%` <- ifelse(doPar, `%dopar%`, `%do%`)
            betak_samples  <- 
                foreach::foreach(1:B, .combine = cbind, .inorder = FALSE, .packages='MixedGraphs') %myfun% {
                    indexes <- sample(n, size = n, replace = TRUE)
                    sample_X <- lapply(X, function(x){x[indexes,]})
                    sample_y <- y[indexes]
                    multi_tmp <- mapply(function(x, y) {x %*% y}, sample_X, beta)
                    o <- rowSums(as.matrix(multi_tmp[, -k]))
                    w <- lambda * runif(pk, 0.5, 1.5)
                    lasso_argv <- list(X = sample_X[[k]], y = sample_y, o = o, family = family, weight = w, 
                                       init.beta = prev_beta[[k]], init.z = prev_beta[[k]], init.u = sign(prev_beta[[k]]) * lambda)
                    as.vector(do.call(glmLasso_impl, c(lasso_argv, lasso.control)))[-1]
                }
            scores[[k]] <- rowSums(abs(sign(betak_samples)))/B
            support_indexes <- which(scores[[k]] >= tau)
            if (length(support_indexes) == 0) {
                support_indexes <- sample(pk, max(min(as.integer(0.2 * n), as.integer(0.5 *pk)), 1))
            }

            # Calculate the support set coefficients using newton solver
            betak_support <- beta[[k]][support_indexes]
            multi_tmp <-mapply(function(x, y) {x %*% y}, X, beta)
            o <- rowSums(as.matrix(multi_tmp[, -k]))
            ridge_argv <- list(X = X[[k]][,support_indexes], y = y, o = o, beta.init = prev_beta[[k]][support_indexes], family = family, lambda = 1e-4)
            sub_beta <- as.vector(do.call(glmRidge_impl, c(ridge_argv, ridge.control)))[-1]
            beta[[k]] <- numeric(pk)
            beta[[k]][support_indexes] <- sub_beta
        }
        # check if the support set stay unchanged
        if (check_stop_criteria(prev_beta, beta)) {
            break
        }
    }
    return(list("coefficients" = beta, "scores" = scores))
}