
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

#' BRAIL is Block-RAndomized Adaptive Iterative Lasso
#'
#' @param X is a list containing n matrices, each of the matrices has n rows. They are the ‘blocks’.
#' @param y is the response vector of length n.
#' @param family is a description of the error distribution and link function to be used in the model. In our package, "binomial", "gaussian" and  "poisson" are available.
#' @param tau is the stability cut-off, which is used for choosing the support features in each block.
#' @param B is the number of bootstraps in each iteration.
#' @param cores is the number of cores used to parallellize the bootstrap computation. The default value is the detected core numbers using “detectCores()” in parallel package.
#' @param threads is the number of threads used for paralleling the ADMM process, the default value would be the available cores on a machine. 
#' @param lasso.KLB is the maximum iteration number which the support set unchanged in glmLasso
#' @param lasso.thresh is the precision when the solver to stop optimizing in glmLasso
#' @param lasso.max.iter is the maximum number of iterations to be performed for the optimization in glmLasso
#' @param ridge.lambda is a single value of ridge penalty in glmRidge
#' @param ridge.max.iter is the maximum number of iterations to be performed for the optimization in glmRidge
#' @param ridge.thresh is the indicates when to stop the solver in glmRidge
#'
#' @return a list containing n vectors, each of the vectors is the coefficients of corresponding bloc
#'
#' @examples
#' X1 <- matrix(rnorm(500), ncol = 10)
#' X2 <- matrix(rnorm(500), ncol = 10)
#' y <- rbinom(50, 1, 0.6)
#' X <- list(X1, X2)
#' BRAIL(X, y, family = “logistic”, tau = 0.8, B = 300, cores = 4, lasso.KLB = 10, lasso.max.iter = 1e6, lasso.thresh = 0.005, ridge.lambda = 1, ridge.max.iter = 1e4, ridge.thresh = 0.001)


BRAIL <- function(X, y, family = "gaussian", tau = 0.8, B = 200, cores = NULL, lasso.KLB = NULL, lasso.max.iter = 1e8, lasso.thresh = NULL, ridge.lambda = 0.25, ridge.max.iter = 1e8, ridge.thresh  = 1e-5) {
    beta <- lapply(X, function(x){rep(0, ncol(x))})
    K <- length(X)
    n <- length(y)
    family <- tolower(family)
    if (family == "binomial")
        family <- "logistic"

    library(doParallel)
    cl <- makeCluster(1)
    registerDoParallel(cl)
    while (TRUE) {
        pre_beta <- beta
        for (k in 1 : K) {
            pk <- ncol(X[[k]])
            beta_norm0 <- sum(abs(sign(beta[[k]])))
            c <- get_c(X[[k]], beta[[k]], family)
            tmp_lambda <- c / n * norm(as.matrix(beta[[k]]), "2") * sqrt(log(pk) * beta_norm0)
            lambda <- sapply(beta[[k]], function(x){if(x!=0) {tmp_lambda} else {2*tmp_lambda}})
            betak_samples <- 
            foreach(b = 1:B, .combine = cbind, .inorder = FALSE, .packages='MixedGraphs') %do% {
                indexes <- sample(nrow(X[[k]]), size = nrow(X[[k]]), replace = TRUE)
                sample_X <- lapply(X, function(x){x[indexes,]})
                sample_y <- y[indexes]
                multi_tmp <- mapply(function(x, y) {x %*% y}, sample_X, beta)
                o <- rowSums(as.matrix(multi_tmp[, -k]))
                w <- lambda * runif(pk, 0.5, 1.5)
                glmLasso(sample_X[[k]], sample_y, o = o, family = family, lambda = w, KLB = lasso.KLB, max.iter = lasso.max.iter, thresh = lasso.thresh)$Coef[-1]
            }
            support_indexes <- which(rowSums(sign(betak_samples))/B >= tau)
            betak_support <- beta[[k]][support_indexes]
            multi_tmp <-mapply(function(x, y) {x %*% y}, X, beta)
            o <- rowSums(as.matrix(multi_tmp[, -k]))
            sub_beta <- glmRidge(X[[k]][,support_indexes], y, o = o, family = family, lambda = ridge.lambda, max.iter = ridge.max.iter, thresh = ridge.thresh)$Coef[-1]
            beta[[k]] <- rep(0, length(beta[[k]]))
            mapply(function(index, value){beta[[k]][index] <<- value}, support_indexes, sub_beta)
        }
        if (check_stop_criteria(pre_beta, beta)) {
            break
        }
    }
    stopCluster(cl)
    return(beta)
}