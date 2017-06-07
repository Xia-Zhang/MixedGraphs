#include "NewtonSolver.h"

NewtonSolver::NewtonSolver(const arma::mat &X,
                           const arma::vec &y,
                           const arma::vec &o,
                           double epsilon,
                           uint32_t maxIter,
                           double sigma) {
    this->X = X;
    this->y = y;
    if (o.empty())
        this->o = arma::vec(X.n_rows, arma::fill::zeros);
    else this->o = o;
    this->epsilon = epsilon;
    this->maxIter = maxIter;
    this->sigma = sigma;
}

void NewtonSolver::setSolver(const arma::mat &X,
                           const arma::vec &y,
                           const arma::vec &o,
                           double epsilon,
                           uint32_t maxIter,
                           double sigma) {
    this->X = X;
    this->y = y;
    if (o.empty())
        this->o = arma::vec(X.n_rows, arma::fill::zeros);
    else this->o = o;
    this->epsilon = epsilon;
    this->maxIter = maxIter;
    this->sigma = sigma;
}

arma::vec NewtonSolver::solve(const arma::mat &X,
                    const arma::vec &y,
                    const arma::vec &o,
                    double epsilon,
                    uint32_t maxIter,
                    double sigma){
    setSolver(X, y, o, epsilon, maxIter, sigma);
    return solve();
}

arma::vec NewtonLogistic::getGradient(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    double denominator;
    arma::vec gradient(p, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        denominator = 1 + exp( - (o[i] + arma::as_scalar(X.row(i) * beta)));
        gradient += X.row(i).t() / denominator - y(i) * X.row(i).t();
    }
    gradient = gradient / n + epsilon * beta;
    return gradient;
}

arma::mat NewtonLogistic::getHessian(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    double factor, denominator;
    arma::mat hessian(p, p, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        factor = exp( - (o[i] + arma::as_scalar(X.row(i) * beta)));
        denominator = pow(1 + factor, 2);
        for (uint32_t j = 0; j < p; j++) {
            for (uint32_t k = 0; k < p; k++) {
                hessian(j, k) += X(i, j) * X(i, k) * factor / denominator;
            }
        }
    }
    hessian /= n;
    hessian += epsilon * arma::eye<arma::mat>(p, p);
    return hessian;
}

arma::vec NewtonLogistic::solve() {
    uint32_t k = 0;
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    arma::vec beta(X.n_cols, arma::fill::zeros), d_beta;
    while (k < maxIter) {
        d_beta = getHessian(beta).i() * getGradient(beta);
        beta = beta - d_beta;
        if (sum(arma::abs(d_beta)) < sigma) {
            break;
        }
        k++;
    }
    Rcpp::Rcout << k << std::endl;
    return beta;
}

arma::vec NewtonPoisson::getGradient(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    double factor;
    arma::vec gradient(p, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        factor = exp(o[i] + arma::as_scalar(X.row(i) * beta));
        gradient += - y(i) * X.row(i).t() + X.row(i).t() * factor;
    }
    gradient = gradient / n + epsilon * beta;
    return gradient;
}

arma::mat NewtonPoisson::getHessian(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    double factor;
    arma::mat hessian(p, p, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        factor = exp(o[i] + arma::as_scalar(X.row(i) * beta));
        for (uint32_t j = 0; j < p; j++)
            for (uint32_t k = 0; k < p; k++) {
                hessian(j, k) += X(i, j) * X(i, k) * factor;
            }
    }
    hessian += epsilon * arma::eye<arma::mat>(p, p);
    return hessian;
}

arma::vec NewtonPoisson::solve() {
    uint32_t k = 0;
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    arma::vec beta(X.n_cols, arma::fill::zeros), d_beta;
    while (k < maxIter) {
        d_beta = getHessian(beta).i() * getGradient(beta);
        beta = beta - d_beta;
        if ( std::sqrt(arma::as_scalar(d_beta.t() * d_beta)) < sigma) {
            break;
        }
        k++;
    }
    return beta;
}

arma::vec NewtonGaussian::solve() {
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    uint32_t n = X.n_rows, p = X.n_cols;
    return (X.t() * X / n + epsilon * arma::eye<arma::mat>(p, p)).i() * X.t() * (y - o) / n;
}