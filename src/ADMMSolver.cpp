#include "ADMMSolver.h"

ADMMSolver::ADMMSolver( const arma::mat &X,
                        const arma::vec &y,
                        const arma::vec &o,
                        const arma::vec beta,
                        const arma::vec u) {
    this->X = X;
    this->y = y;
    this->o = o;
    this->beta = beta;
    this->u = u;
}

arma::vec ADMMSolver::solve(const arma::vec &z) {
    updateU(z);
    updateBeta(z);
    return beta + u;
}

void ADMMSolver::updateU(const arma::vec &z) {
    u = u + beta - z;
}

arma::vec ADMMLogistic::getGradient(const arma::vec &z) {
    uint32_t n = X.n_rows, p = X.n_cols;
    arma::vec vecP(n);
    for (uint32_t i = 0; i < n; i++) {
        vecP[i] = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
    }
    return X.t() * (vecP - y) / n + beta - z + u;
}

arma::mat ADMMLogistic::getHessian() {
    uint32_t n = X.n_rows, p = X.n_cols;
    double doubleP;
    arma::mat W(n, n, arma::fill::zeros);
    
    for (uint32_t i = 0; i < n; i++) {
        doubleP = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
        W(i, i) = doubleP * (1 - doubleP);
    }
    return X.t() * W * X / n + arma::eye<arma::mat>(p, p);
}

void ADMMLogistic::updateBeta(const arma::vec &z) {
    beta = beta - arma::solve(getHessian(), getGradient(z));
}

arma::vec ADMMPoisson::getGradient(const arma::vec &z) {
    uint32_t n = X.n_rows, p = X.n_cols;
    arma::vec vecV(n);

    for (uint32_t i = 0; i < n; i++) {
        vecV[i] = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * (vecV - y) / n + beta - z + u;
}

arma::mat ADMMPoisson::getHessian() {
    uint32_t n = X.n_rows, p = X.n_cols;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        W(i, i) = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * W * X / n + arma::eye<arma::mat>(p, p);
}

void ADMMPoisson::updateBeta(const arma::vec &z) {
    beta = beta - arma::solve(getHessian(), getGradient(z));
}

void ADMMGaussian::updateBeta(const arma::vec &z) {
    uint32_t n = X.n_rows, p = X.n_cols;
    beta = (X.t() * X / n + arma::eye<arma::mat>(p, p)).i() * X.t() * (y - o) / n;
    beta = beta + (X.t() * X / n + arma::eye<arma::mat>(p, p)).i() * (z - u);
}