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
    uint64_t n = X.n_rows;
    arma::vec vecP(n);
    for (uint64_t i = 0; i < n; i++) {
        vecP[i] = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
    }
    return X.t() * (vecP - y) / n + beta - z + u;
}

arma::mat ADMMLogistic::getHessianInv() {
    uint64_t n = X.n_rows, p = X.n_cols;
    double doubleP;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint64_t i = 0; i < n; i++) {
        doubleP = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
        W(i, i) = doubleP * (1 - doubleP);
    }
    W /= n;
    if (XX.empty()) XX = X * X.t();
    return arma::eye<arma::mat>(p, p) - X.t() * (W.i() + XX).i() * X;
}

void ADMMLogistic::updateBeta(const arma::vec &z) {
    beta = beta - getHessianInv() * getGradient(z);
}

arma::vec ADMMPoisson::getGradient(const arma::vec &z) {
    uint64_t n = X.n_rows;
    arma::vec vecV(n);

    for (uint64_t i = 0; i < n; i++) {
        vecV[i] = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * (vecV - y) / n + beta - z + u;
}

arma::mat ADMMPoisson::getHessianInv() {
    uint64_t n = X.n_rows, p = X.n_cols;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint64_t i = 0; i < n; i++) {
        W(i, i) = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    W /= n;
    if (XX.empty()) XX = X * X.t();
    return arma::eye<arma::mat>(p, p) - X.t() * (W.i() + XX).i() * X;
}

void ADMMPoisson::updateBeta(const arma::vec &z) {
    beta = beta - getHessianInv() * getGradient(z);
}

void ADMMGaussian::updateBeta(const arma::vec &z) {
    uint64_t n = X.n_rows, p = X.n_cols;
    if (tmpInv.empty()) {
        tmpInv = (X.t() * X / n + arma::eye<arma::mat>(p, p)).i();
        betaRidge = tmpInv * X.t() * (y - o) / n;
    }
    beta = betaRidge + tmpInv * (z - u);
}