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

void ADMMLogistic::updateBeta(const arma::vec &z) {
    uint64_t n = X.n_rows, p = X.n_cols;
    double doubleP;
    arma::vec vecP(n), gradient;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint64_t i = 0; i < n; i++) {
        vecP[i] = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
    }
    gradient = X.t() * (vecP - y) / n + beta - z + u;

    W.diag() = vecP % (1 - vecP);
    if (XX.empty()) XX = X * X.t();

    if (n <= 2*p) beta = beta - gradient + X.t() * (arma::diagmat(1 / W.diag()) * n + XX).i() * X * gradient;
    else beta = beta - arma::solve(arma::eye(p, p) + X.t() * W / n * X, gradient);
}

void ADMMPoisson::updateBeta(const arma::vec &z) {
    uint64_t n = X.n_rows, p = X.n_cols;
    arma::vec vecV(n), gradient;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint64_t i = 0; i < n; i++) {
        vecV[i] = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    gradient = X.t() * (vecV - y) / n + beta - z + u;

    W.diag() = vecV;
    if (XX.empty()) XX = X * X.t();

    if (n <= 2*p) beta = beta - gradient + X.t() * (arma::diagmat(1 / W.diag()) * n + XX).i() * X * gradient;
    else beta = beta - arma::solve(arma::eye(p, p) + X.t() * W / n * X, gradient);
}

void ADMMGaussian::updateBeta(const arma::vec &z) {
    uint64_t n = X.n_rows, p = X.n_cols;
    if (tmpInv.empty()) {
        if (n <= 2*p) {
            tmpInv = arma::eye<arma::mat>(p, p) - X.t() * (arma::eye<arma::mat>(n, n) * n + X * X.t()).i() * X;
        }
        else {
            tmpInv = (X.t() * X / n + arma::eye<arma::mat>(p, p)).i();
        }
        betaRidge = tmpInv * X.t() * (y - o) / n;
    }
    beta = betaRidge + tmpInv * (z - u);
}