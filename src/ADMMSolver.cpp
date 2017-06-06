#include "ADMMSolver.h"
#include <cmath>

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
    double denominator;
    arma::vec gradient(p, arma::fill::zeros);
    for (uint32_t i = 0; i < n; i++) {
        denominator = 1 + exp( - (o[i] + arma::as_scalar(X.row(i) * beta)));
        gradient += X.row(i).t() / denominator - y(i) * X.row(i).t();
    }
    gradient = gradient / n  + beta - z + u;
    return gradient;
}

arma::mat ADMMLogistic::getHessian() {
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
    hessian += arma::eye<arma::mat>(p, p);
    return hessian;
}

void ADMMLogistic::updateBeta(const arma::vec &z) {
    beta = beta - getHessian().i() * getGradient(z);
}

arma::vec ADMMPoisson::getGradient(const arma::vec &z) {
    uint32_t n = X.n_rows, p = X.n_cols;
    double factor;
    arma::vec gradient(p, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        factor = exp(o[i] + arma::as_scalar(X.row(i) * beta));
        gradient += - y(i) * X.row(i).t() + X.row(i).t() * factor;
    }
    gradient = gradient / n + beta - z + u;
    return gradient;
}

arma::mat ADMMPoisson::getHessian() {
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
    hessian += arma::eye<arma::mat>(p, p);
    return hessian;
}

void ADMMPoisson::updateBeta(const arma::vec &z) {
    beta = beta - getHessian().i() * getGradient(z);
}

void ADMMGaussian::updateBeta(const arma::vec &z) {
    uint32_t n = X.n_rows, p = X.n_cols;
    beta = (X.t() * X / n + arma::eye<arma::mat>(p, p)).i() * X.t() * (y - o) / n;
    beta = beta + (X.t() * X / n + arma::eye<arma::mat>(p, p)).i() * (z - u);
}