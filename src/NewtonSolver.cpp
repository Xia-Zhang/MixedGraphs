#include "NewtonSolver.h"

NewtonSolver::NewtonSolver(const arma::mat &X,
                           const arma::vec &y,
                           const arma::vec &o,
                           double epsilon,
                           uint32_t maxIter,
                           double sigma) {
    setSolver(X, y, o, epsilon, maxIter, sigma);
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

arma::vec NewtonSolver::fit(std::string method) {
    NewtonSolver * solver;
    std::string lowerMethod(method);

    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "logistic") {
        solver = new NewtonLogistic();
    }
    else if (lowerMethod == "poisson") {
        solver = new NewtonPoisson();
    }
    else if (lowerMethod == "gaussian"){
        solver = new NewtonGaussian();
    }

    solver->setSolver(this->X, this->y, this->o, this->epsilon);
    return solver->solve();
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

void NewtonSolver::setEpsilon(const double epsilon) {
    this->epsilon = epsilon;
}

arma::vec NewtonLogistic::getGradient(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    arma::vec vecP(n);
    for (uint32_t i = 0; i < n; i++) {
        vecP[i] = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
    }
    return X.t() * (vecP - y) / n + epsilon * beta;
}

arma::mat NewtonLogistic::getHessian(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    double doubleP;
    arma::mat W(n, n, arma::fill::zeros);
    
    for (uint32_t i = 0; i < n; i++) {
        doubleP = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
        W(i, i) = doubleP * (1 - doubleP);
    }
    return X.t() * W * X / n + epsilon * arma::eye<arma::mat>(p, p);
}

arma::vec NewtonLogistic::solve() {
    uint32_t k = 0;
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    arma::vec beta(X.n_cols, arma::fill::zeros), d_beta;
    while (k < maxIter) {
        d_beta = arma::solve(getHessian(beta), getGradient(beta));
        beta = beta - d_beta;
        if (sum(arma::sqrt(d_beta)) < sigma) {
            break;
        }
        k++;
    }
    return beta;
}

arma::vec NewtonPoisson::getGradient(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    arma::vec vecV(n);

    for (uint32_t i = 0; i < n; i++) {
        vecV[i] = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * (vecV - y) / n + epsilon * beta;
}

arma::mat NewtonPoisson::getHessian(const arma::vec &beta) {
    uint32_t n = X.n_rows, p = X.n_cols;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint32_t i = 0; i < n; i++) {
        W(i, i) = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * W * X / n + epsilon * arma::eye<arma::mat>(p, p);
}

arma::vec NewtonPoisson::solve() {
    uint32_t k = 0;
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    arma::vec beta(X.n_cols, arma::fill::zeros), d_beta;
    while (k < maxIter) {
        d_beta = arma::solve(getHessian(beta), getGradient(beta));
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