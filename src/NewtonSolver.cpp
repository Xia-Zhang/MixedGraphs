#include "NewtonSolver.h"

NewtonSolver::NewtonSolver(const arma::mat &X,
                           const arma::vec &y,
                           const arma::vec &o,
                           const arma::vec betaWS,
                           const double lambda,
                           const uint64_t maxIter,
                           const double thresh) {
    setSolver(X, y, o, betaWS, lambda, maxIter, thresh);
}

void NewtonSolver::setSolver(const arma::mat &X,
                             const arma::vec &y,
                             const arma::vec &o,
                             const arma::vec betaWS,
                             const double lambda,
                             const uint64_t maxIter,
                             const double thresh) {
    this->X = X;
    this->y = y;
    if (o.empty())
        this->o = arma::vec(X.n_rows, arma::fill::zeros);
    else this->o = o;
    if (betaWS.empty())
        this->betaWS = arma::vec(X.n_cols, arma::fill::zeros);
    else this->betaWS = betaWS;
    this->lambda = lambda;
    this->maxIter = maxIter;
    this->thresh = thresh;
}

arma::vec NewtonSolver::fit(std::string family) {
    NewtonSolver * solver = NULL;
    std::string lowerMethod(family);

    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "binomial") {
        solver = new NewtonLogistic();
    }
    else if (lowerMethod == "poisson") {
        solver = new NewtonPoisson();
    }
    else if (lowerMethod == "gaussian"){
        solver = new NewtonGaussian();
    }

    solver->setSolver(this->X, this->y, this->o, this->betaWS, this->lambda);
    arma::vec result = solver->solve();

    delete solver;
    return result;
}

arma::vec NewtonSolver::solve(const arma::mat &X,
                              const arma::vec &y,
                              const arma::vec &o,
                              const arma::vec betaWS,
                              const double lambda,
                              const uint64_t maxIter,
                              const double thresh){
    setSolver(X, y, o, betaWS, lambda, maxIter, thresh);
    return solve();
}

void NewtonSolver::setLambda(const double lambda) {
    this->lambda = lambda;
}

arma::vec NewtonLogistic::getGradient(const arma::vec &beta) {
    uint64_t n = X.n_rows;
    arma::vec vecP(n);
    for (uint64_t i = 0; i < n; i++) {
        vecP[i] = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
    }
    return X.t() * (vecP - y) / n + lambda * beta;
}

arma::mat NewtonLogistic::getHessian(const arma::vec &beta) {
    uint64_t n = X.n_rows, p = X.n_cols;
    double doubleP;
    arma::mat W(n, n, arma::fill::zeros);
    
    for (uint64_t i = 0; i < n; i++) {
        doubleP = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
        W(i, i) = doubleP * (1 - doubleP);
    }
    return X.t() * W * X / n + lambda * arma::eye<arma::mat>(p, p);
}

arma::vec NewtonLogistic::solve() {
    uint64_t k = 0;
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    arma::vec beta, d_beta;
    beta = betaWS;
    while (k < maxIter) {
        d_beta = arma::solve(getHessian(beta), getGradient(beta));
        beta = beta - d_beta;
        if (std::sqrt(sum(arma::square(d_beta))) < thresh) {
            break;
        }
        k++;
    }
    return beta;
}

arma::vec NewtonPoisson::getGradient(const arma::vec &beta) {
    uint64_t n = X.n_rows;
    arma::vec vecV(n);

    for (uint64_t i = 0; i < n; i++) {
        vecV[i] = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * (vecV - y) / n + lambda * beta;
}

arma::mat NewtonPoisson::getHessian(const arma::vec &beta) {
    uint64_t n = X.n_rows, p = X.n_cols;
    arma::mat W(n, n, arma::fill::zeros);

    for (uint64_t i = 0; i < n; i++) {
        W(i, i) = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    return X.t() * W * X / n + lambda * arma::eye<arma::mat>(p, p);
}

arma::vec NewtonPoisson::solve() {
    uint64_t k = 0;
    if (X.empty() || y.empty()) {
        Rcpp::stop("The input matrix and response vector shouldn't be empty!");
    }
    arma::vec beta, d_beta;
    beta = betaWS;
    while (k < maxIter) {
        d_beta = arma::solve(getHessian(beta), getGradient(beta));
        beta = beta - d_beta;
        if ( std::sqrt(arma::as_scalar(d_beta.t() * d_beta)) < thresh) {
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
    uint64_t n = X.n_rows, p = X.n_cols;
    return (X.t() * X / n + lambda * arma::eye<arma::mat>(p, p)).i() * X.t() * (y - o) / n;
}

// [[Rcpp::export]]
Rcpp::List glmRidgeCPP(const arma::mat& X, const arma::vec& y, const arma::vec& o, const arma::vec& betaWS, const double &lambda, const std::string family, const double thresh, const uint64_t maxIter) {
    NewtonSolver newton(X, y, o, betaWS, lambda, maxIter, thresh);
    arma::vec coef = newton.fit(family);
    return Rcpp::List::create(Rcpp::Named("Coef") = coef);
}