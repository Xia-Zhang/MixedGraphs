#include "NewtonSolver.h"

NewtonSolver::NewtonSolver(const arma::mat &X,
                           const arma::vec &y,
                           const arma::vec &o,
                           const arma::vec betaWS,
                           const double lambda,
                           const uint64_t maxIter,
                           const double thresh,
                           const bool intercept) {
    setSolver(X, y, o, betaWS, lambda, maxIter, thresh, intercept);
}

void NewtonSolver::setSolver(const arma::mat &X,
                             const arma::vec &y,
                             const arma::vec &o,
                             const arma::vec betaWS,
                             const double lambda,
                             const uint64_t maxIter,
                             const double thresh,
                             const bool intercept) {
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
    this->intercept = intercept;
}

arma::vec NewtonSolver::fit(std::string family) {
    NewtonSolver * solver = NULL;
    std::string lowerMethod(family);

    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "binomial") {
        solver = new NewtonLogistic(*this);
    }
    else if (lowerMethod == "poisson") {
        solver = new NewtonPoisson(*this);
    }
    else if (lowerMethod == "gaussian"){
        solver = new NewtonGaussian(*this);
    }

    // solver->setSolver(this->X, this->y, this->o, this->betaWS, this->lambda);
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
                              const double thresh,
                              const bool intercept){
    setSolver(X, y, o, betaWS, lambda, maxIter, thresh, intercept);
    return solve();
}

void NewtonSolver::setLambda(const double lambda) {
    this->lambda = lambda;
}

arma::vec NewtonLogistic::getGradient(const arma::vec &beta) {
    uint64_t n = X.n_rows;
    arma::vec vecP(n), lambdaBeta;
    for (uint64_t i = 0; i < n; i++) {
        vecP[i] = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
    }
    lambdaBeta = lambda * beta;
    if (intercept) lambdaBeta[0] = 0;
    return X.t() * (vecP - y) / n + lambdaBeta;
}

arma::mat NewtonLogistic::getHessian(const arma::vec &beta) {
    uint64_t n = X.n_rows, p = X.n_cols;
    double doubleP;
    arma::mat W(n, n, arma::fill::zeros), lambdaI;
    
    for (uint64_t i = 0; i < n; i++) {
        doubleP = 1 / (1 + exp(-(o[i] + arma::as_scalar(X.row(i) * beta))));
        W(i, i) = doubleP * (1 - doubleP);
    }
    lambdaI = lambda * arma::eye<arma::mat>(p, p);
    if (intercept) lambdaI(0, 0) = 0;
    return X.t() * W * X / n + lambdaI;
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
    arma::vec vecV(n), lambdaBeta;

    for (uint64_t i = 0; i < n; i++) {
        vecV[i] = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    lambdaBeta = lambda * beta;
    if (intercept) lambdaBeta[0] = 0;
    return X.t() * (vecV - y) / n + lambdaBeta;
}

arma::mat NewtonPoisson::getHessian(const arma::vec &beta) {
    uint64_t n = X.n_rows, p = X.n_cols;
    arma::mat W(n, n, arma::fill::zeros), lambdaI;

    for (uint64_t i = 0; i < n; i++) {
        W(i, i) = exp(o[i] + arma::as_scalar(X.row(i) * beta));
    }
    lambdaI = lambda * arma::eye<arma::mat>(p, p);
    if (intercept) lambdaI(0, 0) = 0;
    return X.t() * W * X / n + lambdaI;
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
    arma::mat lambdaI = lambda * arma::eye<arma::mat>(p, p);
    if (intercept) {
        lambdaI(0, 0) = 0;
    }
    return (X.t() * X / n + lambdaI).i() * X.t() * (y - o) / n;
}

// [[Rcpp::export]]
Rcpp::NumericVector glmRidgeCPP(const arma::mat& X, 
                                const arma::vec& y, 
                                const arma::vec& o, 
                                const arma::vec& betaWS, 
                                const double &lambda, 
                                const std::string family, 
                                const double thresh, 
                                const uint64_t maxIter,
                                const bool intercept) {
    NewtonSolver newton(X, y, o, betaWS, lambda, maxIter, thresh, intercept);
    return Rcpp::wrap(newton.fit(family));
}