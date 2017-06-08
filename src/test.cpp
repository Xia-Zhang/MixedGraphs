// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "ADMM.h"
#include "NewtonSolver.h"
using namespace Rcpp;
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
//

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List test(const arma::mat& X, const arma::vec& y, const std::string method = "Gaussian", const double lambda = 0.5) {
    ADMM testADMM(X, y);
    testADMM.setWeight(lambda);
    arma::vec result = arma::vec(1);
    Rcpp::Rcout << testADMM.fit(method);
    //testADMM.setThreadNumber(4);
    //Rcpp::Rcout << testADMM.fit(method);
    return List::create(Named("Result") = result);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List testNewton(const arma::mat& X, const arma::vec& y, const std::string method = "Gaussian", const double lambda = 0.5) {
    NewtonSolver *solver;
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

    arma::vec o(X.n_rows, arma::fill::zeros);
    solver->setSolver(X, y, o, lambda);
    arma::vec result = solver->solve();
    return List::create(Named("Result") = result);
}