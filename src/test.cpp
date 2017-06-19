// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "ADMM.h"
#include "NewtonSolver.h"
#include "BRAIL.h"
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
    testADMM.setKLB(8);
    testADMM.setWeight(lambda);
    arma::vec result;
    arma::uvec indices;
    indices << 2 << 3;
    result = testADMM.fit(method);
    Rcout << result.size();
    Rcout << result.elem(indices);
    //testADMM.setThreadNumber(4);
    //Rcpp::Rcout << testADMM.fit(method);
    return List::create(Named("Result") = result);
}

// [[Rcpp::export]]
List testADMM(const arma::mat& X, const arma::vec& y, const std::string method = "Gaussian", const double lambda = 0.5) {
    ADMM testADMM(X, y);
    testADMM.setWeight(lambda);
    arma::vec result = testADMM.fit(method);
    return List::create(Named("Result") = result);
}

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::export]]
List testNewton(const arma::mat& X, const arma::vec& y, const std::string method = "Gaussian", const double lambda = 0.5) {
    NewtonSolver newtonSolver(X, y);
    newtonSolver.setEpsilon(lambda);
    arma::vec result = newtonSolver.fit(method);
    return List::create(Named("Result") = result);
}

// [[Rcpp::export]]
List testBRAIL() {
    arma::mat i_mat(5, 5);
    i_mat.fill(1.1);
    List X = List::create(i_mat, i_mat);
    arma::vec y(5);
    y.fill(2.1);
    arma::field<arma::mat> test_X(X.size());
    for (uint64_t i = 0; i < X.size(); i++) {
        test_X[i] = as<arma::mat>(X[i]);
    }
    BRAIL(test_X, y, "poisson");
    return X;
}