// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"
#include "ADMM.h"
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
List test(const arma::mat& X, const arma::colvec& y) {
    Rcout << "Hello World!";

    ADMM testADMM(X, y);
    arma::vec result = testADMM.fit("logistic");
    for (unsigned int i = 0; i < result.n_elem; i++) {
        Rcout << result[i] << " ";
    }
    return List::create(Named("Result") = result);
}