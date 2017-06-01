#ifndef ADMM_H_
#define ADMM_H_

#include <cstdint>
#include <RcppArmadillo.h>
#include "ADMMSolver.h"

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

class ADMM {
private:
    arma::mat X;
    arma::mat sumUBeta;
    arma::vec y;
    arma::vec o;
    arma::vec betaWS;
    arma::vec zWS;
    arma::vec uWS;
    arma::vec w;
    arma::vec z;
    arma::vec preSupport;
    uint32_t supportIter;
    uint32_t KLB;
    uint32_t maxIter;
    uint32_t threadNum;

    void setVec(arma::vec &target, const arma::vec &source, const uint32_t num);
    bool stopCriteria();
    void softThreashold(const arma::vec &sum, arma::vec &value);
    arma::vec updateUBeta(ADMMSolver *solver);
    void updateZ();

public:
    ADMM(const arma::mat &X, 
         const arma::vec &y, 
         const arma::vec &o = arma::vec(),
         const arma::vec &betaWS = arma::vec(),
         const arma::vec &zWS = arma::vec(),
         const arma::vec &uWS = arma::vec(),
         const uint32_t KLB = 5,
         const uint32_t maxIter = 10e5,
         const uint32_t threadNum = 1);
    void reset();
    void clear();
    arma::vec fit(const std::string method);
    ~ADMM();
};
#endif