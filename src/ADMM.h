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
    arma::vec preZ;
    double thresh;
    uint64_t supportIter;
    uint64_t KLB;
    uint64_t maxIter;
    uint64_t threadNum;

    template <typename T>
    void initializeSolver(std::vector<ADMMSolver *> &solvers);
    void deleteSolver(std::vector<ADMMSolver *> &solvers);
    arma::vec updateUBeta(std::vector<ADMMSolver *> &solvers);
    void updateZ();
    void softThreashold(const arma::vec &sum, arma::vec &value);
    bool stopCriteria();
    void setVec(arma::vec &target, const arma::vec &source, const uint64_t num);

public:
    ADMM(const arma::mat &X,
         const arma::vec &y,
         const arma::vec &o = arma::vec(),
         const arma::vec &w = arma::vec(),
         const arma::vec &betaWS = arma::vec(),
         const arma::vec &zWS = arma::vec(),
         const arma::vec &uWS = arma::vec(),
         double thresh = 0.0,
         uint64_t KLB = 0,
         uint64_t maxIter = 1e8,
         uint64_t threadNum = 1);
    void reset(const arma::mat &X, 
               const arma::vec &y, 
               const arma::vec &o = arma::vec(),
               const arma::vec &w = arma::vec(),
               const arma::vec &betaWS = arma::vec(),
               const arma::vec &zWS = arma::vec(),
               const arma::vec &uWS = arma::vec(),
               double thresh = 0.0,
               uint64_t KLB = 0,
               uint64_t maxIter = 1e8,
               uint64_t threadNum = 1);
    void clear();
    arma::vec fit(const std::string method);
    void setWeight(const double lambda);
    void setWeight(const arma::vec &weight = arma::vec());
    void setWarmStartPara(const arma::vec &zWS, const arma::vec &uWS, const arma::vec &w);
    void setInitialBeta(const arma::vec &beta);
    void setThresh(double thresh);
    void setKLB(uint64_t KLB);
    void setMaxIterator(uint64_t maxIter);
    void setThreadNumber(uint64_t number);
    ~ADMM(){};
};
#endif