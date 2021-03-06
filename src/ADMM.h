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
    arma::vec weight;
    arma::vec lambdas;
    arma::vec z;
    arma::vec preZ;
    double thresh;
    uint64_t supportIter;
    uint64_t support_stability;
    uint64_t maxIter;
    bool intercept;

    template <typename T>
    void initializeSolver(std::vector<ADMMSolver *> &solvers, const uint64_t partitions = 1);
    void deleteSolver(std::vector<ADMMSolver *> &solvers);
    arma::vec updateUBeta(std::vector<ADMMSolver *> &solvers);
    void updateZ();
    arma::vec softThreashold(const arma::vec &x);
    bool stopCriteria();
    void setVec(arma::vec &target, const arma::vec &source, const uint64_t num);

public:
    ADMM(const arma::mat &X,
         const arma::vec &y,
         const arma::vec &o = arma::vec(),
         const arma::vec &weight = arma::vec(),
         const arma::vec &lambdas = arma::vec(),
         const double thresh = 0.0,
         const uint64_t support_stability = 0,
         const uint64_t maxIter = 1e8,
         const arma::vec &betaWS = arma::vec(),
         const arma::vec &zWS = arma::vec(),
         const arma::vec &uWS = arma::vec(),
         const bool intercept = TRUE);
    void reset(const arma::mat &X, 
               const arma::vec &y, 
               const arma::vec &o = arma::vec(),
               const arma::vec &weight = arma::vec(),
               const arma::vec &lambdas = arma::vec(),
               const double thresh = 0.0,
               const uint64_t support_stability = 0,
               const uint64_t maxIter = 1e8,
               const arma::vec &betaWS = arma::vec(),
               const arma::vec &zWS = arma::vec(),
               const arma::vec &uWS = arma::vec(),
               const bool intercept = TRUE);
    void clear();
    arma::mat fit(const std::string method);
    void setWeight(const double weight);
    void setWeight(const arma::vec &weight = arma::vec());
    void setWarmStartPara(const arma::vec &zWS, const arma::vec &uWS, const arma::vec &weight);
    void setInitialBeta(const arma::vec &beta);
    void setThresh(double thresh);
    void setSupportStability(uint64_t support_stability);
    void setMaxIterator(uint64_t maxIter);
    ~ADMM(){};
};

Rcpp::NumericMatrix glmLassoCPP (const arma::mat &X,
                                 const arma::vec &y,
                                 const arma::vec &o = arma::vec(),
                                 const arma::vec &weight = arma::vec(),
                                 const arma::vec &lambdas = arma::vec(),
                                 const std::string family = "gaussian",
                                 const uint64_t support_stability = 0,
                                 const double thresh = 0.0,
                                 const uint64_t maxIter = 1e8,
                                 const arma::vec &betaWS = arma::vec(),
                                 const arma::vec &zWS = arma::vec(),
                                 const arma::vec &uWS = arma::vec(),
                                 const bool intercept = TRUE);

#endif