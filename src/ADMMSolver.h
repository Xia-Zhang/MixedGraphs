#ifndef ADMMSolver_H_
#define ADMMSolver_H_

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

class ADMMSolver {
protected:
    arma::mat X;
    arma::vec y;
    arma::vec o;
    arma::vec beta;
    arma::vec u;

public:
    ADMMSolver(){};
    ADMMSolver(const arma::mat &X, 
               const arma::vec &y, 
               const arma::vec &o,
               const arma::vec beta,
               const arma::vec u);
    arma::vec solve(const arma::vec &z);
    virtual void updateBeta(const arma::vec &z) = 0;
    void updateU(const arma::vec &z);
    virtual ~ADMMSolver(){};
};

class ADMMLogistic : public ADMMSolver {
public:
    ADMMLogistic(const arma::mat &X, 
                 const arma::vec &y, 
                 const arma::vec &o,
                 const arma::vec beta,
                 const arma::vec z):ADMMSolver(X, y, o, beta, z) {}
    arma::vec getGradient(const arma::vec &z);
    arma::mat getHessian();
    void updateBeta(const arma::vec &z);
    ~ADMMLogistic(){};
    
};

class ADMMPoisson : public ADMMSolver {
public:
    ADMMPoisson(const arma::mat &X, 
                const arma::vec &y, 
                const arma::vec &o,
                const arma::vec beta,
                const arma::vec z):ADMMSolver(X, y, o, beta, z) {}
    arma::vec getGradient(const arma::vec &z);
    arma::mat getHessian();
    arma::vec solve(const arma::vec &z);
    void updateBeta(const arma::vec &z);
    ~ADMMPoisson(){};
    
};

class ADMMGaussian : public ADMMSolver {
private:
    arma::vec betaRidge;
    arma::mat tmpInv;

public:
    ADMMGaussian(const arma::mat &X, 
                 const arma::vec &y, 
                 const arma::vec &o,
                 const arma::vec beta,
                 const arma::vec z):ADMMSolver(X, y, o, beta, z) {}
    arma::vec solve(const arma::vec &z);
    void updateBeta(const arma::vec &z);
    ~ADMMGaussian(){};
    
};

Rcpp::List glmLasso(const arma::mat& X, const arma::vec& y, const arma::vec& o, const arma::vec &lambda, const std::string family, const double thresh, const uint64_t maxIter);
#endif