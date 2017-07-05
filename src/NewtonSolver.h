#ifndef NewtonSolver_H_
#define NewtonSolver_H_

#include <RcppArmadillo.h>

class NewtonSolver{
protected:
    arma::mat X;
    arma::vec y;
    arma::vec o;
    arma::vec betaWS;
    double lambda;
    uint64_t maxIter;
    double thresh;
public:
    NewtonSolver(){};
    NewtonSolver(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o = arma::vec(),
                 const arma::vec betaWS = arma::vec(),
                 const double lambda = 0.25,
                 const uint64_t maxIter = 1e8,
                 const double thresh = 1e-8);
    void setSolver(const arma::mat &X,
                   const arma::vec &y,
                   const arma::vec &o = arma::vec(),
                   const arma::vec betaWS = arma::vec(),
                   const double lambda = 0.25,
                   const uint64_t maxIter = 1e8,
                   const double thresh = 1e-8);
    arma::vec fit(std::string family);
    virtual arma::vec solve() {return arma::vec();};
    arma::vec solve(const arma::mat &X,
                    const arma::vec &y,
                    const arma::vec &o = arma::vec(),
                    const arma::vec betaWS = arma::vec(),
                    const double lambda = 0.25,
                    const uint64_t maxIter = 1e8,
                    const double thresh = 1e-8);
    void setLambda(const double lambda = 0.25);
    virtual ~NewtonSolver(){};
};

class NewtonLogistic : public NewtonSolver{
private:
    arma::vec getGradient(const arma::vec &beta);
    arma::mat getHessian(const arma::vec &beta);
public:
    arma::vec solve();
    NewtonLogistic(){};
    ~NewtonLogistic(){};
};

class NewtonPoisson : public NewtonSolver{
private:
    arma::vec getGradient(const arma::vec &beta);
    arma::mat getHessian(const arma::vec &beta);
public:
    arma::vec solve();
    NewtonPoisson(){};
    ~NewtonPoisson(){};
};

class NewtonGaussian : public NewtonSolver{
public:
    arma::vec solve();
    NewtonGaussian(){};
    ~NewtonGaussian(){};
};

Rcpp::List glmRidgeCPP(const arma::mat& X, const arma::vec& y, const arma::vec& o, const arma::vec& betaWS, const double &lambda, const std::string family, const double thresh, const uint64_t maxIter);
#endif