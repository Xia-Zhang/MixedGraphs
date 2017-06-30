#ifndef NewtonSolver_H_
#define NewtonSolver_H_

#include <RcppArmadillo.h>

class NewtonSolver{
protected:
    arma::mat X;
    arma::vec y;
    arma::vec o;
    double lambda;
    uint64_t maxIter;
    double thresh;  // The convergence threshhold
public:
    NewtonSolver(){};
    NewtonSolver(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o = arma::vec(),
                 double lambda = 0.25,
                 uint64_t maxIter = 1e8,
                 double thresh = 1e-8);
    void setSolver(const arma::mat &X,
                   const arma::vec &y,
                   const arma::vec &o = arma::vec(),
                   double lambda = 0.25,
                   uint64_t maxIter = 1e8,
                   double thresh = 1e-8);
    arma::vec fit(std::string method);
    virtual arma::vec solve() {return arma::vec();};
    arma::vec solve(const arma::mat &X,
                    const arma::vec &y,
                    const arma::vec &o = arma::vec(),
                    double lambda = 0.25,
                    uint64_t maxIter = 1e8,
                    double thresh = 1e-8);
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

Rcpp::List glmRidgeCPP(const arma::mat& X, const arma::vec& y, const arma::vec& o, const double &lambda, const std::string family, const double thresh, const uint64_t maxIter);
#endif