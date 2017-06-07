#ifndef NewtonSolver_H_
#define NewtonSolver_H_

#include <RcppArmadillo.h>

class NewtonSolver{
protected:
    arma::mat X;
    arma::vec y;
    arma::vec o;
    double epsilon;
    uint32_t maxIter;
    double sigma;  // The convergence threshhold
public:
    NewtonSolver(){};
    NewtonSolver(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o = arma::vec(),
                 double epsilon = 0.25,
                 uint32_t maxIter = 1e5,
                 double sigma = 1e-5);
    void setSolver(const arma::mat &X,
              const arma::vec &y,
              const arma::vec &o = arma::vec(),
              double epsilon = 0.25,
              uint32_t maxIter = 1e5,
              double sigma = 1e-5);
    virtual arma::vec solve(){};
    arma::vec solve(const arma::mat &X,
               const arma::vec &y,
               const arma::vec &o = arma::vec(),
               double epsilon = 0.25,
               uint32_t maxIter = 1e5,
               double sigma = 1e-5);
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


#endif