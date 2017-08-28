#ifndef NewtonSolver_H_
#define NewtonSolver_H_

#include <RcppArmadillo.h>

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

class NewtonSolver{
protected:
    arma::mat X;
    arma::vec y;
    arma::vec o;
    arma::vec betaWS;
    arma::vec lambda;
    double epsilon;
    uint64_t maxIter;
    double thresh;
    bool intercept;
public:
    NewtonSolver(){};
    NewtonSolver(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o = arma::vec(),
                 const arma::vec betaWS = arma::vec(),
                 const arma::vec lambda = arma::vec(),
                 const uint64_t maxIter = 1e8,
                 const double thresh = 1e-8,
                 const bool intercept = true);
    void setSolver(const arma::mat &X,
                   const arma::vec &y,
                   const arma::vec &o = arma::vec(),
                   const arma::vec betaWS = arma::vec(),
                   const arma::vec lambda = arma::vec(),
                   const uint64_t maxIter = 1e8,
                   const double thresh = 1e-8,
                   const bool intercept = true);
    arma::mat fit(std::string family);
    virtual arma::vec solve() {return arma::vec();};
    arma::vec solve(const arma::mat &X,
                    const arma::vec &y,
                    const arma::vec &o = arma::vec(),
                    const arma::vec betaWS = arma::vec(),
                    const arma::vec lambda = arma::vec(),
                    const uint64_t maxIter = 1e8,
                    const double thresh = 1e-8,
                    const bool intercept = true);
    void setLambda(const arma::vec lambda);
    virtual ~NewtonSolver(){};
};

class NewtonLogistic : public NewtonSolver{
private:
    arma::mat XX;
    arma::vec getBetaUpdate(const arma::vec &beta);
public:
    arma::vec solve();
    NewtonLogistic(const NewtonSolver &solver) : NewtonSolver(solver){};
    ~NewtonLogistic(){};
};

class NewtonPoisson : public NewtonSolver{
private:
    arma::mat XX;
    arma::vec getBetaUpdate(const arma::vec &beta);
public:
    arma::vec solve();
    NewtonPoisson(const NewtonSolver &solver) : NewtonSolver(solver){};
    ~NewtonPoisson(){};
};

class NewtonGaussian : public NewtonSolver{
private:
    arma::mat XX;
    arma::vec commonVec;
    arma::vec commonVecX;
public:
    arma::vec solve();
    NewtonGaussian(const NewtonSolver &solver) : NewtonSolver(solver){};
    ~NewtonGaussian(){};
};

Rcpp::NumericMatrix glmRidgeCPP(const arma::mat& X, 
                                const arma::vec& y, 
                                const arma::vec& o, 
                                const arma::vec betaWS, 
                                const arma::vec lambda,
                                const std::string family, 
                                const double thresh, 
                                const uint64_t maxIter);
#endif