#ifndef ADMMSolver_H_
#define ADMMSolver_H_

class ADMMSolver {
private:
	arma::mat X;
	arma::vec y;
	arma::vec o;
	arma::vec z;
	arma::vec u;

public:
	ADMMSolver();
	virtual arma::vec getGradient(const arma::vec &beta, const arma::vec &z, const arma::vec &u);
	virtual arma::mat getHessian(const arma::vec &beta);
	virtual void solveBeta(arma::vec &beta, const arma::vec &z, const arma::vec &u);
	~ADMMSolver();
};

class ADMMLogistic : public ADMMSolver {
public:
	ADMMLogistic();
	~ADMMLogistic();
	
};

class ADMMPoisson : public ADMMSolver {
public:
	ADMMPoisson();
	~ADMMPoisson();
	
};

class ADMMGaussian : public ADMMSolver {
public:
	ADMMGaussian();
	~ADMMGaussian();
	
};
#endif