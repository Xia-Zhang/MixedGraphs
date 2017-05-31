#include "ADMMsolver.h"
#include <cmath>

ADMMsolver::ADMMsolver(	const arma::mat &X, 
						const arma::vec &y, 
						const arma::vec &o = arma::vec() ) {
	this->X = X;
	this->y = y;
	if (o.empty()) {
		this->o = arma::vec(X.n_rows, arma::fill::zeros);
	}
}

arma::vec ADMMLogist::getGradient(const arma::vec &beta, const arma::vec &z, const arma::vec &u) {
	int n = X.n_rows, p = X.n_cols;
	double denominator;
	arma::colvec gradient(p, arma::fill::zeros);

	for (int i = 0; i < n; i++) {
		denominator = 1 + exp( - (o[i] + X.row(i).t * beta));
		gradient += (-1) * X.row(i) / denominator + y(i) * X.row(i);
	}
	gradient = gradient / n  + beta + z + u;
	return gradient;
}

arma::mat ADMMLogist::getHessian(const arma::vec &beta) {
	int n = X.n_rows, p = X.n_cols;
	double factor, denominator;
	arma::mat hessian(p, p, arma::fill::zeros);

	for (int i = 0; i < n; i++) {
		factor = exp( - (o[i] + X.row(i).t * beta));
		denominator = pow(1 + factor, 2);
		for (int j = 0; j < p; j++) {
			for (int k = 0; j < p; k++) {
				hessian(i, j) += X(i, j) * X(i, k) * factor / denominator;
			}
		}
	}
	hessian /= n;
	hessian += arma::eye<arma::mat>(p, p);
	return hessian;
}

void ADMMLogist::solveBeta(arma::vec &beta, const arma::vec &z, const arma::vec &u) {
	beta = beta - getHessian(beta).i() * getGradient(beta, z, u);
}

arma::vec ADMMPoisson::getGradient(const arma::vec &beta, const arma::vec &z, const arma::vec &u) {
	int n = X.n_rows, p = X.n_cols;
	double factor;
	arma::colvec gradient(p, arma::fill::zeros);

	for (int i = 0; i < n; i++) {
		factor = exp(o[i] + X.row(i).t * beta);
		gradient += y(i) * X.row(i) + X.row(i) * factor;
	}
	gradient = gradient / n + beta + z + u;
	return gradient;
}

arma::mat ADMMPoisson::getHessian(const arma::vec &beta) {
	int n = X.n_rows, p = X.n_cols;
	double factor;
	arma::mat hessian(p, p, arma::fill::zeros);

	for (int i = 0; i < n; i++) {
		factor = exp(o[i] + X.row(i).t * beta);
		for (int j = 0; j < p; j++)
			for (int k = 0; k < p; k++) {
				hessian(j, k) += X(i, j) * X(i, k) * factor;
			}
	}
	hessian += arma::eye<arma::mat>(p, p);
	return hessian;
}

void ADMMPoisson::solveBeta(arma::vec &beta, const arma::vec &z, const arma::vec &u) {
	beta = beta - getHessian(beta).i() * getGradient(beta, z, u);
}

void ADMMGaussian::solveBeta(arma::vec &beta, const arma::vec &z, const arma::vec &u) {
	int n = X.n_rows, p = X.n_cols;
	beta = (X.t() * X / n + arma::eye<arma::mat>(p, p)).i() * X.t() * (y - z) / n;
	beta = beta - (X.t() * X / n + arma::eye<arma::mat>(p, p)).i() * (z + u);
}