#ifndef BRAIL_H_
#define BRAIL_H_

#include "RcppArmadillo.h"
using namespace arma;

field<vec> BRAIL(field<mat> &X, const vec &y, std::string family, double tau = 0.8, uint64_t B = 200);

#endif