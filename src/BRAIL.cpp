#include <stdlib.h>
#include <random>
#include "ADMM.h"
#include "BRAIL.h"
#include "NewtonSolver.h"

bool checkStopCriteria(const field<vec> &pre_beta, const field<vec> &beta) {
    for (uint64_t i = 0; i < beta.size(); i++) {
        if (any(sign(pre_beta[i]) != sign(beta[i]))) {
            return false;
        }
    }
    return true;
}

field<vec> initializaBeta(const field<mat> &X) {
    field<vec> beta(X.n_elem);
    for (uint64_t i = 0; i < X.n_elem; i++) {
        beta[i] = vec(X[i].n_cols, fill::zeros);
    }
    return beta;
}

double getC(const mat &X, const vec &beta, std::string &method) {
    std::string lowerMethod(method);
    uint64_t n = X.n_rows;
    mat W(n, n);

    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "gaussian") {
        return 1.0;
    }
    if (lowerMethod == "logistic") {
        double doubleP;
        for (uint64_t i = 0; i < n; i++) {
            doubleP = 1 / (1 + exp(-(as_scalar(X.row(i) * beta))));
            W(i, i) = doubleP * (1 - doubleP) / n;
        }
    }
    else if (lowerMethod == "poisson") {
        for (uint64_t i = 0; i < n; i++) {
            W(i, i) = exp(as_scalar(X.row(i) * beta)) / n;
        }
    }
    return max(eig_sym(X.t() * W * X)) / max(eig_sym(X.t() * X));
}

// [[Rcpp::plugins(cpp11)]]
vec uniformRandom(uint64_t p) {
    std::random_device rd;
    std::uniform_real_distribution<double> uniform(0.5, 1.5);
    vec gamma(p);
    for (uint64_t i = 0; i < p; i++) {
        gamma[i] = uniform(rd);
    }
    return gamma;
}

void bootstrapSample(const field<mat> &X, const vec &y, field<mat> &sample_X, vec &sample_y) {
    uint64_t n = y.size(), K = X.n_elem, rand_num;
    sample_X = X;
    sample_y.set_size(n);
    std::random_device rd;
    for (uint64_t i = 0; i < n; i++) {
        rand_num = rd() % n;
        for (uint64_t k = 0; k < K; k++) {
            sample_X[k].row(i) = X[k].row(rand_num);
        }
        sample_y[i] = y[rand_num];
    }
}

// [[Rcpp::plugins(cpp11)]]
field<vec> BRAIL(field<mat> &X, const vec &y, std::string family, double tau, uint64_t B) {
    field<vec> beta = initializaBeta(X), pre_beta;
    uint64_t n = y.n_elem, K = X.n_elem, p, c, tmp_lambda;
    mat betaK_mat;
    vec lambda;

    while (true) {
        pre_beta = beta;
        for (uint64_t k = 0; k < K; k++) {
            p = X[k].n_cols;
            lambda.set_size(p);
            betaK_mat.set_size(B, p);
            c = getC(X[k], beta[k], family);
            tmp_lambda = c / n * (std::sqrt(sum(square(beta[k]))) );
            tmp_lambda *=  sqrt(log(p) * sum(abs(sign(beta[k]))) );
            lambda.fill(tmp_lambda);
            for (uint64_t j = 0; j < p; j++) {
                if (beta[k][j] == 0) {
                    lambda *= 2;
                }
            }
            //#pragma omp parallel for
            for (uint64_t b = 0; b < B; b++) {
                field<mat> sample_X(X);
                vec sample_y(n), o(n, fill::zeros), w;
                bootstrapSample(X, y, sample_X, sample_y);
                for (uint64_t l = 0; l < K; l++) {
                    if (l == k) {
                        continue;
                    }
                    o += sample_X[l] * beta[l];
                }
                w = uniformRandom(p) % lambda;
                mat cbind_tmp = mat(n, 1, fill::ones);
                ADMM admm(join_rows(cbind_tmp, sample_X[k]), sample_y, o, w);
                rowvec fit_result = admm.fit(family).t();
                betaK_mat.row(b) = fit_result.subvec(1, p);
            }
            uvec supportIndex = find(sum(sign(betaK_mat))/ static_cast<double> (B) >= tau);
            mat sub_X = X[k].cols(supportIndex);
            vec sub_o(n, fill::zeros);
            for (uint64_t l = 0; l < K; l++) {
                if (l == k) {
                    continue;
                }
                sub_o  += X[l] * beta[l];
            }
            mat cbind_tmp = mat(n, 1, fill::ones);
            NewtonSolver newton(join_rows(cbind_tmp, sub_X), y, sub_o);
            vec sub_beta = supportIndex.size() >= 1 ? newton.fit(family).subvec(1, supportIndex.size()) : vec();
            beta[k].zeros();
            uint64_t beta_i = 0;

            for (uint64_t index : supportIndex) {
                beta[k][index] = sub_beta[beta_i++];
            }
        }
        if (checkStopCriteria(pre_beta, beta)) {
            break;
        }
    }
    return beta;
}