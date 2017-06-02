#include "ADMM.h"

ADMM::ADMM( const arma::mat &X, 
            const arma::vec &y, 
            const arma::vec &o,
            const arma::vec &betaWS,
            const arma::vec &zWS,
            const arma::vec &uWS,
            const arma::vec &w,
            const uint32_t KLB,
            const uint32_t maxIter,
            const uint32_t threadNum) {
    this->X = X;
    this->y = y;
    uint32_t n = X.n_rows, p = X.n_cols;
    setVec(this->o, o, n);
    setVec(this->betaWS, betaWS, p);
    setVec(this->zWS, zWS, p);
    setVec(this->uWS, uWS, p);
    setWeight(w);
    this->KLB = KLB;
    this->maxIter = maxIter;
    this->threadNum = threadNum;

    sumUBeta = arma::mat(p, threadNum, arma::fill::zeros);
}


arma::vec ADMM::fit(const std::string method) {
    uint32_t k = 1;
    ADMMSolver *solver;
    std::string lowerMethod(method);
    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "logistic") {
        solver = new ADMMLogistic(X, y, o, betaWS, uWS);
    }
    else if (lowerMethod == "poisson") {
        solver = new ADMMPoisson(X, y, o, betaWS, uWS);
    }
    else if (lowerMethod == "gaussian"){
        solver = new ADMMGaussian(X, y, o, betaWS, uWS);
    }

    z = zWS;

    while (k <= maxIter) {
        sumUBeta.col(0) = updateUBeta(solver);
        updateZ();
        if (stopCriteria()) {
            break;
        }
        k++;
    }
    Rcpp::Rcout << "The total iterations: " << k << std::endl;
    return z;
}


arma::vec ADMM::updateUBeta(ADMMSolver *solver) {
    // TODO: parallelization
    return solver->solve(z);
}

void ADMM::updateZ() {
    softThreashold(sum(sumUBeta, 1) / sumUBeta.n_cols, this->z);
}

void ADMM::softThreashold(const arma::vec &sum, arma::vec &value) {
    uint32_t p = X.n_cols;
    for (uint32_t i = 0; i < p; i++) {
        if (sum[i] > w[i]) {
            value[i] = sum[i] - w[i];
        }
        else if (sum[i] < (-1) * w[i]) {
            value[i] = sum[i] + w[i];
        }
        else {
            value[i] = 0;
        }
    }
}

bool ADMM::stopCriteria() {
    if (preSupport.empty()) {
        preSupport = arma::vec(z.n_elem, arma::fill::zeros);
        supportIter = -1;
    }
    bool remainSame = true;
    for (uint32_t i = 0; i < z.n_elem; i++) {
        if (remainSame) {
            if (preSupport[i] == !z[i]) {
                remainSame = false;
                preSupport[i] = (z[i] == 0) ? 0 : 1;
            }
        }
        else {
            preSupport[i] = (z[i] == 0) ? 0 : 1;
        }
    }
    if (remainSame) {
        supportIter += 1;
    }
    else {
        supportIter = 0;
    }
    return supportIter >= KLB;
}

void ADMM::clear() {
    X.clear();
    sumUBeta.clear();
    y.clear();
}

void ADMM::reset() {

}

void ADMM::setVec(arma::vec &target, const arma::vec &source, const uint32_t num) {
    if (source.empty()) {
        target = arma::vec(num, arma::fill::zeros);
    }
    else {
        target = source;
    }
}

void ADMM::setWeight(const double lambda) {
    this->w = arma::vec(X.n_cols);
    this->w.fill(lambda);
    this->w[0] = 0;
}

void ADMM::setWeight(const arma::vec &weight) {
    if (weight.empty()) {
        this->w = arma::vec(X.n_cols, arma::fill::ones);
        this->w[0] = 0;
    }
    else if (weight.n_elem == X.n_cols){
        this->w = weight;
    }
    else if (weight.n_elem == X.n_cols - 1) {
        this-> w = arma::vec(X.n_cols, arma::fill::zeros);
        for (uint32_t i = 0; i < weight.n_elem; i++) {
            this->w[i + 1] = weight[i];
        }
    }
}