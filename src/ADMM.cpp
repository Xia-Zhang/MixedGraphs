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
}

template <typename T>
void ADMM::initialSolver(std::vector<ADMMSolver *> &solvers) {
    uint32_t n = X.n_rows;
    double interval = n / static_cast<double>(threadNum);
    for (uint32_t i = 0; i < threadNum; i++) {
        uint32_t col0 = i * interval, col1 = (i + 1) * interval - 1;
        arma::mat subX = X.rows(col0, col1);
        solvers[i] = new T(X.rows(col0, col1), y.subvec(col0, col1), o, betaWS, uWS);
    }
}

arma::vec ADMM::fit(const std::string method) {
    uint32_t k = 1, p = X.n_cols, n = X.n_rows;
    std::vector<ADMMSolver *> solvers(threadNum);
    std::string lowerMethod(method);
    sumUBeta = arma::mat(p, threadNum, arma::fill::zeros);

    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "logistic") {
        initialSolver<ADMMLogistic>(solvers);
    }
    else if (lowerMethod == "poisson") {
        initialSolver<ADMMPoisson>(solvers);
    }
    else if (lowerMethod == "gaussian"){
        initialSolver<ADMMGaussian>(solvers);
    }

    z = zWS;

    while (k <= maxIter) {
        updateUBeta(solvers);
        updateZ();
        if (stopCriteria()) {
            break;
        }
        k++;
    }
    Rcpp::Rcout << "The total iterations: " << k << std::endl;
    return z;
}

arma::vec ADMM::updateUBeta(std::vector<ADMMSolver *> &solvers) {
    uint32_t p = X.n_cols;
    #pragma omp parallel for
    for (uint32_t i = 0; i < threadNum; i++) {
         sumUBeta.col(i) = solvers[i]->solve(z);
    }
    return arma::sum(sumUBeta, 1) / threadNum;
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
    o.clear();
    betaWS.clear();
    zWS.clear();
    uWS.clear();
    w.clear();
    z.clear();
    preSupport.clear();
    KLB = 1e3;
    maxIter = 1e5;
    threadNum = 1;
}

void ADMM::reset(const arma::mat &X, 
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

void ADMM::setThreadNumber(const int number) {
    if (number <= 0) {
        Rcpp::stop("The thread number should not be less than 1.");
    }
    omp_set_num_threads(number);
    this->threadNum = number;
}

void ADMM::setMaxIterator(const int maxIter) {
    this->maxIter = maxIter;
}

void ADMM::setKLB(const int KLB) {
    this->KLB = KLB;
}