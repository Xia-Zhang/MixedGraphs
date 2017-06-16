#include "ADMM.h"

ADMM::ADMM (const arma::mat &X,
            const arma::vec &y,
            const arma::vec &o,
            const arma::vec &w,
            const arma::vec &betaWS,
            const arma::vec &zWS,
            const arma::vec &uWS,
            const uint32_t KLB,
            const uint32_t maxIter,
            const uint32_t threadNum) {
    reset(X, y, o, w, betaWS, zWS, uWS, KLB, maxIter, threadNum);
}

void ADMM::reset(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o,
                 const arma::vec &w,
                 const arma::vec &betaWS,
                 const arma::vec &zWS,
                 const arma::vec &uWS,
                 const uint32_t KLB,
                 const uint32_t maxIter,
                 const uint32_t threadNum) {
    uint32_t n = X.n_rows, p = X.n_cols;
    this->X = X;
    this->y = y;
    setVec(this->o, o, n);
    setWeight(w);
    setVec(this->betaWS, betaWS, p);
    setVec(this->zWS, zWS, p);
    setVec(this->uWS, uWS, p);
    this->KLB = KLB;
    this->maxIter = maxIter;
    this->threadNum = threadNum;
    supportBeta = arma::vec(p, arma::fill::zeros);
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
    supportBeta.clear();
    KLB = 0;
    maxIter = 0;
    threadNum = 0;
}

arma::vec ADMM::fit(const std::string method) {
    uint32_t k = 1, p = X.n_cols, n = X.n_rows;
    std::vector<ADMMSolver *> solvers(threadNum);
    std::string lowerMethod(method);
    sumUBeta = arma::mat(p, threadNum, arma::fill::zeros);
    
    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "logistic") {
        initializeSolver<ADMMLogistic>(solvers);
    }
    else if (lowerMethod == "poisson") {
        initializeSolver<ADMMPoisson>(solvers);
    }
    else if (lowerMethod == "gaussian"){
        initializeSolver<ADMMGaussian>(solvers);
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
    deleteSolver(solvers);
    return z;
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

void ADMM::setWarmStartPara(const arma::vec &zWS, const arma::vec &uWS, const arma::vec &w) {
    uint32_t p = X.n_cols;
    setVec(this->betaWS, betaWS, p);
    setVec(this->zWS, zWS, p);
    setVec(this->uWS, uWS, p);
}

void ADMM::setInitialBeta(const arma::vec &beta) {
    setVec(this->betaWS, betaWS, X.n_cols);
}

void ADMM::setKLB(const int KLB) {
    if (KLB < 0) {
        Rcpp::stop("The KLB should not be less than 0.");
    }
    this->KLB = KLB;
}

void ADMM::setMaxIterator(const int maxIter) {
    if (maxIter < 0) {
        Rcpp::stop("The maxIter should not be less than 0.");
    }
    this->maxIter = maxIter;
}

void ADMM::setThreadNumber(const int number) {
    if (number <= 0) {
        Rcpp::stop("The thread number should not be less than 1.");
    }
    omp_set_num_threads(number);
    this->threadNum = number;
}


template <typename T>
void ADMM::initializeSolver(std::vector<ADMMSolver *> &solvers) {
    uint32_t n = X.n_rows;
    double interval = n / static_cast<double>(threadNum);
    for (uint32_t i = 0; i < threadNum; i++) {
        uint32_t col0 = i * interval, col1 = (i + 1) * interval - 1;
        arma::mat subX = X.rows(col0, col1);
        solvers[i] = new T(X.rows(col0, col1), y.subvec(col0, col1), o, betaWS, uWS);
    }
}

void ADMM::deleteSolver(std::vector<ADMMSolver *> &solvers) {
    for (uint32_t i = 0; i < threadNum; i++) {
        delete solvers[i];
    }
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
    bool remainSame = arma::all(sign(supportBeta) == sign(z));
    if (remainSame) {
        supportIter += 1;
    }
    else {
        supportIter = 0;
        supportBeta = z;
    }
    return supportIter >= KLB;
}

void ADMM::setVec(arma::vec &target, const arma::vec &source, const uint32_t num) {
    if (source.empty()) {
        target = arma::vec(num, arma::fill::zeros);
    }
    else {
        target = source;
    }
}