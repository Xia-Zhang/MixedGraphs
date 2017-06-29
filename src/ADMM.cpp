#include "ADMM.h"
#include <RcppParallel.h>
using namespace RcppParallel;

ADMM::ADMM (const arma::mat &X,
            const arma::vec &y,
            const arma::vec &o,
            const arma::vec &lambda,
            double thresh,
            uint64_t KLB,
            uint64_t maxIter,
            uint64_t threadNum,
            const arma::vec &betaWS,
            const arma::vec &zWS,
            const arma::vec &uWS) {
    reset(X, y, o, lambda, thresh, KLB, maxIter, threadNum, betaWS, zWS, uWS);
}

void ADMM::reset(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o,
                 const arma::vec &lambda,
                 double thresh,
                 uint64_t KLB,
                 uint64_t maxIter,
                 uint64_t threadNum,
                 const arma::vec &betaWS,
                 const arma::vec &zWS,
                 const arma::vec &uWS) {
    uint64_t n = X.n_rows, p = X.n_cols;
    this->X = X;
    this->y = y;
    setVec(this->o, o, n);
    setWeight(lambda);
    this->thresh = thresh;
    this->KLB = KLB;
    this->maxIter = maxIter;
    this->threadNum = threadNum;
    setVec(this->betaWS, betaWS, p);
    setVec(this->zWS, zWS, p);
    setVec(this->uWS, uWS, p);
    preZ = arma::vec(p, arma::fill::zeros);
}

void ADMM::clear() {
    X.clear();
    sumUBeta.clear();
    y.clear();
    o.clear();
    betaWS.clear();
    zWS.clear();
    uWS.clear();
    lambda.clear();
    z.clear();
    preZ.clear();
    thresh = 0.0;
    KLB = 0;
    maxIter = 0;
    threadNum = 0;
}

arma::vec ADMM::fit(const std::string method) {
    uint64_t k = 1, p = X.n_cols;
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
    this->lambda = arma::vec(X.n_cols);
    this->lambda.fill(lambda);
    this->lambda[0] = 0;
}

void ADMM::setWeight(const arma::vec &weight) {
    if (weight.empty()) {
        this->lambda = arma::vec(X.n_cols, arma::fill::ones);
        this->lambda[0] = 0;
    }
    else if (weight.n_elem == X.n_cols){
        this->lambda = weight;
    }
    else if (weight.n_elem == X.n_cols - 1) {
        this-> lambda = arma::vec(X.n_cols, arma::fill::zeros);
        for (uint64_t i = 0; i < weight.n_elem; i++) {
            this->lambda[i + 1] = weight[i];
        }
    }
}

void ADMM::setWarmStartPara(const arma::vec &zWS, const arma::vec &uWS, const arma::vec &lambda) {
    uint64_t p = X.n_cols;
    setVec(this->betaWS, betaWS, p);
    setVec(this->zWS, zWS, p);
    setVec(this->uWS, uWS, p);
}

void ADMM::setInitialBeta(const arma::vec &beta) {
    setVec(this->betaWS, betaWS, X.n_cols);
}

void ADMM::setThresh(const double thresh) {
    if (thresh < 0) {
        Rcpp::stop("The thresh should not be negative!");
    }
    this->thresh = thresh;
}

void ADMM::setKLB(uint64_t KLB) {
    if (KLB < 0) {
        Rcpp::stop("The KLB should not be negative!");
    }
    this->KLB = KLB;
}

void ADMM::setMaxIterator(uint64_t maxIter) {
    if (maxIter < 0) {
        Rcpp::stop("The maxIter should not be less than 0.");
    }
    this->maxIter = maxIter;
}

void ADMM::setThreadNumber(uint64_t number) {
    if (number <= 0) {
        Rcpp::stop("The thread number should not be less than 1.");
    }
    //RcppParallel::setThreadOptions(numThreads = 4);
    this->threadNum = number;
}

template <typename T>
void ADMM::initializeSolver(std::vector<ADMMSolver *> &solvers) {
    uint64_t n = X.n_rows;
    double interval = n / static_cast<double>(threadNum);
    for (uint64_t i = 0; i < threadNum; i++) {
        uint64_t col0 = i * interval, col1 = (i + 1) * interval - 1;
        arma::mat subX = X.rows(col0, col1);
        solvers[i] = new T(X.rows(col0, col1), y.subvec(col0, col1), o.subvec(col0, col1), betaWS, uWS);
    }
}

void ADMM::deleteSolver(std::vector<ADMMSolver *> &solvers) {
    for (uint64_t i = 0; i < threadNum; i++) {
        delete solvers[i];
    }
}

struct admmParallel: public Worker {
    const std::vector<ADMMSolver *> input;
    const arma::vec &z;
    arma::mat &sum;
    admmParallel(const std::vector<ADMMSolver *> &input, const arma::vec &z, arma::mat &sum) : input(input), z(z), sum(sum){}
    void operator()(std::size_t begin, std::size_t end) {
        for (uint64_t i = begin; i < end; i++) {
            sum.col(i) = input[i]->solve(z);
        }
    }
};

arma::vec ADMM::updateUBeta(std::vector<ADMMSolver *> &solvers) {
    admmParallel admmP(solvers, z, sumUBeta);
    parallelFor(0, threadNum, admmP);
    return arma::sum(sumUBeta, 1) / threadNum;
}

void ADMM::updateZ() {
    softThreashold(sum(sumUBeta, 1) / sumUBeta.n_cols, this->z);
}

void ADMM::softThreashold(const arma::vec &sum, arma::vec &value) {
    uint64_t p = X.n_cols;
    for (uint64_t i = 0; i < p; i++) {
        if (sum[i] > lambda[i]) {
            value[i] = sum[i] - lambda[i];
        }
        else if (sum[i] < (-1) * lambda[i]) {
            value[i] = sum[i] + lambda[i];
        }
        else {
            value[i] = 0;
        }
    }
}

bool ADMM::stopCriteria() {
    bool result = true;
    if (!KLB && thresh) {
        result = (std::sqrt(arma::sum(arma::square(z - preZ))) <= thresh);
    }
    else {
        uint64_t tmpKLB = 5;
        bool remainSame = arma::all(sign(preZ) == sign(z));
        if (KLB) {
            tmpKLB = KLB;
        }
        if (remainSame) {
            supportIter += 1;
        }
        else {
            supportIter = 0;
        }
        result = (supportIter >= tmpKLB);
    }
    preZ = z;
    return result;
}

void ADMM::setVec(arma::vec &target, const arma::vec &source, const uint64_t num) {
    if (source.empty()) {
        target = arma::vec(num, arma::fill::zeros);
    }
    else {
        target = source;
    }
}

// [[Rcpp::export]]
Rcpp::List glmLasso(const arma::mat& X, const arma::vec& y, const arma::vec& o, const arma::vec &lambda, const std::string family, const uint64_t KLB, const double thresh, const uint64_t maxIter, const uint64_t threads) {
    ADMM admm(X, y, o, lambda, thresh, KLB, maxIter, threads);
    arma::vec coef = admm.fit(family);
    return Rcpp::List::create(Rcpp::Named("Coef") = coef); 
}