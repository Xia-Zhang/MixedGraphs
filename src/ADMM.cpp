#include "ADMM.h"
#include <RcppParallel.h>
using namespace RcppParallel;

ADMM::ADMM (const arma::mat &X,
            const arma::vec &y,
            const arma::vec &o,
            const arma::vec &weight,
            const arma::vec &lambdas,
            const double thresh,
            const uint64_t support_stability,
            const uint64_t maxIter,
            const arma::vec &betaWS,
            const arma::vec &zWS,
            const arma::vec &uWS,
            const bool intercept) {
    reset(X, y, o, weight, lambdas, thresh, support_stability, maxIter, betaWS, zWS, uWS);
}

void ADMM::reset(const arma::mat &X,
                 const arma::vec &y,
                 const arma::vec &o,
                 const arma::vec &weight,
                 const arma::vec &lambdas,
                 const double thresh,
                 const uint64_t support_stability,
                 const uint64_t maxIter,
                 const arma::vec &betaWS,
                 const arma::vec &zWS,
                 const arma::vec &uWS,
                 const bool intercept) {
    uint64_t n = X.n_rows, p = X.n_cols;
    this->X = X;
    this->y = y;
    setVec(this->o, o, n);
    if (lambdas.empty() && weight.empty()) {
        Rcpp::stop("Please input lambda!");
    }
    if (!lambdas.empty()) {
        this->lambdas = lambdas;
    }
    else {
        setWeight(weight);
    }
    this->thresh = thresh;
    this->support_stability = support_stability;
    this->maxIter = maxIter;
    setVec(this->betaWS, betaWS, p);
    setVec(this->zWS, zWS, p);
    setVec(this->uWS, uWS, p);
    this->intercept = intercept;
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
    weight.clear();
    lambdas.clear();
    z.clear();
    preZ.clear();
    thresh = 0.0;
    support_stability = 0;
    maxIter = 0;
}

arma::mat ADMM::fit(const std::string family) {
    uint64_t k = 1, p = X.n_cols;
    // TODO: consensus ADMM partition number
    uint64_t partitions = 1;
    std::vector<ADMMSolver *> solvers(partitions);
    std::string lowerMethod(family);
    sumUBeta = arma::mat(p, partitions, arma::fill::zeros);
    arma::mat result;

    std::transform(lowerMethod.begin(), lowerMethod.end(), lowerMethod.begin(), ::tolower);
    if (lowerMethod == "binomial") {
        initializeSolver<ADMMLogistic>(solvers);
    }
    else if (lowerMethod == "poisson") {
        initializeSolver<ADMMPoisson>(solvers);
    }
    else if (lowerMethod == "gaussian"){
        initializeSolver<ADMMGaussian>(solvers);
    }
    z = zWS;

    if (lambdas.empty()) {
        while (k <= maxIter) {
            Rcpp::checkUserInterrupt();
            updateUBeta(solvers);
            updateZ();
            if (stopCriteria()) {
                break;
            }
            k++;
        }
        result = arma::conv_to<arma::mat>::from(z);
    }
    else {
        for (double lambda_tmp : lambdas) {
            setWeight(lambda_tmp);
            while (k <= maxIter) {
                Rcpp::checkUserInterrupt();
                updateUBeta(solvers);
                updateZ();
                if (stopCriteria()) {
                    supportIter = 0;
                    break;
                }
                k++;
            }
            if (result.empty()) {
                result = arma::conv_to<arma::mat>::from(z);
            }
            else {
                result = arma::join_rows(result, arma::conv_to<arma::mat>::from(z));
            }
        }
    }
    deleteSolver(solvers);
    return result;
}

void ADMM::setWeight(const double weight) {
    this->weight = arma::vec(X.n_cols);
    this->weight.fill(weight);
    if (intercept) this->weight[0] = 0;
}

void ADMM::setWeight(const arma::vec &weight) {
    if (weight.empty()) {
        Rcpp::stop("Can't set empty lambda!");
    }
    else if (weight.n_elem == X.n_cols){
        this->weight = weight;
    }
    else if (weight.n_elem == X.n_cols - 1) {
        if (intercept) {
            this->weight = arma::vec(X.n_cols, arma::fill::zeros);
        }
        else {
            this->weight = arma::vec(X.n_cols, arma::fill::ones);
        }
        for (uint64_t i = 0; i < weight.n_elem; i++) {
            this->weight[i + 1] = weight[i];
        }
    }
}

void ADMM::setWarmStartPara(const arma::vec &zWS, const arma::vec &uWS, const arma::vec &weight) {
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

void ADMM::setSupportStability(uint64_t support_stability) {
    if (support_stability < 0) {
        Rcpp::stop("The support_stability should not be negative!");
    }
    this->support_stability = support_stability;
}

void ADMM::setMaxIterator(uint64_t maxIter) {
    if (maxIter < 0) {
        Rcpp::stop("The maxIter should not be less than 0.");
    }
    this->maxIter = maxIter;
}

template <typename T>
void ADMM::initializeSolver(std::vector<ADMMSolver *> &solvers, const uint64_t partitions) {
    uint64_t n = X.n_rows;
    double interval = n / static_cast<double>(partitions);
    for (uint64_t i = 0; i < partitions; i++) {
        uint64_t col0 = i * interval, col1 = (i + 1) * interval - 1;
        arma::mat subX = X.rows(col0, col1);
        solvers[i] = new T(X.rows(col0, col1), y.subvec(col0, col1), o.subvec(col0, col1), betaWS, uWS);
    }
}

void ADMM::deleteSolver(std::vector<ADMMSolver *> &solvers) {
    for (uint64_t i = 0; i < solvers.size(); i++) {
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
    uint64_t partitions = 1;
    admmParallel admmP(solvers, z, sumUBeta);
    parallelFor(0, partitions, admmP);
    return arma::sum(sumUBeta, 1) / partitions;
}

void ADMM::updateZ() {
    this->z = softThreashold(sum(sumUBeta, 1) / sumUBeta.n_cols);
}

arma::vec ADMM::softThreashold(const arma::vec &x) {
    return sign(x) % arma::max(abs(x) - weight, arma::vec(x.n_elem, arma::fill::zeros));
}

bool ADMM::stopCriteria() {
    bool result = true;
    if (!support_stability && thresh) {
        result = (std::sqrt(arma::sum(arma::square(z - preZ))) <= thresh);
    }
    else {
        bool remainSame = arma::all(sign(preZ) == sign(z));
        if (!support_stability) {
            support_stability = 5;
        }
        if (remainSame) {
            supportIter += 1;
        }
        else {
            supportIter = 0;
        }
        result = (supportIter >= support_stability);
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
Rcpp::NumericMatrix glmLassoCPP (const arma::mat &X, 
                                 const arma::vec &y, 
                                 const arma::vec &o, 
                                 const arma::vec &weight, 
                                 const arma::vec &lambdas,
                                 const std::string family, 
                                 const uint64_t support_stability, 
                                 const double thresh, 
                                 const uint64_t maxIter, 
                                 const arma::vec& betaWS,
                                 const arma::vec &zWS,
                                 const arma::vec &uWS,
                                 const bool intercept) {
    ADMM admm(X, y, o, weight, lambdas, thresh, support_stability, maxIter, betaWS, zWS, uWS, intercept);
    return Rcpp::wrap(admm.fit(family));
}