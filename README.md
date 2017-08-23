Graphical Models for Mixed Multi Modal Data
=====
[![Build Status](https://travis-ci.com/Xia-Zhang/MixedGraphs.svg?token=oYxg4uPnDpxizy9yT9x8&branch=master)](https://travis-ci.com/Xia-Zhang/MixedGraphs) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

The is the [project](https://summerofcode.withgoogle.com/projects/#5375151708307456) of [GSoC 2017](https://summerofcode.withgoogle.com/projects/). You can see more detail at [homepage](http://xia-zhang.github.io/MixedGraphs).

### [Background](https://github.com/rstats-gsoc/gsoc2017/wiki/Graphical-Models-for-Mixed-Multi-Modal-Data)

Graphical models provide a powerful and flexible framework for understanding complex multivariate data. These models, sometimes also referred to as network models, capture dependencies in multivariate data, allowing statisticians to discover underlying connections among measured variables. These models have been widely used in applied statistics and machine learning, with particular success in genetics, neuroscience, and finance. 

In big-data domain, the observed variables usually consist of multiple types, including binary, count-valued, continuous, categorical, bounded, etc. But classical graphical models, typically assume that the data are generated from a multivariate Gaussian and that each variable is marginally Gaussian. While mathematically tractable, this assumption is plainly inappropriate for mixed multi-modal data, so they cannot apply to multi modal data directly. 

Recently, mentors proposed a block randomized adaptive iterative lasso ("BRAIL") procedure to fit the mixed graphical models. In this project, we propose a new package to make graphical models for mixed multi-modal data readily available to a wide audience. The proposed package will allow for fitting, simulating from, and visualizing mixed graphical models. 

### Implementations
Five main works of the summer
- [ADMM](http://stanford.edu/~boyd/admm.html) framework for l1-penalized Gaussian, Logistic and Poisson regression with warm start and early stopping based on support convergence in C++
- Newton with l2-penalized Gaussian, Logistic and Poisson regression in C++
- the BRAIL algorithm with [foreach](https://cran.r-project.org/web/packages/foreach) parallelization in R
- the MixedGraph fitting routine with [foreach](https://cran.r-project.org/web/packages/foreach/) parallelization in R
- the plotting of MixedGraph object using [igraph](http://igraph.org/r/), [Cytoscape](http://www.cytoscape.org/) and [Cytoscape.js](http://js.cytoscape.org/)

**Note**: In package, we use [Rcpp](http://www.rcpp.org/) and [RcppArmadillo](https://cran.r-project.org/web/packages/RcppArmadillo/) to integrate R and C++ code. For the algorithm detail, you can get from [here](https://xia-zhang.github.io/MixedGraphs/articles/glmLasso.html).


### Installation
Install from github
```r
library(devtools)
install_github("Xia-Zhang/MixedGraphs")
```

### Usage
- glmLasso
    ```r
    X <- matrix(rnorm(10 * 200), 10, 200)
    y <- rbinom(10, 1, 0.6)
    glmLasso(X, y, lambda = 0.5, family = "binomial", support_stability = 10) 
    ```
- glmRidge
    ```r
    glmRidge(X, y, lambda = 0.5, family = "binomial", thresh = 0.005)
    ```
- BRAIL
    ```r
    X <- lapply(1:2, function(x) {matrix(rnorm(10 * 200), 10, 200)})
    y <- rnorm(10)
    BRAIL(X, y, family = "gaussian", tau = 0.8, B = 20, doPar = TRUE)
    ```

- MixedGraph
    ```r
    X <- lapply(1 : 3, function(x){matrix(rnorm(12), nrow = 4)})
    crf_structure = matrix(c(1, 0, 1, 1, 1, 1, 0, 0, 1), 3, 3)
    brail_control <- list(B = 5, tau = 0.6)
    G <- MixedGraph(X, crf_structure, brail_control = brail_control)
    ```
- plot.MixedGraph
    ```r
    plot(G, method = "igraph",  weighted = TRUE)
    plot(G, method = "cytoscape", layout = "")
    plot(G, method = "cytoscape.js", "attributes-layout")
    ```
### Student
[Xia Zhang](https://github.com/Xia-Zhang)  
Department of Computer Science and Technology, Peking University

### Mentors
[Genevera Allen](http://www.stat.rice.edu/~gallen)
Departments of Statistics and ECE, Rice University  
Jan and Dan Duncan Neurological Research Institute, Baylor College of Medicine and Texas Children’s Hospital

[Michael Weylandt](https://github.com/michaelweylandt)  
Department of Statistics, Rice University

### License
GPL (>= 2)
