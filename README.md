Graphical Models for Mixed Multi Modal Data
=====
[![Build Status](https://travis-ci.com/Xia-Zhang/MixedGraphs.svg?token=oYxg4uPnDpxizy9yT9x8&branch=master)](https://travis-ci.com/Xia-Zhang/MixedGraphs) [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html)

An R package for the mixed graphical models. You can see more detail at [homepage](http://xia-zhang.github.io/MixedGraphs).

### Background
The is the project of [GSoC 2017](https://summerofcode.withgoogle.com/projects/).

In big-data domain, the observed variables usually consist of multiple types, including binary, count-valued, continuous, categorical, bounded, etc. But traditional graphical models, always assume the data are in the same type, so they cannot apply to multi modal data directly. 

Recently, mentors proposed a block randomized adaptive iterative lasso ("BRAIL") procedure to fit the mixed graphical models. In order to make the powerful new method available to a wide audience, we develop the "MixedGraphs" package.

### Implementations
Five main works of the summer
- [ADMM](http://stanford.edu/~boyd/admm.html) framework for l1-penalized Gaussian, Logistic and Poisson regression with warm start and early stopping based on support convergence in C++
- Newton with l2-penalized Gaussian, Logistic and Poisson regression in C++
- the BRAIL algorithm with [foreach](https://cran.r-project.org/web/packages/foreach) parallelization in R
- the MixedGraph fitting routine with [foreach](https://cran.r-project.org/web/packages/foreach/) parallelization in R
- the plotting of MixedGraph object using [igraph](http://igraph.org/r/), [Cytoscape](http://www.cytoscape.org/) and [Cytoscape.js](http://js.cytoscape.org/)

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
