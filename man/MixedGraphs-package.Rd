\name{MixedGraphs-package}
\alias{MixedGraphs-package}
\alias{MixedGraphs}
\docType{package}
\title{
\packageTitle{MixedGraphs}
}
\description{
\packageDescription{MixedGraphs}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{MixedGraphs}
\packageIndices{MixedGraphs}
% ~~ An overview of how to use the package, including the most important functions ~~
}
\author{
Xia Zhang

Maintainer: Xia Zhang <zhangxia9403@gmail.com>
}
% \references{
% % ~~ Literature or other references for background information ~~
% }
% ~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
% \keyword{ package }
% \seealso{
% ~~ Optional links to other man pages, e.g. ~~
% ~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
% }
\examples{
    X <- matrix(rnorm(500), ncol = 10)
    y <- rpois(50, 3)
    glm(y ~ X, family = "poisson")
    test(cbind(1, X), y, "poisson", 0.0)
}
