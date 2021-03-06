% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SIRS.R
\name{cv.sirs}
\alias{cv.sirs}
\title{CV for SIRS}
\usage{
cv.sirs(X, y, k = 10, a = seq(0, 1, 0.1))
}
\arguments{
\item{X}{the design matrix}

\item{y}{the response}

\item{k}{the number of folds}

\item{a}{where the parameter a can be chosen from}
}
\value{
the best parameter a which minimizes the loss function \code{a.best}
}
\description{
CV for SIRS to select the best parameter a
}
\examples{
\dontrun{
library(MASS)
set.seed(12345)
n <- 35
p <- 200
beta_0 <- numeric(p)
beta_0[1:7] <- c(1, -0.5, 0.7,-1.2, -0.9, 0.3, 0.55)
r <- 0.1
gamma_r <- matrix(r, p, p) + diag(1 - r, p, p)
X <- mvrnorm(n, numeric(p), gamma_r)
y <- X\%*\%beta_0
k <- 3
a <- seq(0.2, 0.4, 0.1)
a.best <- cv.sirs(X, y, k, a)
print(a.best)
}
}
