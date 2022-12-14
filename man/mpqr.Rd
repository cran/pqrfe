% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PQRFE_main.R
\name{mpqr}
\alias{mpqr}
\title{Multiple penalized quantile regression}
\usage{
mpqr(x, y, subj, tau = 1:9/10, effect = "simple", c = 0)
}
\arguments{
\item{x}{Numeric matrix, covariates}

\item{y}{Numeric vector, outcome.}

\item{subj}{Numeric vector, identifies the unit to which the observation belongs.}

\item{tau}{Numeric vector, identifies the percentiles.}

\item{effect}{Factor, "simple" simple regression, "fixed" regression with fixed effects, "lasso" penalized regression with fixed effects.}

\item{c}{Numeric, 0 is quantile, Inf is expectile, any number between zero and infinite is M-quantile.}
}
\value{
Beta Numeric array, with three dimmensions: 1) tau, 2) coef., lower bound, upper bound, 3) exploratory variables.

Beta array with dimension (ntau, 3, d), where Beta[i,1,k] is the i-th tau estimation
of beta_k, Beta[i,2,k] is the i-th tau lower bound 95\% confidence of beta_k, and Beta[i,3,k] 
is the i-th tau lower bound 95\% confidence of beta_k.
}
\description{
Estimate penalized quantile regression for several taus
}
\examples{
n = 10
m = 5
d = 4
N = n*m
L = N*d
x = matrix(rnorm(L), ncol=d, nrow=N)
subj = rep(1:n, each=m)
alpha = rnorm(n)
beta = rnorm(d)
eps = rnorm(N)
y = as.vector(x \%*\% beta + rep(alpha, each=m) + eps)

Beta = mpqr(x,y,subj,tau=1:9/10, effect="fixed", c = 1.2)
Beta

}
