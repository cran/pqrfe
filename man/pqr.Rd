% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PQRFE_main.R
\name{pqr}
\alias{pqr}
\title{Penalized quantile regression with fixed effects}
\usage{
pqr(x, y, subj, tau = 0.5, effect = "simple", c = 1)
}
\arguments{
\item{x}{Numeric matrix, covariates}

\item{y}{Numeric vector, outcome.}

\item{subj}{Numeric vector, identifies the unit to which the observation belongs.}

\item{tau}{Numeric scalar between zero and one, identifies the percentile.}

\item{effect}{Factor, "simple" simple regression, "fixed" regression with fixed effects, "lasso" penalized regression with fixed effects.}

\item{c}{Numeric, 0 is quantile, Inf is expectile, any number between zero and infinite is M-quantile.}
}
\value{
alpha       Numeric vector, intercepts' coefficients.

beta        Numeric vector, exploratory variables' coefficients.

lambda      Numeric, estimated lambda.

res         Numeric vector, percentile residuals.

tau         Numeric scalar, the percentile.

penalty     Numeric scalar, indicate the chosen effect.

c           Numeric scalar, indicate the chosen c.

sig2_alpha  Numeric vector, intercepts' standard errors.

sig2_beta   Numeric vector, exploratory variables' standard errors.

Tab_alpha   Data.frame, intercepts' summary.

Tab_beta    Data.frame, exploratory variables' summary.

Mat_alpha   Numeric matrix, intercepts' summary.

Mat_beta    Numeric matrix, exploratory variables' summary.
}
\description{
Estimate parameters and tuning parameter.
}
\examples{
n = 10
m = 5
d = 4
N = n*m
x = matrix(rnorm(d*N), ncol=d, nrow=N)
subj = rep(1:n, each=m)
alpha = rnorm(n)
beta = rnorm(d)
eps = rnorm(N)
y = as.vector(x \%*\% beta + rep(alpha, each=m) + eps)
m1 = pqr(x=x, y=y, subj=subj, tau=0.75, effect="lasso", c = 0)
m1$Tab_beta

}
\references{
Koenker, R. (2004) "Quantile regression for longitudinal data", J. Multivar. Anal., 91(1): 74-89, <doi:10.1016/j.jmva.2004.05.006>
}
