% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{loss_qr}
\alias{loss_qr}
\title{Loss quantile regression}
\usage{
loss_qr(beta, x, y, tau, N, d)
}
\arguments{
\item{beta}{initial values}

\item{x}{design matrix}

\item{y}{vector output}

\item{tau}{percentile}

\item{N}{sample size}

\item{d}{columns of x}
}
\value{
eta Numeric, sum of quantile regression
}
\description{
This function returns the core of quantile regression to be minimized
}
