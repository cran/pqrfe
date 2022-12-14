% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{loss_er}
\alias{loss_er}
\title{Loss expectile regression}
\usage{
loss_er(beta, x, y, tau, N, d)
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
eta Numeric, sum of expectile regression
}
\description{
This function returns the core of expectile regression to be minimized
}
