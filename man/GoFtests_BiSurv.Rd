% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GoFtests_BiSurv.R
\name{GoFtests_BiSurv}
\alias{GoFtests_BiSurv}
\title{PIOS and IR tests of copula function for bivariate survival time}
\usage{
GoFtests_BiSurv(
  data,
  copula.fam,
  control = list(yes.PIOSexact = FALSE, yes.boot = TRUE, nboot = 500, seed1 = 1234)
)
}
\arguments{
\item{data}{a data frame with four columns, in which the first two columns are censored event times, and the next two columns are censoring indicator variables}

\item{copula.fam}{a character indicating which one of the following copula families: "clayton", "frank", "gumbel", and "normal"}

\item{control}{a list of the following components: \code{yes.PIOSexact}, \code{yes.boot}, \code{nboot}, \code{seed1}. \code{yes.PIOSexact} is a logical value indicating whether to calculate the exact PIOS test statistic, and the default value is \code{FALSE}. \code{yes.boot} is a logical value indicating whether to implement the bootstrap procedure. \code{nboot} is the number of bootstrap samples.  \code{seed1} is the seed for generating the bootstrap samples.}
}
\value{

}
\description{

}
