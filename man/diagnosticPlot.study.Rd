% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_study-fn.R
\name{diagnosticPlot.study}
\alias{diagnosticPlot.study}
\title{Wrapper function for performance::check_model().
 Used to be a custom funciton, but the check_model function just seems so much better.
 Why reinvent the wheel?}
\usage{
diagnosticPlot.study(object, ...)
}
\arguments{
\item{object}{study object}

\item{...}{forwarded to performance::check_model()}
}
\value{
output of performance::check_model()
}
\description{
Wrapper function for performance::check_model().
 Used to be a custom funciton, but the check_model function just seems so much better.
 Why reinvent the wheel?
}
\author{
Stefan Thoma
}
