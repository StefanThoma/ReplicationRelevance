% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/2_study-fn.R
\name{create_pub_tables.study}
\alias{create_pub_tables.study}
\title{creates latex table for publication.}
\usage{
create_pub_tables.study(object, type = "wald", path = NULL)
}
\arguments{
\item{object}{a study object}

\item{type}{which type of difference table should be returned.
defaults to "wald"}

\item{path}{where table should be saved.}
}
\value{
last called table in latex format. It does save BOTH tables (if available) to path.
}
\description{
The resulting table may still need to be adjusted, but serves as a sceleton.
}
