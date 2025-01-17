% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/R_matrix.R
\name{R_matrix}
\alias{R_matrix}
\title{Calculate \eqn{R} matrix}
\usage{
R_matrix(X.model)
}
\arguments{
\item{X.model}{numeric model matrix}
}
\value{
Numeric matrix \eqn{R}
}
\description{
Given model matrix \eqn{X.model}, return matrix \eqn{R} with elements \eqn{r_{ij} = a^2_{ij}/(a^2_{ii}a_{jj})},
where \eqn{a_{ij}} are the information matrix elements.
}
