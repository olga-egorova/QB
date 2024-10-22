#' Calculate \eqn{R} matrix
#' @description Given model matrix \eqn{X.model}, return matrix \eqn{R} with elements \eqn{r_{ij} = a^2_{ij}/(a^2_{ii}a_{jj})},
#' where \eqn{a_{ij}} are the information matrix elements.
#'
#' @param X.model numeric model matrix
#' @return Numeric matrix \eqn{R}
#' @export

R_matrix <- function(X.model) {
  XtX <- crossprod(as.matrix(X.model))
  npar <- nrow(XtX)
  R <- XtX^2
  dXtX <- matrix(rep(diag(XtX), npar), nrow = npar, bycol = TRUE)
  R <- R/t(dXtX)/(dXtX^2)  # divide by diagonal and squared diagonal elements

  return(R)
}
