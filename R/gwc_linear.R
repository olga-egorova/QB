#' Calculate generalized word count for linear effects only, i.e. \eqn{b_k(k)}
#'
#' @description
#' Given design matrix \eqn{X.design} and a number of linear effects,
#' return the generalized word count. Special case of gwc function.
#' @param \eqn{X.design} Design matrix
#' @param nlinear Number of linear effects in the polynomial effect
#' @return Generalized word count
#' @export
#' @example


gwc_linear <- function(X.design, nlinear) {

  nfactors <- ncol(X.design)
  nlinear <- as.numeric(nlinear)

  if (nlinear <= 0) return(NA)
  if (nlinear > nfactors) {
    stop("Number of effects should be no more than the number of factors")
  }

  # combinations of factors that form linear effects
  idx_linear <- matrix(combn(x = 1:nfactors, m = nlinear), nrow = nlinear)

  # sum across all combinations of indices
  bk = sum(apply(idx_linear, 2, function(x) sum(apply(as.matrix(X.design[,x]), 1, prod)))^2)/
       (nrow(X.design))^2

  return (bk)

}
