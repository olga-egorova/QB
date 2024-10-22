#' Calculate Q_B criterion for the 1st order maximal model
#' @description Returns the \eqn{Q_B} value for the 1st order maximal model, as well as weights \eqn{\xi}
#' and generalized word count (\eqn{gwc})
#' @param \eqn{X.design} Design matrix
#' @param \eqn{pi.linear} Prior probability of linear terms being in the best model
#' @export
#' @return

QB_I <- function(X.design, pi.linear) {

  nf <- ncol(X.design)
  xi.1st.order <- xi_1(nfactors = nf, pi = pi.linear)
  gwc.1st.order <- vapply(c(1,2), FUN = gwc_linear, X.design = X.design,
                          FUN.VALUE = numeric(1))

  qb <- sum(c(xi.1st.order$xi.1, 2*xi.1st.order$xi.2)*gwc.1st.order)

  return(list(xi = xi.1st.order, gwc = gwc.1st.order, qb = qb))

}
