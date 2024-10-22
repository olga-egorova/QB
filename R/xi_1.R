#' Calculate \eqn{\xi} for the 1st order maximal model
#' @description Evaluate \eqn{\xi} - sum of probabilities of each model being the best, over the submodels of the
#' first-order maximal model containing (1) main effects of each factor (2) pairs of main effects
#' @param nfactors Number of factors, each is of two levels
#' @param pi Probability of each main effect being present in the best model
#' @returns A list of weights: \eqn{\xi.1} for main effects and {\xi.2} for pairs of main effects
#' @export

xi_1 <- function(nfactors, pi){

  num.of.terms <- 1:nfactors     # number of linear terms in the submodels of the 1st order maximal model
  probs <- vapply(num.of.terms,
                  FUN = function(x, pi, nf) pi^x*(1-pi)^(nf - x), pi = pi, nf = nfactors,
                  FUN.VALUE = numeric(1))
  num.of.models.1 <- vapply(num.of.terms,                # num. of models containing one ME
                            FUN = function(x, nf) nf*choose(n = nf - 1, k = x - 1),
                            nf = nfactors,
                            FUN.VALUE = numeric(1))
  xi.1 <- sum(probs*num.of.models.1)

  num.of.models.2 <- vapply(num.of.terms[-1],           # num. of models containing two MEs
                            FUN = function(x, nf)
                              choose(nf,2)*choose(n = nf - 2, k = x - 2), nf = nfactors,
                            FUN.VALUE = numeric(1))
  xi.2 <- sum(probs[-1]*num.of.models.2)

  return (list(xi.1 = xi.1, xi.2 = xi.2))
}
