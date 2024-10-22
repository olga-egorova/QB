#' Calculate generalised word count, for linear and quadratic effects
#'
#' @description
#' Given design matrix \eqn{X.design} and a named list of polynomial effects,
#' return the generalised word count.
#' Generalised case where quadratic effects of all factors can be included in the gwc.
#' @param \eqn{X.design} Design matrix
#' @param effects Named list of numbers of polynomial effects: 'linear', 'quadratic'
#' @return Generalised word count -- a number
#' @export
#' @example



# Add NaN checks, NA (e.g. for 2-level designs only)

gwc <- function(X.design, effects) {

  nfactors <- ncol(X.design)
  k <- sum(unlist(effects))  # k -- number of factors in the gwc

  if (k > nfactors) {
    stop("Sum of number of effects should be no more than the number of factors")
    }
  nlinear <- as.numeric(effects$linear)
  nquadratic <- as.numeric(effects$quadratic)

  if (nquadratic == 0) {return (gwc_linear(X.design, nlinear))}

  # ## Option 1.  try with gtools: all permutations of all subsets of k out of nfactors. Takes longer, a bit more elegant.
  # idx_effects <- gtools::permutations(n = choose(n = nfactors, k = k), r = k)  #gtools::
  #
  # # order factors within the 'linear' and 'quadratic' parts of each of the subset. Adding higher terms is easier here -- just another "order(..)"
  # idx_effects <- apply(idx_effects, 1, function(x) c(sort(x[1:nlinear]), sort(x[-(1:nlinear)])))
  # idx_effects <- unique(idx_effects, MARGIN=2)
  # idx_linear <- idx_effects[1:nlinear,]
  # idx_quadratic <- idx_effects[nlinear+(1:nquadratic),]


  ## Option 2. Expanding the grid of all combinations

  lvls <- expand.grid("linear" = 1:choose(n = nfactors, k = nlinear),
                      "quadratic" = 1:choose(n = nfactors - nlinear, k = nquadratic))
  if (nlinear > 0){
    idx_linear <- matrix(combn(x = 1:nfactors, m = nlinear)[, lvls$linear], nrow = nlinear)
  } else {
    idx_linear <-  matrix(0, nrow = 1, ncol = nrow(lvls))
  }

  idx_avail <- matrix(0, ncol = ncol(idx_linear), nrow = nfactors - nlinear)
  idx_avail[] <- apply(idx_linear, 2, FUN = function(x) setdiff(1:nfactors, x))   # factors available for higher order effects
  idx_quadratic <- matrix(0, ncol = ncol(idx_linear), nrow = nquadratic)     # indices of the available factors to form quadratic effects (from the available ones)
  # from each set of available indices choose the ones indicated by lvls$quadratic
  lvl_quadratic <- matrix(0, nrow = nquadratic, ncol = ncol(idx_avail))
  lvl_quadratic[] <- combn(x = 1:(nfactors - nlinear), m = nquadratic)[, lvls$quadratic]
  idx_quadratic[] <- apply(rbind(idx_avail, lvl_quadratic), 2, FUN = function(x) x[x[-(1:nrow(idx_avail))]])

  # sum across all combinations of indices. check for either linear or quadratic equal to 1

  if (nlinear > 0) {
    bk <- sum(
      apply(apply(idx_linear, 2,
                  function(x) apply(as.matrix(sqrt(3/2)*X.design[,x]), 1, prod))*         # products of linear and
              apply(idx_quadratic, 2,
                    function(x) (apply(as.matrix(sqrt(1/2)*(3*X.design[,x]^2 - 2)), 1, prod))^2),  # quadratic terms for each index
            2, function(x) (sum(x))^2)                                                           # squared sum = j(R1,.., Rp)^2 across indeces
    ) /(nrow(X.design))^2                                                                      # final sum, divided by Nruns^2
  } else {
    bk <- sum(
      apply(apply(idx_quadratic, 2,
                    function(x) (apply(as.matrix(sqrt(1/2)*(3*X.design[,x]^2 - 2)), 1, prod))^2),  # quadratic terms for each index
              2, function(x) (sum(x))^2)                                                           # squared sum = j(R1,.., Rp)^2 across indeces
    ) /(nrow(X.design))^2
  }


  return (bk)

}
