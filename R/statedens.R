#' State-dependent densities of one data stream
#'
#' Builds the matrix of state-dependent densities (or probabilities) of one data stream,
#' as used in the \code{allprobs} matrix of \code{\link{forward}}.
#' Several data streams that are conditionally independent given the state are combined by multiplying their matrices.
#'
#' @details
#' Each parameter in \code{...} can be
#' \itemize{
#'   \item a single value, used for all states,
#'   \item a vector of length \code{nStates}, one value per state, or
#'   \item a matrix of dimension \code{c(nObs, nStates)}, one value per observation and state, e.g. from a regression on covariates.
#' }
#' Missing observations contribute nothing to the likelihood, so their rows are 1.
#' As each data stream gets its own matrix, a missing value in one stream does not remove the observations of the others.
#'
#' Writing the loop over states by hand is equally valid and gives the same result. \code{statedens()} is just a shortcut.
#'
#' @param x vector of observations of length \code{nObs}, possibly containing \code{NA}s
#' @param dens density function, e.g. \code{dnorm}, \code{\link{dgamma2}} or \code{\link{dvm}}, with the observations as first argument
#' @param ... parameters of \code{dens}, see \sQuote{Details}
#' @param nStates number of states. Only needed if all parameters are single values.
#'
#' @return matrix of dimension \code{c(nObs, nStates)}
#' @export
#'
#' @examples
#' step <- trex$step[1:200]
#' angle <- trex$angle[1:200]
#' mu <- c(0.3, 2.5); sigma <- c(0.2, 1.5); kappa <- c(0.2, 1.5)
#'
#' # one data stream
#' allprobs <- statedens(step, dgamma2, mu, sigma)
#'
#' # the same by hand
#' allprobs2 <- matrix(1, length(step), 2)
#' ind <- which(!is.na(step))
#' for(j in 1:2) allprobs2[ind, j] <- dgamma2(step[ind], mu[j], sigma[j])
#' all.equal(allprobs, allprobs2)
#'
#' # two data streams, mean of the turning angle the same in both states
#' allprobs <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)
statedens <- function(x, dens, ..., nStates = NULL) {
  "[<-" <- ADoverload("[<-")
  pars <- list(...)
  # number of states from the parameters: columns of matrices, else the longest vector
  if(is.null(nStates)) nStates <- max(sapply(pars, function(p) if(is.matrix(p)) ncol(p) else length(p)))

  ind <- which(!is.na(x))
  allprobs <- matrix(1, length(x), nStates)
  for(j in 1:nStates) {
    # parameters of state j at the observed rows
    pj <- lapply(pars, function(p) if(is.matrix(p)) p[ind, j] else if(length(p) == 1) p else p[j])
    allprobs[ind, j] <- do.call(dens, c(list(x[ind]), pj))
  }
  allprobs
}
