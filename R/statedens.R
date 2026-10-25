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
#' Missing observations contribute nothing to the likelihood, so their rows are 1 (0 with \code{log = TRUE}).
#' As each data stream gets its own matrix, a missing value in one stream does not remove the observations of the others.
#'
#' With \code{log = TRUE}, the log-densities are returned, to be added over data streams and used with \code{forward(..., logspace = TRUE)}.
#' They are computed by \code{dens(..., log = TRUE)} if \code{dens} has an argument \code{log}, otherwise as \code{log(dens(...))}
#' (e.g. for densities from \code{\link{zero_inflate}}).
#'
#' In \code{obj$simulate()} of \code{RTMB}, \code{x} is registered as data and simulated by \code{\link{forward}}, see \code{\link{simulate.LaMaModel}}.
#'
#' Writing the loop over states by hand is equally valid and gives the same result. \code{statedens()} is just a shortcut.
#'
#' @param x vector of observations of length \code{nObs}, possibly containing \code{NA}s
#' @param dens density function, e.g. \code{dnorm}, \code{\link{dgamma2}} or \code{\link{dvm}}, with the observations as first argument
#' @param ... parameters of \code{dens}, see \sQuote{Details}
#' @param nStates number of states. Only needed if all parameters are single values.
#' @param log logical; if \code{TRUE}, log-densities are returned
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
#'
#' # on the log scale, for forward(..., logspace = TRUE)
#' logprobs <- statedens(step, dgamma2, mu, sigma, log = TRUE) +
#'   statedens(angle, dvm, 0, kappa, log = TRUE)
statedens <- function(x, dens, ..., nStates = NULL, log = FALSE) {
  # name of the data stream for pseudo_res(mod), which swaps dens for its CDF
  s <- .lama_env$stream <- stream_name(substitute(x))
  on.exit(.lama_env$stream <- NULL)
  # register x as data for obj$simulate(), as OBS(x) would: there it returns a simref to be filled by forward()
  xobs <- eval(call("OBS", as.name(s)), list2env(setNames(list(x), s)))
  if(inherits(xobs, "simref")) x <- xobs
  pars <- list(...)
  # number of states from the parameters: columns of matrices, else the longest vector
  if(is.null(nStates)) nStates <- max(sapply(pars, function(p) if(is.matrix(p)) ncol(p) else length(p)))

  # obj$simulate(): random effects that are not simulated yet cannot enter the parameters.
  # Their transformations stay "simref", but cbind() of them gives a matrix without one row per observation.
  unsimulated <- vapply(pars, function(p) inherits(p, "simref") ||
                          (inherits(x, "simref") && is.matrix(p) && nrow(p) != length(x)), logical(1))
  if(any(unsimulated)) {
    stop("Automatic simulation requires the likelihood to follow the model hierarchy: random effects first, ",
         "then data given those random effects. Evaluate the density of each random effect ",
         "(u %~% dnorm(...) or nll <- nll - sum(dnorm(u, ...))) before using it in the parameters of statedens().")
  }

  # stream for plot(mod), reported by forward()
  .lama_env$streams[[s]] <- list(x = x, dens = dens, name = paste(deparse(substitute(dens)), collapse = ""),
                                 pars = pars, nStates = nStates)
  # simulation by obj$simulate(): forward() simulates x, nothing to evaluate here
  if(inherits(x, "simref")) return(matrix(if(log) 0 else 1, length(x), nStates))
  dens_matrix(x, dens, pars, nStates, log)
}

# matrix of state-dependent densities of x, as returned by statedens()
dens_matrix <- function(x, dens, pars, nStates, log = FALSE) {
  "[<-" <- ADoverload("[<-")
  # log-densities from dens itself if it can, else log() afterwards
  dlog <- log && "log" %in% names(formals(dens))

  ind <- which(!is.na(x))
  allprobs <- matrix(if(log) 0 else 1, length(x), nStates)
  for(j in 1:nStates) {
    d <- do.call(dens, c(list(x[ind]), state_pars(pars, j, ind), if(dlog) list(log = TRUE)))
    allprobs[ind, j] <- if(log && !dlog) base::log(d) else d
  }
  allprobs
}

# parameters of state j at the given rows: matrices have one column per state, vectors one entry
state_pars <- function(pars, j, rows) {
  lapply(pars, function(p) if(is.matrix(p)) p[rows, j] else if(length(p) == 1) p else p[j])
}
