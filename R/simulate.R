# Simulation from a fitted model, on the observed design: tracks, their lengths and
# the reported delta, Gamma and parameters are kept fixed

# states from the output of prepare_decoding():
# first state of each track ~ delta, then Z_t ~ Gamma[Z_(t-1), , gidx[t]]
sim_states <- function(inp) {
  nObs <- length(inp$gidx)
  N <- ncol(inp$delta)
  # cumulative probabilities: cumGamma[j, i, k] = P(next state <= j | current state i, slice k)
  cumGamma <- apply(inp$Gamma, c(1, 3), cumsum)
  cumDelta <- t(apply(inp$delta, 1, cumsum))
  isStart <- seq_len(nObs) %in% inp$start
  track <- cumsum(isStart)

  u <- runif(nObs)
  Z <- integer(nObs)
  for(t in seq_len(nObs)) {
    p <- if(isStart[t]) cumDelta[track[t], ] else cumGamma[, Z[t - 1], inp$gidx[t]]
    Z[t] <- min(sum(u[t] > p) + 1L, N) # min() guards against rounding in the cumulative sums
  }
  Z
}

# random number generator of a stream: rdens if given, else r<family> for the density d<family> from statedens()
get_rdens <- function(st, rdens = NULL) {
  if(!is.null(rdens)) return(match.fun(rdens))
  r <- paste0("r", sub("^d", "", sub("^.*::", "", st$name)))
  env <- environment(st$dens)
  f <- get0(r, envir = if(is.null(env)) globalenv() else env, mode = "function")
  if(is.null(f)) stop("No random number generator ", r, "() found for ", st$name,
                      ", provide one via rdens.")
  f
}

# observations of one stream given the states: draws per state with the parameters of the rows in that state
sim_stream <- function(st, Z, rdens = NULL) {
  r <- get_rdens(st, rdens)
  x <- rep(NA_real_, length(Z))
  for(j in seq_len(st$nStates)) {
    rows <- which(Z == j)
    if(length(rows)) x[rows] <- do.call(r, c(list(length(rows)), state_pars(st$pars, j, rows)))
  }
  # keep the missing observations of the data, unless the stream is all NA (placeholder to simulate from)
  if(!all(is.na(st$x))) x[is.na(st$x)] <- NA
  x
}

#' Simulate from a fitted model
#'
#' Simulates states and observations from a model fitted with \code{\link{forward}} and \code{\link{statedens}},
#' on the observed design: the tracks and their lengths, the initial distribution \code{delta},
#' the transition probabilities \code{Gamma} and the parameters of the state-dependent densities
#' are those reported by the model, i.e. conditional on the observed covariates.
#' Simulating for new covariates is not possible; instead, rebuild the likelihood function for the new data
#' and call \code{\link{report}} on the resulting object.
#'
#' @details
#' The states are drawn as in the decoding functions, so that all types of \code{Gamma} that \code{\link{forward}} accepts
#' (homogeneous, track-specific, time-varying, periodic) work, as does a continuous-time model with one transition matrix
#' per time interval (the simulation is then conditional on the observed time points).
#' Each data stream is drawn from the random number generator belonging to its density, e.g. \code{rgamma2} for \code{dgamma2}.
#' Missing observations of the data stay missing, except for a stream that is missing completely: such a placeholder
#' (e.g. \code{step = rep(NA, n)}) lets you simulate from a likelihood function without data, see the examples. Models with random effects are simulated conditional on the estimated random effects.
#'
#' @param object model object from \code{\link{report}} or \code{\link{qreml}} that was built with \code{\link{statedens}}
#' @param nsim number of simulations
#' @param seed optional seed for the random number generator
#' @param rdens optional named list of random number generators of the data streams, for densities without one named
#' by replacing the leading \code{d} with \code{r} (e.g. from \code{\link{zero_inflate}}). The first argument is the number of draws,
#' the others are the parameters as in the density.
#' @param ... currently not used
#'
#' @return data frame with the columns \code{trackID} (if used), \code{state} and one column per data stream, for \code{nsim = 1},
#' else a list of \code{nsim} such data frames.
#' @export
#' @importFrom stats simulate
#'
#' @examples
#' nll <- function(par) {
#'   getAll(par, dat)
#'   Gamma <- tpm(eta)
#'   allprobs <- statedens(step, dgamma2, exp(log_mu), exp(log_sigma)) *
#'     statedens(angle, dvm, 0, exp(log_kappa))
#'   -forward(stationary(Gamma), Gamma, allprobs)
#' }
#' dat <- list(step = trex$step[1:500], angle = trex$angle[1:500])
#' par <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)),
#'             log_kappa = log(c(0.2, 1.5)))
#' obj <- MakeADFun(nll, par, silent = TRUE)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' mod <- report(obj)
#'
#' sim <- simulate(mod, seed = 1)
#' head(sim)
#'
#' # simulate without data: all-NA streams as placeholders, parameters of your choice, then refit
#' dat$step <- dat$angle <- rep(NA, 500)
#' obj0 <- MakeADFun(nll, par, silent = TRUE)
#' sim <- simulate(report(obj0), seed = 1)
#' dat$step <- sim$step; dat$angle <- sim$angle
#' obj1 <- MakeADFun(nll, par, silent = TRUE)
#' opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
simulate.LaMaModel <- function(object, nsim = 1, seed = NULL, rdens = NULL, ...) {
  # keep the random number stream of the session, as stats::simulate() methods do
  if(!is.null(seed)) {
    if(!exists(".Random.seed", envir = globalenv())) runif(1)
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  st <- get_streams(object)
  inp <- prepare_decoding(mod = object)

  sims <- lapply(seq_len(nsim), function(i) {
    Z <- sim_states(inp)
    x <- lapply(names(st), function(s) sim_stream(st[[s]], Z, rdens[[s]]))
    names(x) <- names(st)
    out <- data.frame(state = Z, x)
    if(!is.null(object$trackID)) out <- cbind(trackID = object$trackID, out)
    out
  })
  if(nsim == 1) sims[[1]] else sims
}
