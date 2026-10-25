# Simulation of states and data streams from models built with forward() and statedens():
# - simulate(mod): on the observed design, with delta, Gamma and parameters as reported
# - simulate(mod, random = TRUE) and obj$simulate(): forward() simulates inside RTMB's simulation
#   (sim_forward()), after RTMB has drawn the random effects

# states from the output of prepare_decoding():
# first state of each track ~ delta, then Z_t ~ Gamma[Z_(t-1), , gidx[t]]
sim_states <- function(inp) {
  N <- ncol(inp$delta)
  isStart <- seq_along(inp$gidx) %in% inp$start
  track <- cumsum(isStart)
  Z <- integer(length(inp$gidx))
  for(t in seq_along(Z)) {
    prob <- if(isStart[t]) inp$delta[track[t], ] else inp$Gamma[Z[t - 1], , inp$gidx[t]]
    Z[t] <- sample.int(N, 1, prob = prob)
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
sim_stream <- function(st, Z, rdens = NULL, fill_na = FALSE) {
  r <- get_rdens(st, rdens)
  x <- rep(NA_real_, length(Z))
  for(j in seq_len(st$nStates)) {
    rows <- which(Z == j)
    if(length(rows)) x[rows] <- do.call(r, c(list(length(rows)), state_pars(st$pars, j, rows)))
  }
  # keep the missing observations of the data, unless asked to fill them or the stream is all NA (placeholder to simulate from)
  if(!fill_na && !all(is.na(st$x))) x[is.na(st$x)] <- NA
  x
}

# forward() in obj$simulate(): simulate the states and fill the data of the streams (simref objects)
# with draws from the streams, given the delta, Gamma and parameters of this call
sim_forward <- function(streams, delta, Gamma, allprobs, trackID, logspace, tod) {
  inp <- prepare_forward(delta, Gamma, allprobs, trackID, logspace, gidx = tod)
  states <- sim_states(inp)
  for(s in names(streams)) {
    st <- streams[[s]]
    x <- st$x
    if(!inherits(x, "simref")) next
    assign(s, x) # RTMB finds the simref by its name in the frames of the call stack, to return the simulated data
    st$x <- x$orig # the data itself, so that sim_stream() keeps its missing values
    x$reverse.update(sim_stream(st, states, .lama_env$sim$rdens[[s]], isTRUE(.lama_env$sim$fill_na)), seq_along(x))
  }
  REPORT(states)
  0
}

#' Simulate from a fitted model
#'
#' Simulates states and observations from a model fitted with \code{\link{forward}} and \code{\link{statedens}}.
#' By default, this is done on the observed design: the tracks and their lengths, the initial distribution \code{delta},
#' the transition probabilities \code{Gamma} and the parameters of the state-dependent densities
#' are those reported by the model, i.e. conditional on the observed covariates.
#' Simulating for new covariates is not possible; instead, rebuild the likelihood function for the new data
#' and call \code{\link{lama}} on the resulting object.
#'
#' @details
#' The states are drawn as in the decoding functions, so that all types of \code{Gamma} that \code{\link{forward}} accepts
#' (homogeneous, track-specific, time-varying, periodic) work, as does a continuous-time model with one transition matrix
#' per time interval (the simulation is then conditional on the observed time points).
#' Each data stream is drawn from the random number generator belonging to its density, e.g. \code{rgamma2} for \code{dgamma2}.
#'
#' Missing observations of the data stay missing (unless \code{fill_na = TRUE}), except for a stream that is missing completely: such a placeholder
#' (e.g. \code{step = rep(NA, n)}) lets you simulate from a likelihood function without data, see the examples.
#'
#' Models with random effects are by default simulated conditional on the estimated random effects.
#' With \code{random = TRUE}, \code{RTMB} draws the random effects first, which requires that the density of each random effect
#' is evaluated in the likelihood function \emph{before} the random effect is used, either by \code{u \%~\% dnorm(...)}
#' or by \code{nll <- nll - sum(dnorm(u, ..., log = TRUE))} (or another density).
#' The states and data are then simulated given these draws, so that every simulation has new random effects.
#' This also requires that \code{allprobs} is the product of the data streams built by \code{\link{statedens}},
#' with each data stream passed once.
#'
#' @param object model object from \code{\link{lama}} or \code{\link{qreml}} that was built with \code{\link{statedens}}
#' @param nsim number of simulations
#' @param seed optional seed for the random number generator
#' @param rdens optional named list of random number generators of the data streams, for densities without one named
#' by replacing the leading \code{d} with \code{r} (e.g. from \code{\link{zero_inflate}}). The first argument is the number of draws,
#' the others are the parameters as in the density.
#' @param random logical; if \code{TRUE}, \code{obj$simulate()} of \code{RTMB} is used, so that random effects are drawn
#' from their distribution first, see \sQuote{Details}. The data streams are those passed to \code{\link{statedens}}.
#' Default is \code{FALSE}, which keeps the reported (estimated) random effects.
#' @param fill_na logical; if \code{TRUE}, missing observations of the data are filled in by simulated values
#' (their parameters need to be defined, e.g. covariates must not be missing). Default is \code{FALSE}, which keeps them missing.
#' Streams that are missing completely are always simulated.
#' @param ... currently not used
#'
#' @return data frame with the columns \code{trackID} (index of the track, if used), \code{state} and one column per data stream,
#' for \code{nsim = 1}, else a list of \code{nsim} such data frames.
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
#' mod <- lama(obj)
#'
#' sim <- simulate(mod, seed = 1)
#' head(sim)
#'
#' # simulate without data: all-NA streams as placeholders, parameters of your choice, then refit
#' dat$step <- dat$angle <- rep(NA, 500)
#' obj0 <- MakeADFun(nll, par, silent = TRUE)
#' sim <- simulate(lama(obj0), seed = 1)
#' dat$step <- sim$step; dat$angle <- sim$angle
#' obj1 <- MakeADFun(nll, par, silent = TRUE)
#' opt1 <- nlminb(obj1$par, obj1$fn, obj1$gr)
#'
#' # random effects: RTMB draws them first, then the states and data are simulated given them
#' nll_re <- function(par) {
#'   getAll(par, dat)
#'   u %~% dnorm(0, exp(log_su)) # track-level random effect on the mean step length in state 1
#'   mu <- exp(cbind(log_mu[1] + u[id], log_mu[2]))
#'   Gamma <- tpm(eta)
#'   .nll <- .nll - forward(stationary(Gamma), Gamma, statedens(step, dgamma2, mu, exp(log_sigma)), id)
#'   .nll
#' }
#' dat <- list(step = trex$step[1:500], id = rep(1:5, each = 100))
#' par_re <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)),
#'                log_su = log(0.3), u = rep(0, 5))
#' obj_re <- MakeADFun(nll_re, par_re, random = "u", silent = TRUE)
#' sim_re <- simulate(lama(obj_re), random = TRUE, seed = 1)
simulate.LaMaModel <- function(object, nsim = 1, seed = NULL, rdens = NULL, random = FALSE, fill_na = FALSE, ...) {
  # keep the random number stream of the session, as stats::simulate() methods do
  if(!is.null(seed)) {
    if(!exists(".Random.seed", envir = globalenv())) runif(1)
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  # one simulation: list of the states and the data streams
  if(random) {
    obj <- get_obj(object)
    if(!length(obj$env$obs)) stop("No data streams found: build allprobs with statedens() to simulate with random = TRUE.")
    .lama_env$sim <- list(rdens = rdens, fill_na = fill_na) # for sim_forward(), called by forward() inside obj$simulate()
    on.exit(.lama_env$sim <- NULL, add = TRUE)
    sim_once <- function() {
      # RTMB draws the random effects, forward() the states and data given them
      sim <- obj$simulate(obj$env$last.par.best)
      if(is.null(sim$states)) stop("No states were simulated: the likelihood needs forward() with statedens().")
      list(states = sim$states, streams = sim[intersect(names(obj$env$obs), names(sim))])
    }
  } else {
    st <- get_streams(object)
    inp <- prepare_decoding(mod = object)
    sim_once <- function() {
      Z <- sim_states(inp)
      list(states = Z, streams = sapply(names(st), function(s) sim_stream(st[[s]], Z, rdens[[s]], fill_na), simplify = FALSE))
    }
  }

  sims <- lapply(seq_len(nsim), function(i) {
    sim <- sim_once()
    out <- data.frame(state = sim$states, sim$streams)
    if(!is.null(object$trackID)) out <- cbind(trackID = object$trackID, out)
    out
  })
  if(nsim == 1) sims[[1]] else sims
}
