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
  x[is.na(st$x)] <- NA # keep the missing observations of the data
  x
}
