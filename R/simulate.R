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
