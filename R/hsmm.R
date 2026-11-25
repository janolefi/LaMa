# Internal functions for HSMMs, approximated by HMMs on an extended state space
#
# Notation as in Koslik (2025, CSDA 209, 108171):
# - N states, state i is represented by an aggregate of Ni[i] states, M = sum(Ni)
# - I_low[i], I_high[i]: lowest and highest state of aggregate i
# - dm[[i]]: dwell-time distribution of state i, column r: probability that a stay lasts r time points
# - being in the r-th state of aggregate i means: r-th time point of this stay
#   (the highest state also holds all longer stays)


# Probabilities 1 - c_i(r) to stay for one more time point, c_i(r): hazard rates, eqs. (3) and (7)
# d: dwell-time distribution(s) of one state: vector, or matrix with one row for each time of entry
# F_before: F_i(r - 1), probability that the stay is over before its r-th time point (computed if not given)
# returns a matrix with the same rows and columns
hsmm_stay <- function(d, eps = 1e-10, F_before = NULL){
  if(is.null(dim(d))) d <- matrix(d, nrow = 1)
  Ni <- ncol(d)
  if(is.null(F_before)) {
    before <- matrix(0, Ni, Ni)
    before[upper.tri(before)] <- 1
    F_before <- d %*% before
  }
  hazard <- max2(d, eps) / (1 - F_before + eps / 2)
  max2(1 - hazard, 0)
}

# Stay probabilities of all M states for the steps into the time points `into` (one row each)
# the r-th state of an aggregate at time into - 1 belongs to a stay entered at into - r
# stay: list of matrices from hsmm_stay(), one row for each time of entry
# L: cycle length, for periodic models
hsmm_lag <- function(stay, into, L = NULL){
  columns <- list()
  for(stay_i in stay) {
    for(r in seq_len(ncol(stay_i))) {
      entry <- into - r
      if(!is.null(L)) entry <- (entry - 1) %% L + 1
      columns <- c(columns, list(stay_i[entry, r]))
    }
  }
  do.call(cbind, columns)
}

# Checks that omega (matrix or array) has zero diagonals: states are only left, not re-entered
# possible without AD and, while taping, for constant omega (the usual omega from tpm_emb() is AD then)
check_omega <- function(omega){
  if(inherits(omega, "advector")) return(invisible())
  diagonals <- omega[slice.index(omega, 1) == slice.index(omega, 2)] # of all slices
  if(anyNA(diagonals) || any(diagonals != 0)) {
    stop("omega must have zeros on the diagonal (conditional transition probabilities, given that the state is left), see tpm_emb().")
  }
}

# One step of the forward algorithm, phi %*% Gamma without building Gamma
# phi: vector of length M, (scaled) forward variables at t - 1
# stay: vector of length M, probabilities to stay for one more time point
# omega: conditional transition probabilities (zero diagonal)
hsmm_step <- function(phi, stay, omega, Ni){
  "[<-" <- ADoverload("[<-")
  N <- length(Ni)
  I_high <- cumsum(Ni)
  I_low <- I_high - Ni + 1

  kept <- phi * stay        # mass that stays in its aggregate
  gone <- phi * (1 - stay)  # mass whose stay ends

  # mass that leaves each aggregate, and where it arrives
  left <- numeric(N)
  for(i in 1:N) left[i] <- sum(gone[I_low[i]:I_high[i]])
  arrive <- as.vector(left %*% omega)

  phi_new <- numeric(sum(Ni))
  for(i in 1:N) {
    low <- I_low[i]
    high <- I_high[i]
    phi_new[low] <- arrive[i]                                      # new stays enter the lowest state
    if(high > low) phi_new[(low + 1):high] <- kept[low:(high - 1)] # the others move up by one
    phi_new[high] <- phi_new[high] + kept[high]                    # the highest state waits
  }
  phi_new
}

# Transition matrix of the extended state space (M x M), eqs. (1)-(2) and (5)-(6)
# only needed for stationary distributions and for users, not for the forward algorithm
# stay, omega, Ni: as in hsmm_step()
hsmm_Gamma <- function(stay, omega, Ni){
  "[<-" <- ADoverload("[<-")
  M <- sum(Ni)
  state <- rep(seq_along(Ni), Ni) # aggregate of each of the M states
  I_high <- cumsum(Ni)
  I_low <- I_high - Ni + 1

  Gamma <- matrix(0, M, M)
  # the stay ends: to the lowest state of another aggregate
  Gamma[, I_low] <- (1 - stay) * omega[state, , drop = FALSE]
  # the stay continues: one state up, or wait in the highest state
  up <- pmin(seq_len(M) + 1, I_high[state])
  Gamma[cbind(seq_len(M), up)] <- stay
  Gamma
}
