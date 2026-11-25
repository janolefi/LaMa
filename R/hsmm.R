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
# first: first row that may be used (e.g. start of the track), used in place of earlier ones
hsmm_lag <- function(stay, into, L = NULL, first = 1){
  columns <- list()
  for(stay_i in stay) {
    for(r in seq_len(ncol(stay_i))) {
      entry <- into - r
      if(is.null(L)) {
        entry <- pmax(entry, first)
      } else {
        entry <- (entry - 1) %% L + 1
      }
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


# Preparation of the forward algorithm -------------------------------------

# Brings all inputs into the form of the forward algorithm, like prepare_forward() for HMMs
# - dm: list of N dwell-time distributions, each a vector or a matrix with 1, nTracks, nObs or L rows
# - omega: matrix or array with 1, nTracks, nObs, nObs - 1 or L slices, as Gamma in forward()
#   (L: cycle length, if tod is given; a stay entered at t has omega[, , t] and row t of dm)
# - delta: NULL or "stationary", or initial distribution of the N states, with nTracks rows if needed
# - use_all_obs, startInd: only for dm varying over time, where max(Ni) - 1 observations are dropped by default
# Returns a list with
# - stay[t, ]: stay probabilities of all M states for the step into t
# - omega[, , oidx[t]]: conditional transition probabilities for the step into t
# - delta: initial distribution of the M states (nTracks x M)
# - start, end: first and last observation used in each track
# - allprobs, shift, trackID, type: as in prepare_forward()
prepare_hsmm <- function(dm, omega, allprobs, trackID = NULL, delta = NULL, tod = NULL,
                         use_all_obs = FALSE, startInd = NULL, logspace = FALSE, eps = 1e-10){

  if(!is.matrix(allprobs)) stop("allprobs needs to be a matrix of dimension c(n, N).")
  nObs <- nrow(allprobs)
  N <- ncol(allprobs)
  tracks <- prepare_tracks(trackID, nObs)

  ## dwell times
  if(!is.list(dm) || length(dm) != N) stop("dm must be a list with one dwell-time distribution for each state.")
  dm <- lapply(dm, function(d) if(is.null(dim(d))) matrix(d, nrow = 1) else d)
  Ni <- vapply(dm, ncol, 1L)
  dm_index <- lapply(dm, function(d) hsmm_index(nrow(d), nObs, tracks, tod, "rows of dm"))
  dm_type <- vapply(dm_index, function(x) x$type, "")

  ## conditional transition probabilities
  check_omega(omega)
  omega <- AD(as.array(omega))
  if(length(dim(omega)) == 2) omega <- AD(array(omega, dim = c(dim(omega), 1)))
  if(!all(dim(omega)[1:2] == N)) stop("omega must be square with dimensions (N, N).")
  omega_index <- hsmm_index(dim(omega)[3], nObs, tracks, tod, "slices of omega")
  oidx <- omega_index$index

  ## periodic models: cycle length
  types <- c(dm_type, omega_index$type)
  periodic <- types == "periodic"
  L <- unique(c(vapply(dm, nrow, 1L), dim(omega)[3])[periodic])
  if(length(L) > 1) stop("All periodic elements of dm and omega must have the same cycle length.")
  if("periodic" %in% dm_type) check_tod(tod, L, tracks)

  ## stay probabilities for the step into each time point
  track_start <- rep(tracks$start, tracks$end - tracks$start + 1)
  stay <- hsmm_stay_into(dm, dm_type, lapply(dm_index, function(x) x$index), L, track_start, eps)

  ## first observation used in each track
  stationary <- is.null(delta) || is.character(delta)
  start <- tracks$start
  if("inhomogeneous" %in% dm_type) {
    if(use_all_obs && stationary) {
      stop("use_all_obs = TRUE needs delta: the stationary start would require dwell-time distributions ",
           "from before the first observation.")
    }
    if(!use_all_obs) {
      if(is.null(startInd)) startInd <- max(Ni)
      if(startInd < max(Ni)) stop("startInd must be at least the maximum aggregate size")
      start <- start + startInd - 1
      # the stationary distribution needs the step into the next observation
      if(any(start + stationary > tracks$end)) {
        stop("Too few observations in at least one track: with dm varying over time, ",
             "each track starts at the maximum aggregate size (or startInd).")
      }
    }
  }

  ## initial distribution
  if(!stationary) {
    delta <- hsmm_delta(delta, Ni, tracks$nTracks)
  } else if(is.character(delta) && !identical(delta, "stationary")) {
    stop("delta must be numeric or \"stationary\".")
  } else if(any(periodic)) {
    # transition matrices of the whole cycle, as in forward() with tod
    cycle <- function(nSlices) if(nSlices == 1) rep(1, L) else 1:L
    stay_cycle <- hsmm_stay_into(dm, dm_type, lapply(dm, function(d) cycle(nrow(d))), L, 1, eps)
    Gamma <- hsmm_Gammas(stay_cycle, omega[, , cycle(dim(omega)[3]), drop = FALSE], Ni)
    delta <- stationary_delta(Gamma, tod[start], periodic = TRUE)
  } else {
    # stationary for the first transition of each track, out of observation start
    into <- pmin(start + 1, tracks$end)
    if(all(types == "homogeneous")) into <- into[1] # all tracks are the same
    Gamma <- hsmm_Gammas(stay[into, , drop = FALSE], omega[, , oidx[into], drop = FALSE], Ni)
    delta <- stationary_delta(Gamma, rep_len(seq_along(into), tracks$nTracks))
  }

  # most general type of dm and omega
  ranking <- c("periodic", "inhomogeneous", "track-specific", "homogeneous")
  type <- ranking[min(match(types, ranking))]

  scaled <- scale_allprobs(allprobs, logspace)
  list(stay = stay, omega = omega, oidx = oidx, delta = delta, dm = dm, Ni = Ni,
       allprobs = scaled$allprobs, shift = scaled$shift,
       start = start, end = tracks$end, trackID = tracks$trackID, type = type)
}

# Which row of dm or slice of omega belongs to each observation? See slice_index()
# a single row or slice is also possible in periodic models
hsmm_index <- function(nSlices, nObs, tracks, tod, what){
  if(nSlices == 1) tod <- NULL
  slices <- slice_index(nSlices, nObs, tracks, tod, what)
  # no nObs - 1 for dm: row t belongs to stays entered at t
  if(is.null(slices) || (what == "rows of dm" && slices$type == "inhomogeneous" && nSlices != nObs)) {
    stop("Wrong number of ", what, ": must be 1, nTracks (track-specific), nObs (time-varying),",
         " or the cycle length if tod is given.")
  }
  slices
}

# tod must move on by one in each row of a track, as dwell times count rows
check_tod <- function(tod, L, tracks){
  expected <- c(0, tod[-length(tod)] %% L + 1)
  inside <- !seq_along(tod) %in% tracks$start # not the first row of a track
  if(any(tod[inside] != expected[inside])) {
    stop("With periodic dm, tod must advance by one in each row within a track. ",
         "Keep rows with missing observations (set allprobs to 1).")
  }
}

# Stay probabilities of all M states, one row for each step (as in hsmm_lag())
# type, index: type and row index of each element of dm, see slice_index()
hsmm_stay_into <- function(dm, type, index, L, first, eps){
  stay <- lapply(seq_along(dm), function(i) {
    stay_i <- hsmm_stay(dm[[i]], eps)
    if(type[i] == "inhomogeneous") {
      hsmm_lag(list(stay_i), index[[i]], first = first)
    } else if(type[i] == "periodic") {
      hsmm_lag(list(stay_i), index[[i]], L)
    } else { # the same for all stays of a track
      stay_i[index[[i]], , drop = FALSE]
    }
  })
  do.call(cbind, stay)
}

# Array of transition matrices, for the rows of stay and the slices of omega
hsmm_Gammas <- function(stay, omega, Ni){
  "[<-" <- ADoverload("[<-")
  M <- sum(Ni)
  Gamma <- array(0, dim = c(M, M, nrow(stay)))
  for(k in seq_len(nrow(stay))) Gamma[, , k] <- hsmm_Gamma(stay[k, ], omega[, , k], Ni)
  Gamma
}

# Initial distribution of the N states as matrix (nTracks x M): stays start at the first observation,
# so all mass is on the lowest states
hsmm_delta <- function(delta, Ni, nTracks){
  "[<-" <- ADoverload("[<-")
  N <- length(Ni)
  delta <- as.matrix(delta)
  if(ncol(delta) != N) delta <- t(delta)
  if(ncol(delta) != N || !nrow(delta) %in% c(1, nTracks)) {
    stop("delta must be a vector of length N or a matrix of dimension c(nTracks, N).")
  }
  if(nrow(delta) == 1) delta <- delta[rep(1, nTracks), , drop = FALSE]
  I_low <- cumsum(Ni) - Ni + 1
  delta_M <- matrix(0, nTracks, sum(Ni))
  delta_M[, I_low] <- delta
  delta_M
}
