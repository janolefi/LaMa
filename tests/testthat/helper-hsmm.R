# notation as in the paper (Koslik, 2025, CSDA 209, 108171)
# N states, aggregate sizes Ni, M = sum(Ni) states of the approximating chain
# state aggregate I_i runs from its lowest state I_low[i] to its highest state I_high[i]
# dm[[i]][t, r]: probability that a stay in state i entered at time t lasts r time points

# ---- exact references, written directly from the equations in the paper ----

# hazard rates c_i(1), ..., c_i(N_i) of one dwell-time distribution d_i, eqs. (3), (7)
ref_hazard <- function(d) {
  F_before <- c(0, cumsum(d))[seq_along(d)] # F_i(r - 1)
  ifelse(F_before < 1, d / (1 - F_before), 1)
}

# block matrix Gamma with diagonal blocks (1), (5) and off-diagonal blocks (2), (6)
# hazard: list of N vectors, hazard[[i]][r] applies to the r-th state of aggregate i
ref_Gamma <- function(omega, hazard) {
  Ni <- lengths(hazard)
  blocks <- lapply(seq_along(Ni), function(i) {
    row_i <- lapply(seq_along(Ni), function(j) {
      block <- matrix(0, Ni[i], Ni[j])
      if(i != j) {
        block[, 1] <- omega[i, j] * hazard[[i]]
      } else {
        for(r in seq_len(Ni[i])) block[r, min(r + 1, Ni[i])] <- 1 - hazard[[i]][r]
      }
      block
    })
    do.call(cbind, row_i)
  })
  do.call(rbind, blocks)
}

# Gamma^(t): state r of aggregate i has the hazard of a stay entered at t - r + 1
# L: cycle length of a periodic model
ref_Gamma_t <- function(t, omega, dm, L = NULL) {
  hazard <- lapply(dm, function(d) {
    sapply(seq_len(ncol(d)), function(r) {
      entry <- t - r + 1
      if(!is.null(L)) entry <- (entry - 1) %% L + 1
      ref_hazard(d[entry, ])[r]
    })
  })
  ref_Gamma(omega, hazard)
}

# initial distribution: stationary or a fresh stay at the start (mass on the lowest states)
ref_delta_stat <- function(Gamma) solve(t(diag(nrow(Gamma)) - Gamma + 1), rep(1, nrow(Gamma)))
ref_delta_fresh <- function(delta, Ni) {
  I_low <- cumsum(Ni) - Ni + 1
  replace(numeric(sum(Ni)), I_low, delta)
}

# ---- direct forward step (to be moved to R/, then delete it here) ----

# probability to stay for one more time point, 1 - c_i(r), same form as d
# d: matrix, one row per entry time, one column per dwell time
hsmm_stay <- function(d, eps = 1e-10) {
  upper <- matrix(0, ncol(d), ncol(d))
  upper[upper.tri(upper)] <- 1
  F_before <- d %*% upper # F_i(r - 1)
  hazard <- max2(d, eps) / (1 - F_before + eps / 2)
  max2(1 - hazard, 0)
}

# stay probabilities of all M states at time t: state r of aggregate i was entered at t - r + 1
# stay: list of matrices from hsmm_stay()
hsmm_stay_t <- function(t, stay, L = NULL) {
  "c" <- ADoverload("c")
  stay_t <- lapply(stay, function(stay_i) {
    r <- seq_len(ncol(stay_i))
    entry <- t - r + 1
    if(!is.null(L)) entry <- (entry - 1) %% L + 1
    stay_i[cbind(entry, r)]
  })
  do.call(c, stay_t)
}

# one step phi %*% Gamma without building Gamma
hsmm_step <- function(phi, stay_t, omega, Ni) {
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  N <- length(Ni)
  I_high <- cumsum(Ni)
  I_low <- I_high - Ni + 1

  kept <- phi * stay_t         # mass that stays in its aggregate
  gone <- phi * (1 - stay_t)   # mass whose stay ends

  # mass that leaves each aggregate and where it arrives
  left <- AD(numeric(N))
  for(i in 1:N) left[i] <- sum(gone[I_low[i]:I_high[i]])
  arrive <- as.vector(left %*% omega)

  phi_new <- AD(numeric(sum(Ni)))
  for(i in 1:N) {
    low <- I_low[i]
    high <- I_high[i]
    phi_new[low] <- arrive[i]                                      # new stays enter the lowest state
    if(high > low) phi_new[(low + 1):high] <- kept[low:(high - 1)] # others move up by one
    phi_new[high] <- phi_new[high] + kept[high]                    # highest state waits
  }
  phi_new
}

# forward algorithm built on hsmm_step()
# stay_into(t), omega_into(t): for the step from t - 1 into t, delta: length M
hsmm_forward <- function(delta, stay_into, omega_into, allprobs, Ni, start = 1) {
  state <- rep(seq_along(Ni), Ni) # aggregate of each of the M states
  foo <- delta * allprobs[start, state]
  l <- log(sum(foo))
  phi <- foo / sum(foo)
  for(t in seq_len(nrow(allprobs))) {
    if(t <= start) next
    foo <- hsmm_step(phi, stay_into(t), omega_into(t), Ni) * allprobs[t, state]
    l <- l + log(sum(foo))
    phi <- foo / sum(foo)
  }
  l
}

# ---- known bugs ----

# skip unless LAMA_KNOWN_BUGS=true, remove the call once the bug is fixed
skip_known_bug <- function(bug) {
  skip_if(Sys.getenv("LAMA_KNOWN_BUGS") != "true", paste("known bug:", bug))
}

# ---- small test data: 3 states, 30 observations ----

Ni_test <- c(7, 5, 4)
state_test <- rep(1:3, Ni_test)
omega_test <- matrix(c(0, 0.7, 0.3,
                       0.2, 0, 0.8,
                       0.5, 0.5, 0), 3, byrow = TRUE) # from the paper
set.seed(2)
logP_hsmm <- matrix(log(runif(90)), 30, 3)
z_test <- rnorm(30)

# shifted Poisson dwell times, mean exp(beta0 + beta1 * z), one row per element of z
dm_poisson <- function(beta0, beta1, z, Ni) {
  lapply(seq_along(Ni), function(i) {
    lambda <- exp(beta0[i] + beta1[i] * z)
    d <- lapply(seq_len(Ni[i]) - 1, function(r) dpois(r, lambda))
    do.call(cbind, d) # one column per dwell time
  })
}
beta0_test <- log(c(3, 2, 1.5))
beta1_test <- c(0.4, -0.3, 0.2)
dm_hom <- lapply(dm_poisson(beta0_test, 0 * beta1_test, 0, Ni_test), as.vector)
dm_inhom <- dm_poisson(beta0_test, beta1_test, z_test, Ni_test)

# periodic: cycle length 6 is shorter than the largest aggregate
L_test <- 6
tod_test <- rep(1:L_test, length.out = 30)
dm_per <- dm_poisson(beta0_test, beta1_test, sin(2 * pi * 1:L_test / L_test), Ni_test)

# time-varying omega
omega_t_test <- array(sapply(z_test, function(z) tpm_emb(c(0.3, -0.2, 0.1) * z)), dim = c(3, 3, 30))
omega_p_test <- omega_t_test[, , 1:L_test]
