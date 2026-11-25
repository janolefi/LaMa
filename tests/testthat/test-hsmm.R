lp <- logP_hsmm
ap <- exp(lp)
lp_M <- lp[, state_test] # state-dependent densities on the extended state space

# the package adds eps = 1e-10 to the hazards, so comparisons are not exact
tol <- 1e-7

expect_close <- function(object, expected, tolerance = tol) {
  expect_equal(as.numeric(object), as.numeric(expected), tolerance = tolerance)
}

# ---- approximating transition matrix ----

test_that("tpm_hsmm matches eqs. (1)-(3)", {
  for(Ni in list(c(7, 5, 4), c(1, 5, 2), c(1, 1, 1))) {
    dm <- lapply(dm_poisson(beta0_test, 0 * beta1_test, 0, Ni), as.vector)
    Gamma <- ref_Gamma(omega_test, lapply(dm, ref_hazard))
    expect_close(as.matrix(tpm_hsmm(omega_test, dm)), Gamma)
    expect_close(tpm_hsmm(omega_test, dm, sparse = FALSE), Gamma)
  }
})

test_that("tpm_ihsmm matches eqs. (5)-(7)", {
  first <- max(Ni_test) # first time point with a transition matrix
  for(omega in list(omega_test, omega_t_test)) {
    Gamma <- tpm_ihsmm(omega, dm_inhom)
    expect_length(Gamma, 30 - first + 1)
    for(t in first:30) {
      om <- if(is.matrix(omega)) omega else omega[, , t]
      expect_close(as.matrix(Gamma[[t - first + 1]]), ref_Gamma_t(t, om, dm_inhom))
    }
  }
})

test_that("tpm_phsmm matches eqs. (5)-(7) and (9)-(10)", {
  for(omega in list(omega_test, omega_p_test)) {
    Gamma <- tpm_phsmm(omega, dm_per)
    expect_length(Gamma, L_test)
    for(t in 1:L_test) {
      om <- if(is.matrix(omega)) omega else omega[, , t]
      expect_close(as.matrix(Gamma[[t]]), ref_Gamma_t(t, om, dm_per, L_test))
    }
  }
})

test_that("tpm_hsmm, tpm_ihsmm and tpm_phsmm: arguments and checks", {
  # sparse by default
  expect_s4_class(tpm_hsmm(omega_test, dm_hom), "sparseMatrix")
  expect_true(is.matrix(tpm_hsmm(omega_test, dm_hom, sparse = FALSE)))
  expect_s4_class(tpm_ihsmm(omega_test, dm_inhom)[[1]], "sparseMatrix")
  expect_s4_class(tpm_phsmm(omega_test, dm_per)[[1]], "sparseMatrix")
  # given distribution functions Fm are used: zero means that hazards equal the dwell-time probabilities
  Fm <- lapply(dm_hom, function(d) 0 * d)
  expect_close(tpm_hsmm(omega_test, dm_hom, Fm = Fm, sparse = FALSE), ref_Gamma(omega_test, dm_hom))
  # omega needs a zero diagonal
  wrong <- omega_test; wrong[2, 2] <- 0.1
  expect_error(tpm_hsmm(wrong, dm_hom), "diagonal")
  expect_error(tpm_ihsmm(wrong, dm_inhom), "diagonal")
  expect_error(tpm_phsmm(wrong, dm_per), "diagonal")
})

# the chain that tpm_ihsmm builds has the dwell times and conditional transition
# probabilities it is supposed to have (propositions 2.1 and 2.2)
test_that("tpm_ihsmm gives the right dwell times and conditional transition probabilities", {
  first <- max(Ni_test)
  Gamma <- tpm_ihsmm(omega_t_test, dm_inhom)
  for(i in 1:3) {
    I_low <- sum(Ni_test[seq_len(i - 1)]) + 1
    entry <- first + 2 # time point at which state i is entered
    p <- replace(numeric(sum(Ni_test)), I_low, 1) # probability to still be in this stay
    for(r in seq_len(Ni_test[i])) {
      t <- entry + r - 1
      p_next <- as.vector(p %*% as.matrix(Gamma[[t - first + 1]]))
      arrive <- as.vector(rowsum(p_next, state_test))
      # dwell time is r: the state is left after r time points
      expect_close(sum(arrive[-i]), dm_inhom[[i]][entry, r])
      # next state given that i is left
      expect_close(arrive[-i] / sum(arrive[-i]), omega_t_test[i, -i, t])
      p <- p_next * (state_test == i)
    }
  }
})

# ---- direct forward step ----

test_that("hsmm_step equals phi %*% Gamma, homogeneous", {
  set.seed(3)
  for(Ni in list(c(7, 5, 4), c(1, 5, 2), c(1, 1, 1), c(2, 2), c(1, 6))) {
    N <- length(Ni)
    omega <- if(N == 3) omega_test else matrix(c(0, 1, 1, 0), 2)
    dm <- dm_poisson(beta0_test[1:N], 0 * beta1_test[1:N], 0, Ni)
    Gamma <- ref_Gamma(omega, lapply(dm, function(d) ref_hazard(as.vector(d))))
    stay_t <- unlist(lapply(dm, hsmm_stay))
    phi <- runif(sum(Ni))
    expect_close(hsmm_step(phi, stay_t, omega, Ni), phi %*% Gamma)
  }
})

test_that("hsmm_step equals phi %*% Gamma, inhomogeneous and periodic", {
  set.seed(4)
  stay <- lapply(dm_inhom, hsmm_stay)
  for(t in max(Ni_test):30) {
    phi <- runif(sum(Ni_test))
    expect_close(hsmm_step(phi, hsmm_lag(stay, t + 1)[1, ], omega_t_test[, , t], Ni_test),
                 phi %*% ref_Gamma_t(t, omega_t_test[, , t], dm_inhom))
  }
  stay <- lapply(dm_per, hsmm_stay)
  for(t in 1:L_test) {
    phi <- runif(sum(Ni_test))
    expect_close(hsmm_step(phi, hsmm_lag(stay, t + 1, L_test)[1, ], omega_test, Ni_test),
                 phi %*% ref_Gamma_t(t, omega_test, dm_per, L_test))
  }
})

test_that("hsmm_stay accepts vectors and matrices", {
  d <- dm_hom[[1]]
  stay <- 1 - ref_hazard(d)
  expect_close(hsmm_stay(d), stay)
  expect_equal(dim(hsmm_stay(d)), c(1, length(d)))
  expect_close(hsmm_stay(rbind(d, d))[2, ], stay)
})

test_that("hsmm_Gamma matches the paper", {
  for(Ni in list(c(7, 5, 4), c(1, 5, 2), c(1, 1, 1), c(2, 2), c(1, 6))) {
    N <- length(Ni)
    omega <- if(N == 3) omega_test else matrix(c(0, 1, 1, 0), 2)
    dm <- dm_poisson(beta0_test[1:N], 0 * beta1_test[1:N], 0, Ni)
    stay <- unlist(lapply(dm, hsmm_stay))
    Gamma <- hsmm_Gamma(stay, omega, Ni)
    expect_equal(dim(Gamma), rep(sum(Ni), 2))
    expect_close(Gamma, ref_Gamma(omega, lapply(dm, function(d) ref_hazard(as.vector(d)))))
    expect_close(rowSums(Gamma), rep(1, sum(Ni)), tolerance = 1e-12)
    # both give the same step
    phi <- runif(sum(Ni))
    expect_close(hsmm_step(phi, stay, omega, Ni), phi %*% Gamma, tolerance = 1e-12)
  }
})

test_that("check_omega stops for non-zero diagonals", {
  expect_silent(check_omega(omega_test))
  expect_silent(check_omega(omega_t_test))
  expect_silent(check_omega(matrix(c(0, 1, 1, 0), 2)))
  wrong <- omega_test; wrong[2, 2] <- 0.1
  expect_error(check_omega(wrong), "diagonal")
  wrong[2, 2] <- NA
  expect_error(check_omega(wrong), "diagonal")
  wrong <- omega_t_test; wrong[3, 3, 17] <- 0.1 # only one time point
  expect_error(check_omega(wrong), "diagonal")
  # while taping: constant omega is checked, AD omega is not
  wrong <- omega_test; wrong[2, 2] <- 0.1
  expect_error(RTMB::MakeTape(function(x) { check_omega(wrong); x }, 1), "diagonal")
  expect_silent(RTMB::MakeTape(function(x) { check_omega(tpm_emb(x)); sum(x) }, c(0.3, -0.2, 0.1)))
})

test_that("hsmm_step and hsmm_Gamma return plain numbers without AD", {
  stay <- unlist(lapply(dm_hom, hsmm_stay))
  expect_identical(class(hsmm_step(rep(1, 16), stay, omega_test, Ni_test)), "numeric")
  expect_true(is.numeric(hsmm_Gamma(stay, omega_test, Ni_test)))
})

# every combination of AD and constant inputs, derivatives against central differences
test_that("hsmm_step and hsmm_Gamma with AD and constant inputs", {
  phi <- 1:16 / 136
  stay_const <- unlist(lapply(dm_hom, hsmm_stay)) # not inside f: unlist() drops the AD class
  f <- function(par, ad_phi, ad_stay, ad_omega) {
    phi_par <- if(ad_phi) phi * exp(par[1]) else phi
    stay <- if(ad_stay) stay_const * exp(-exp(par[2])) else stay_const
    omega <- if(ad_omega) tpm_emb(par[3:5]) else omega_test
    weight <- cos(1:16)
    c(sum(weight * hsmm_step(phi_par, stay, omega, Ni_test)),
      sum(as.vector(phi_par) * as.vector(hsmm_Gamma(stay, omega, Ni_test) %*% weight)))
  }
  par <- c(0.1, -1, 0.3, -0.2, 0.1)
  h <- 1e-5
  combis <- expand.grid(ad_phi = c(TRUE, FALSE), ad_stay = c(TRUE, FALSE), ad_omega = c(TRUE, FALSE))
  for(k in seq_len(nrow(combis))) {
    f_k <- function(par) with(combis[k, ], f(par, ad_phi, ad_stay, ad_omega))
    if(!any(combis[k, ])) next # nothing to differentiate
    # taped at other values than evaluated: no constants are baked in
    F <- RTMB::MakeTape(function(par) { "c" <- ADoverload("c"); f_k(par) }, 0 * par)
    expect_close(F(par), f_k(par), tolerance = 1e-12)
    expect_close(F(par)[1], F(par)[2], tolerance = 1e-12)
    num <- sapply(seq_along(par), function(j) {
      e <- replace(0 * par, j, h)
      (f_k(par + e) - f_k(par - e)) / (2 * h)
    })
    expect_close(F$jacobian(par), num, tolerance = 1e-6)
  }
})

# ---- prepare_hsmm ----

# compare the transition matrices of prepare_hsmm() with Gamma_into(t) (step from t - 1 into t)
# and the initial distributions of all tracks with the rows of delta
expect_prepared <- function(inp, Gamma_into, delta) {
  Ni <- inp$Ni
  for(i in seq_along(inp$start)) {
    for(t in (inp$start[i]:inp$end[i])[-1]) {
      expect_close(hsmm_Gamma(inp$stay[t, ], inp$omega[, , inp$oidx[t]], Ni), Gamma_into(t, i))
    }
  }
  expect_equal(dim(inp$delta), c(length(inp$start), sum(Ni)))
  expect_close(inp$delta, delta, tolerance = 1e-6)
}

two_tracks <- rep(1:2, c(12, 18))
delta_N <- c(0.5, 0.3, 0.2)
delta_M <- ref_delta_fresh(delta_N, Ni_test)

test_that("prepare_hsmm: homogeneous", {
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  inp <- prepare_hsmm(dm_hom, omega_test, ap)
  expect_prepared(inp, function(t, i) Gamma, t(ref_delta_stat(Gamma)))
  expect_identical(inp$type, "homogeneous")
  expect_equal(c(inp$start, inp$end), c(1, 30))
  # tracks and given delta
  inp <- prepare_hsmm(dm_hom, omega_test, ap, two_tracks)
  expect_prepared(inp, function(t, i) Gamma, rbind(ref_delta_stat(Gamma), ref_delta_stat(Gamma)))
  expect_equal(inp$start, c(1, 13)); expect_equal(inp$end, c(12, 30))
  expect_close(prepare_hsmm(dm_hom, omega_test, ap, two_tracks, delta_N)$delta, rbind(delta_M, delta_M))
  Delta <- rbind(delta_N, rev(delta_N))
  expect_close(prepare_hsmm(dm_hom, omega_test, ap, two_tracks, Delta)$delta,
               rbind(delta_M, ref_delta_fresh(rev(delta_N), Ni_test)))
  # column vector, single observation
  expect_close(prepare_hsmm(dm_hom, omega_test, ap, delta = matrix(delta_N))$delta, delta_M)
  expect_close(prepare_hsmm(dm_hom, omega_test, ap[1, , drop = FALSE])$delta, ref_delta_stat(Gamma), tolerance = 1e-6)
})

test_that("prepare_hsmm: track-specific dm and omega", {
  dm2 <- dm_poisson(beta0_test, beta1_test, c(-0.5, 0.8), Ni_test)
  omega2 <- omega_t_test[, , 1:2]
  Gamma <- lapply(1:2, function(i) ref_Gamma(omega2[, , i], lapply(dm2, function(d) ref_hazard(d[i, ]))))
  inp <- prepare_hsmm(dm2, omega2, ap, two_tracks)
  expect_prepared(inp, function(t, i) Gamma[[i]], t(sapply(Gamma, ref_delta_stat)))
  expect_identical(inp$type, "track-specific")
})

test_that("prepare_hsmm: time-varying dm, omega[, , t] for the step into t", {
  first <- max(Ni_test)
  for(omega in list(omega_test, omega_t_test)) {
    om <- function(t) if(is.matrix(omega)) omega else omega[, , t]
    Gamma_into <- function(t, i) ref_Gamma_t(t - 1, om(t), dm_inhom)
    inp <- prepare_hsmm(dm_inhom, omega, ap)
    expect_equal(inp$start, first)
    expect_prepared(inp, Gamma_into, t(ref_delta_stat(Gamma_into(first + 1))))
    inp <- prepare_hsmm(dm_inhom, omega, ap, first_obs = 10)
    expect_equal(inp$start, 10)
    expect_prepared(inp, Gamma_into, t(ref_delta_stat(Gamma_into(11))))
    # given delta
    expect_prepared(prepare_hsmm(dm_inhom, omega, ap, delta = delta_N), Gamma_into, t(delta_M))
    # tracks do not use rows of earlier tracks
    starts <- c(1, 13)
    Gamma_into <- function(t, i) {
      dm_i <- lapply(dm_inhom, function(d) d[starts[i]:nrow(d), ])
      ref_Gamma_t(t - starts[i], om(t), dm_i)
    }
    inp <- prepare_hsmm(dm_inhom, omega, ap, two_tracks)
    expect_equal(inp$start, starts + first - 1)
    expect_prepared(inp, Gamma_into, t(sapply(1:2, function(i) ref_delta_stat(Gamma_into(starts[i] + first, i)))))
  }
  expect_identical(inp$type, "inhomogeneous")
})

test_that("prepare_hsmm: tracks do not depend on each other", {
  # changing dm of track 1 does not change track 2
  dm2 <- lapply(dm_inhom, function(d) { d[1:12, ] <- 0.1; d })
  a <- prepare_hsmm(dm_inhom, omega_test, ap, two_tracks)
  b <- prepare_hsmm(dm2, omega_test, ap, two_tracks)
  rows <- b$start[2]:30
  expect_identical(a$stay[rows, ], b$stay[rows, ])
  expect_identical(a$delta[2, ], b$delta[2, ])
})

test_that("prepare_hsmm: fixed dm with time-varying omega", {
  hazard <- lapply(dm_hom, ref_hazard)
  Gamma_into <- function(t, i) ref_Gamma(omega_t_test[, , t], hazard)
  inp <- prepare_hsmm(dm_hom, omega_t_test, ap)
  expect_equal(inp$start, 1)
  expect_prepared(inp, Gamma_into, t(ref_delta_stat(Gamma_into(2))))
  # without the slice into t = 1
  inp2 <- prepare_hsmm(dm_hom, omega_t_test[, , -1], ap)
  expect_prepared(inp2, Gamma_into, inp$delta)
})

test_that("prepare_hsmm: fixed dwell times in one state, varying ones in the others", {
  dm_mix <- dm_inhom
  dm_mix[[2]] <- dm_hom[[2]]
  dm_full <- dm_mix
  dm_full[[2]] <- matrix(dm_hom[[2]], 30, Ni_test[2], byrow = TRUE)
  a <- prepare_hsmm(dm_mix, omega_test, ap)
  b <- prepare_hsmm(dm_full, omega_test, ap)
  expect_equal(a[c("stay", "delta", "start")], b[c("stay", "delta", "start")])
})

test_that("prepare_hsmm: periodic", {
  tod <- c(3:6, rep(1:6, 4), 1:2) # does not start at 1
  for(case in c("dm", "omega", "both")) {
    dm <- if(case == "omega") dm_hom else dm_per
    omega <- if(case == "dm") omega_test else omega_p_test
    dm_ref <- if(case == "omega") lapply(dm_hom, function(d) matrix(d, 6, length(d), byrow = TRUE)) else dm_per
    om <- function(s) if(is.matrix(omega)) omega else omega[, , s]
    into <- function(s) ref_Gamma_t(s - 1, om(s), dm_ref, L_test) # step into time of day s
    inp <- prepare_hsmm(dm, omega, ap, tod = tod)
    # stationary at time of day 3: cycle 4, 5, 6, 1, 2, 3
    Gamma_cycle <- Reduce(`%*%`, lapply(c(4:6, 1:3), into))
    expect_prepared(inp, function(t, i) into(tod[t]), t(ref_delta_stat(Gamma_cycle)))
    expect_identical(inp$type, "periodic")
    # track-specific start, not all times of day observed
    tod2 <- c(5, 6, 1:4)
    inp <- prepare_hsmm(dm, omega, ap[1:6, ], c(1, 1, 2, 2, 2, 2), tod = tod2)
    Delta <- rbind(ref_delta_stat(Reduce(`%*%`, lapply(c(6, 1:5), into))),
                   ref_delta_stat(Reduce(`%*%`, lapply(c(2:6, 1), into))))
    expect_prepared(inp, function(t, i) into(tod2[t]), Delta)
  }
})

test_that("prepare_hsmm: use_all_obs", {
  inp <- prepare_hsmm(dm_inhom, omega_test, ap, two_tracks, delta_N, use_all_obs = TRUE)
  expect_equal(inp$start, c(1, 13))
  expect_error(prepare_hsmm(dm_inhom, omega_test, ap, use_all_obs = TRUE), "needs delta")
  expect_error(prepare_hsmm(dm_inhom, omega_test, ap, delta = "stationary", use_all_obs = TRUE), "needs delta")
  # no effect for fixed and periodic dm
  expect_identical(prepare_hsmm(dm_hom, omega_test, ap, use_all_obs = TRUE), prepare_hsmm(dm_hom, omega_test, ap))
  expect_identical(prepare_hsmm(dm_per, omega_test, ap, tod = tod_test, use_all_obs = TRUE),
                   prepare_hsmm(dm_per, omega_test, ap, tod = tod_test))
})

test_that("prepare_hsmm: checks", {
  prep <- function(dm = dm_hom, omega = omega_test, allprobs = ap, ...) prepare_hsmm(dm, omega, allprobs, ...)
  expect_error(prep(allprobs = ap[, 1]), "matrix")
  expect_error(prep(dm = dm_hom[1:2]), "each state")
  expect_error(prep(dm = dm_hom[[1]]), "each state")
  expect_error(prep(omega = omega_test + diag(3)), "diagonal")
  expect_error(prep(omega = matrix(c(0, 1, 1, 0), 2)), "square")
  expect_error(prep(omega = omega_t_test[, , 1:5]), "slices of omega")
  expect_error(prep(dm = lapply(dm_inhom, function(d) d[1:5, ])), "rows of dm")
  expect_error(prep(dm = lapply(dm_inhom, function(d) d[-1, ])), "rows of dm") # nObs - 1 only for omega
  expect_error(prep(trackID = rep(1:2, 15)), "contiguous")
  expect_error(prep(trackID = 1:3), "trackID")
  expect_error(prep(delta = c(0.5, 0.5)), "delta")
  expect_error(prep(delta = "stat"), "delta")
  expect_error(prep(trackID = two_tracks, delta = rbind(delta_N, delta_N, delta_N)), "delta")
  # time-varying dm
  expect_error(prep(dm = dm_inhom, first_obs = 6), "at least 7")
  for(wrong in list(0, 1.5, NA, 1:2)) expect_error(prep(first_obs = wrong), "first_obs")
  expect_error(prep(first_obs = 31), "Too few")
  expect_error(prep(trackID = two_tracks, first_obs = 13), "Too few")
  expect_silent(prep(trackID = two_tracks, first_obs = 12))
  expect_error(prep(delta = delta_N, use_all_obs = TRUE, first_obs = 3), "either")
  expect_error(prep(dm = dm_inhom, trackID = rep(1:2, c(7, 23))), "Too few")
  expect_silent(prep(dm = dm_inhom, trackID = rep(1:2, c(8, 22))))
  expect_silent(prep(dm = dm_inhom, trackID = rep(1:2, c(7, 23)), delta = delta_N))
  expect_error(prep(dm = dm_inhom, trackID = rep(1:2, c(6, 24)), delta = delta_N), "Too few")
  # periodic
  expect_error(prep(dm = dm_per, tod = 1:3), "same length")
  expect_error(prep(dm = dm_per, tod = rep(1:10, 3)), "rows of dm")
  expect_error(prep(dm = dm_per, omega = omega_t_test[, , 1:8], tod = tod_test), "cycle length")
  expect_error(prep(dm = dm_per, tod = tod_test[c(1:10, 12:30, 1)]), "advance by one")
  expect_error(prep(dm = dm_per, tod = rev(tod_test)), "advance by one")
  # a new track may start at any time of day, and the rule only holds for dm
  expect_silent(prep(dm = dm_per, tod = c(1:6, 1:6, rep(4, 1), 5, 6, 1:6, 1:6, 1:3), trackID = rep(1:2, c(12, 18))))
  expect_silent(prep(omega = omega_p_test, tod = rep(c(2, 5), 15)))
})

test_that("prepare_hsmm works with AD, also if only some inputs are AD", {
  par <- c(beta0_test, beta1_test, 0.3, -0.2, 0.1, log(c(0.6, 0.4)))
  # sum of everything that depends on parameters
  total <- function(inp) sum(inp$stay[inp$start[1]:30, ]) + sum(inp$delta * cos(1:length(inp$delta))) + sum(inp$omega * 2)
  f <- function(par, ad_dm, ad_omega, ad_delta) {
    "c" <- ADoverload("c")
    dm <- if(ad_dm) dm_poisson(par[1:3], par[4:6], z_test, Ni_test) else dm_inhom
    if(ad_dm) dm[[2]] <- dm_hom[[2]] # constant in one state
    omega <- if(ad_omega) tpm_emb(par[7:9]) else omega_test
    delta <- if(ad_delta) c(exp(par[10:11]), 1) / (1 + sum(exp(par[10:11])))
    total(prepare_hsmm(dm, omega, ap, delta = delta))
  }
  combis <- expand.grid(ad_dm = c(TRUE, FALSE), ad_omega = c(TRUE, FALSE), ad_delta = c(TRUE, FALSE))
  h <- 1e-5
  for(k in seq_len(nrow(combis) - 1)) { # last one: nothing is AD
    f_k <- function(par) f(par, combis$ad_dm[k], combis$ad_omega[k], combis$ad_delta[k])
    F <- RTMB::MakeTape(f_k, 0.9 * par)
    expect_close(F(par), f_k(par), tolerance = 1e-10)
    num <- sapply(seq_along(par), function(j) {
      e <- replace(0 * par, j, h)
      (f_k(par + e) - f_k(par - e)) / (2 * h)
    })
    expect_close(F$jacobian(par), num, tolerance = 1e-5)
  }
})

# ---- forward_hsmm_ad ----

# compare log-likelihood with ref, with and without AD, with and without logspace
expect_loglik <- function(ref, ..., logP = lp) {
  for(logspace in c(FALSE, TRUE)) {
    f <- function(x) new_forward(..., allprobs = if(logspace) x else exp(x), logspace = logspace)
    F <- RTMB::MakeTape(function(x) f(matrix(x, nrow(logP))), 0.9 * as.vector(logP))
    expect_close(F(as.vector(logP)), ref)
    expect_close(f(logP), ref)
  }
}

# log-likelihood of the reference: Gamma_into(t, i), delta (nTracks x M), used observations
ref_loglik <- function(Gamma_into, delta, used = 1:30, ID = rep(1, 30), logP = lp) {
  ref_forward(delta, function(t, i) Gamma_into(used[t], i), logP[used, state_test, drop = FALSE], ID[used])
}

test_that("forward_hsmm_ad: homogeneous", {
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  delta <- t(ref_delta_stat(Gamma))
  into <- function(t, i) Gamma
  expect_loglik(ref_loglik(into, delta), dm = dm_hom, omega = omega_test)
  expect_loglik(ref_loglik(into, t(delta_M)), dm = dm_hom, omega = omega_test, delta = delta_N)
  # tracks, also of length one
  for(ID in list(two_tracks, rep(1:3, c(8, 10, 12)), rep(c("b", "a", "c"), c(1, 28, 1)))) {
    Delta <- delta[rep(1, length(unique(ID))), ]
    expect_loglik(ref_loglik(into, Delta, ID = ID), dm = dm_hom, omega = omega_test, trackID = ID)
  }
  expect_loglik(ref_loglik(into, delta, used = 1), dm = dm_hom, omega = omega_test, logP = lp[1, , drop = FALSE])
})

test_that("forward_hsmm_ad: aggregates of size 1 and 2", {
  for(Ni in list(c(1, 5, 2), c(1, 1, 1))) {
    dm <- lapply(dm_poisson(beta0_test, 0 * beta1_test, 0, Ni), as.vector)
    Gamma <- ref_Gamma(omega_test, lapply(dm, ref_hazard))
    ref <- ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp[, rep(1:3, Ni)])
    expect_loglik(ref, dm = dm, omega = omega_test)
  }
})

test_that("forward_hsmm_ad: track-specific dm, omega and delta", {
  dm2 <- dm_poisson(beta0_test, beta1_test, c(-0.5, 0.8), Ni_test)
  omega2 <- omega_t_test[, , 1:2]
  Gamma <- lapply(1:2, function(i) ref_Gamma(omega2[, , i], lapply(dm2, function(d) ref_hazard(d[i, ]))))
  into <- function(t, i) Gamma[[i]]
  expect_loglik(ref_loglik(into, t(sapply(Gamma, ref_delta_stat)), ID = two_tracks),
                dm = dm2, omega = omega2, trackID = two_tracks)
  Delta <- rbind(delta_N, rev(delta_N))
  expect_loglik(ref_loglik(into, t(apply(Delta, 1, ref_delta_fresh, Ni = Ni_test)), ID = two_tracks),
                dm = dm2, omega = omega2, trackID = two_tracks, delta = Delta)
})

test_that("forward_hsmm_ad: time-varying dm and omega", {
  first <- max(Ni_test)
  for(omega in list(omega_test, omega_t_test)) {
    om <- function(t) if(is.matrix(omega)) omega else omega[, , t]
    into <- function(t, i) ref_Gamma_t(t - 1, om(t), dm_inhom)
    expect_loglik(ref_loglik(into, t(ref_delta_stat(into(first + 1))), first:30), dm = dm_inhom, omega = omega)
    expect_loglik(ref_loglik(into, t(ref_delta_stat(into(11))), 10:30), dm = dm_inhom, omega = omega, first_obs = 10)
    expect_loglik(ref_loglik(into, t(delta_M), first:30), dm = dm_inhom, omega = omega, delta = delta_N)
    # two tracks
    starts <- c(1, 13)
    into <- function(t, i) ref_Gamma_t(t - starts[i], om(t), lapply(dm_inhom, function(d) d[starts[i]:30, ]))
    Delta <- t(sapply(1:2, function(i) ref_delta_stat(into(starts[i] + first, i))))
    used <- c(first:12, (12 + first):30)
    expect_loglik(ref_loglik(into, Delta, used, two_tracks), dm = dm_inhom, omega = omega, trackID = two_tracks)
  }
  # fixed dm
  into <- function(t, i) ref_Gamma(omega_t_test[, , t], lapply(dm_hom, ref_hazard))
  expect_loglik(ref_loglik(into, t(ref_delta_stat(into(2)))), dm = dm_hom, omega = omega_t_test)
})

test_that("forward_hsmm_ad: first_obs for fixed, track-specific and periodic models", {
  # homogeneous: same as without the first observations
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  into <- function(t, i) Gamma
  delta <- t(ref_delta_stat(Gamma))
  expect_loglik(ref_loglik(into, delta, 9:30), dm = dm_hom, omega = omega_test, first_obs = 9)
  expect_loglik(ref_loglik(into, t(delta_M), 9:30), dm = dm_hom, omega = omega_test, first_obs = 9, delta = delta_N)
  expect_loglik(ref_loglik(into, delta[c(1, 1), ], c(5:12, 17:30), two_tracks),
                dm = dm_hom, omega = omega_test, trackID = two_tracks, first_obs = 5)
  # time-varying omega: stationary for the step out of the first observation used
  into <- function(t, i) ref_Gamma(omega_t_test[, , t], lapply(dm_hom, ref_hazard))
  expect_loglik(ref_loglik(into, t(ref_delta_stat(into(10))), 9:30), dm = dm_hom, omega = omega_t_test, first_obs = 9)
  # periodic: stationary for its time of day
  into_tod <- function(s) ref_Gamma_t(s - 1, omega_p_test[, , s], dm_per, L_test)
  delta <- t(ref_delta_stat(Reduce(`%*%`, lapply(c(4:6, 1:3), into_tod)))) # observation 9: time of day 3
  expect_loglik(ref_loglik(function(t, i) into_tod(tod_test[t]), delta, 9:30),
                dm = dm_per, omega = omega_p_test, tod = tod_test, first_obs = 9)
})

test_that("forward_hsmm_ad: use_all_obs is exact for a given delta", {
  starts <- c(1, 13)
  for(value in c(0.01, 0.1)) { # nonsense dm before the first observation of each track has no effect
    into <- function(t, i) {
      dm_i <- lapply(dm_inhom, function(d) rbind(matrix(value, 10, ncol(d)), d[starts[i]:30, ]))
      ref_Gamma_t(t - starts[i] + 10, omega_t_test[, , t], dm_i)
    }
    expect_loglik(ref_loglik(into, rbind(delta_M, delta_M), ID = two_tracks), dm = dm_inhom, omega = omega_t_test,
                  trackID = two_tracks, delta = delta_N, use_all_obs = TRUE)
  }
})

test_that("forward_hsmm_ad: periodic dm and omega", {
  tod <- c(3:6, rep(1:6, 4), 1:2)
  into_tod <- function(s) ref_Gamma_t(s - 1, omega_p_test[, , s], dm_per, L_test)
  into <- function(t, i) into_tod(tod[t])
  delta <- t(ref_delta_stat(Reduce(`%*%`, lapply(c(4:6, 1:3), into_tod))))
  expect_loglik(ref_loglik(into, delta), dm = dm_per, omega = omega_p_test, tod = tod)
  expect_loglik(ref_loglik(into, t(delta_M)), dm = dm_per, omega = omega_p_test, tod = tod, delta = delta_N)
})

test_that("forward_hsmm_ad agrees with forward() on the extended chain, value and gradient", {
  par <- c(beta0_test, beta1_test, 0.3, -0.2, 0.1)
  f_new <- function(par) {
    new_forward(dm_poisson(par[1:3], par[4:6], z_test, Ni_test), tpm_emb(par[7:9]), ap, trackID = two_tracks)
  }
  f_dense <- function(par) {
    inp <- prepare_hsmm(dm_poisson(par[1:3], par[4:6], z_test, Ni_test), tpm_emb(par[7:9]), ap, two_tracks)
    used <- unlist(Map(seq, inp$start, inp$end))
    Gamma <- hsmm_Gammas(inp$stay[used, ], inp$omega[, , inp$oidx[used]], Ni_test)
    forward(inp$delta, Gamma, ap[used, state_test], trackID = two_tracks[used], report = FALSE)
  }
  F_new <- RTMB::MakeTape(f_new, 0.9 * par)
  F_dense <- RTMB::MakeTape(f_dense, 0.9 * par)
  expect_close(F_new(par), F_dense(par), tolerance = 1e-12)
  expect_close(F_new(par), f_new(par), tolerance = 1e-12)
  expect_close(F_new$jacobian(par), F_dense$jacobian(par), tolerance = 1e-10)
  h <- 1e-5
  num <- sapply(seq_along(par), function(j) {
    e <- replace(0 * par, j, h)
    (f_new(par + e) - f_new(par - e)) / (2 * h)
  })
  expect_close(F_new$jacobian(par), num, tolerance = 1e-6)
})

test_that("forward_hsmm_ad: delta as a parameter", {
  F <- RTMB::MakeTape(function(delta) new_forward(dm_hom, omega_test, ap, delta = delta), rev(delta_N))
  f <- function(delta) new_forward(dm_hom, omega_test, ap, delta = delta)
  expect_close(F(delta_N), f(delta_N), tolerance = 1e-12)
  num <- sapply(1:3, function(j) (f(delta_N + 1e-6 * (1:3 == j)) - f(delta_N - 1e-6 * (1:3 == j))) / 2e-6)
  expect_close(F$jacobian(delta_N), num, tolerance = 1e-6)
})

test_that("forward_hsmm_ad: weights of the log-likelihood contributions", {
  first <- max(Ni_test)
  used <- c(first:12, (12 + first):30)
  starts <- c(1, 13)
  into <- function(t, i) ref_Gamma_t(t - starts[i], omega_test, lapply(dm_inhom, function(d) d[starts[i]:30, ]))
  Delta <- t(sapply(1:2, function(i) ref_delta_stat(into(starts[i] + first, i))))
  contrib <- ref_contrib(Delta, function(t, i) into(used[t], i), lp[used, state_test], two_tracks[used])
  set.seed(5)
  w <- runif(length(used))
  for(logspace in c(FALSE, TRUE)) {
    l <- forward_weighted(function(x) new_forward(dm_inhom, omega_test, x, two_tracks, logspace = logspace), lp, w, logspace)
    expect_close(l, sum(w * contrib))
  }
})

test_that("forward_hsmm_ad: dropped observations and missing values", {
  dropped <- 1:(max(Ni_test) - 1)
  changed <- lp; changed[dropped, ] <- -5
  for(logspace in c(FALSE, TRUE)) {
    trans <- if(logspace) identity else exp
    expect_close(new_forward(dm_inhom, omega_test, trans(changed), logspace = logspace),
                 new_forward(dm_inhom, omega_test, trans(lp), logspace = logspace), tolerance = 1e-12)
  }
  missing <- ap; missing[3, 1] <- NA
  expect_error(new_forward(dm_hom, omega_test, missing), "NA")
  expect_silent(new_forward(dm_inhom, omega_test, missing)) # row 3 is dropped
  expect_error(new_forward(dm_hom, omega_test, lp), "logspace")
})

# ---- current forward functions ----

# log-likelihood on the extended state space, Gamma_t(t): Gamma^(t), used for the step t -> t + 1
ref_hsmm <- function(delta, Gamma_t, rows = seq_len(nrow(lp_M))) {
  ref_forward(delta, function(t, i) Gamma_t(rows[t] - 1), lp_M[rows, , drop = FALSE])
}

# value with and without AD
expect_hsmm <- function(f, ref) {
  F <- RTMB::MakeTape(function(x) f(exp(matrix(x, nrow(lp)))), as.vector(lp))
  expect_close(F(as.vector(lp)), ref)
  expect_close(f(ap), ref)
}

test_that("forward_hsmm matches the paper", {
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  delta <- c(0.5, 0.3, 0.2)
  expect_hsmm(function(ap) forward_hsmm(dm_hom, omega_test, ap, report = FALSE),
              ref_hsmm(ref_delta_stat(Gamma), function(t) Gamma))
  expect_hsmm(function(ap) forward_hsmm(dm_hom, omega_test, ap, delta = delta, report = FALSE),
              ref_hsmm(ref_delta_fresh(delta, Ni_test), function(t) Gamma))
})

test_that("forward_hsmm with two tracks matches the paper", {
  ID <- rep(1:2, c(12, 18))
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ref <- ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID)
  expect_close(forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE), ref)
  delta <- c(0.5, 0.3, 0.2)
  ref <- ref_forward(ref_delta_fresh(delta, Ni_test), function(t, i) Gamma, lp_M, ID)
  expect_hsmm(function(ap) forward_hsmm(dm_hom, omega_test, ap, trackID = ID, delta = delta, report = FALSE), ref)
})

test_that("forward_ihsmm matches the paper", {
  first <- max(Ni_test)
  for(omega in list(omega_test, omega_t_test)) {
    Gamma_t <- function(t) {
      om <- if(is.matrix(omega)) omega else omega[, , t]
      ref_Gamma_t(t, om, dm_inhom)
    }
    # starts at the largest aggregate size, or at startInd
    for(start in c(first, first + 4)) {
      delta <- ref_delta_stat(Gamma_t(start))
      f <- function(ap) forward_ihsmm(dm_inhom, omega, ap, startInd = start, report = FALSE)
      expect_hsmm(f, ref_hsmm(delta, Gamma_t, start:30))
    }
    expect_hsmm(function(ap) forward_ihsmm(dm_inhom, omega, ap, report = FALSE),
                ref_hsmm(ref_delta_stat(Gamma_t(first)), Gamma_t, first:30))
    # given delta
    delta <- c(0.5, 0.3, 0.2)
    expect_hsmm(function(ap) forward_ihsmm(dm_inhom, omega, ap, delta = delta, report = FALSE),
                ref_hsmm(ref_delta_fresh(delta, Ni_test), Gamma_t, first:30))
  }
  expect_error(forward_ihsmm(dm_inhom, omega_test, ap, startInd = first - 1), "startInd")
  # no effect for fixed dm, as before
  expect_error(forward_ihsmm(dm_hom, omega_t_test, ap, startInd = first - 1), "startInd")
  expect_identical(forward_ihsmm(dm_hom, omega_t_test, ap, startInd = 10, report = FALSE),
                   forward_ihsmm(dm_hom, omega_t_test, ap, report = FALSE))
})

test_that("forward_ihsmm with fixed dm and time-varying omega matches the paper", {
  Gamma_t <- function(t) ref_Gamma(omega_t_test[, , t], lapply(dm_hom, ref_hazard))
  expect_hsmm(function(ap) forward_ihsmm(dm_hom, omega_t_test, ap, report = FALSE),
              ref_hsmm(ref_delta_stat(Gamma_t(1)), Gamma_t))
})

test_that("forward_phsmm matches the paper", {
  for(omega in list(omega_test, omega_p_test)) {
    Gamma_t <- function(t) {
      om <- if(is.matrix(omega)) omega else omega[, , t]
      ref_Gamma_t(t, om, dm_per, L_test)
    }
    # periodically stationary start at tod 1
    Gamma_cycle <- Reduce(`%*%`, lapply(1:L_test, Gamma_t))
    delta <- ref_delta_stat(Gamma_cycle)
    expect_hsmm(function(ap) forward_phsmm(dm_per, omega, ap, tod_test, report = FALSE),
                ref_hsmm(delta, function(t) Gamma_t(tod_test[t])))
    delta <- c(0.5, 0.3, 0.2)
    expect_hsmm(function(ap) forward_phsmm(dm_per, omega, ap, tod_test, delta = delta, report = FALSE),
                ref_hsmm(ref_delta_fresh(delta, Ni_test), function(t) Gamma_t(tod_test[t])))
  }
})

test_that("forward_hsmm: all types of models", {
  f <- function(...) forward_hsmm(..., report = FALSE)
  expect_identical(f(dm_inhom, omega_t_test, ap), new_forward(dm_inhom, omega_t_test, ap))
  expect_identical(f(dm_inhom, omega_t_test, ap, two_tracks, delta_N, use_all_obs = TRUE),
                   new_forward(dm_inhom, omega_t_test, ap, two_tracks, delta_N, use_all_obs = TRUE))
  expect_identical(f(dm_per, omega_p_test, ap, tod = tod_test), new_forward(dm_per, omega_p_test, ap, tod = tod_test))
  expect_identical(f(dm_hom, omega_test, ap, first_obs = 4), new_forward(dm_hom, omega_test, ap, first_obs = 4))
  # models with different aggregate sizes on the same observations
  dm_large <- dm_poisson(beta0_test, beta1_test, z_test, c(11, 10, 9))
  expect_false(f(dm_inhom, omega_test, ap) == f(dm_inhom, omega_test, ap, first_obs = 11))
  expect_identical(f(dm_large, omega_test, ap), f(dm_large, omega_test, ap, first_obs = 11))
  expect_identical(f(dm_hom, omega_test, lp, logspace = TRUE), new_forward(dm_hom, omega_test, lp, logspace = TRUE))
  expect_false(f(dm_hom, omega_test, ap, eps = 1e-3) == f(dm_hom, omega_test, ap))
  # nested lists for track-specific dm
  expect_error(f(list(dm_hom, dm_hom, dm_hom), omega_test, ap, rep(1:3, each = 10)), "nested list")
})

test_that("forward_hsmm: reported quantities", {
  nll <- function(par) {
    dm <- dm_poisson(par$beta0, par$beta1, z_test, Ni_test)
    -forward_hsmm(dm, tpm_emb(par$eta), lp, trackID = two_tracks, logspace = TRUE)
  }
  par <- list(beta0 = beta0_test, beta1 = beta1_test, eta = c(0.3, -0.2, 0.1))
  obj <- RTMB::MakeADFun(nll, par, silent = TRUE)
  rep <- obj$report()
  expect_equal(dim(rep$delta), c(2, sum(Ni_test)))
  expect_close(rep$dm[[2]], dm_inhom[[2]])
  expect_close(rep$omega, tpm_emb(par$eta))
  expect_equal(rep$Ni, Ni_test)
  expect_close(rep$allprobs, ap)
  expect_equal(rep$trackID, two_tracks)
  expect_null(RTMB::MakeADFun(function(par) -forward_hsmm(dm_hom, tpm_emb(par$eta), ap, report = FALSE), par["eta"], silent = TRUE)$report()$dm)
  # model object
  mod <- lama(obj)
  expect_output(print(mod), "3 states")
  expect_output(print(summary(mod)))
  expect_close(logLik(mod), -nll(par))
})

test_that("forward_hsmm: delta_sparse and Gamma_sparse are still reported", {
  par <- list(beta0 = beta0_test, eta = c(0.3, -0.2, 0.1))
  report <- function(...) {
    nll <- function(par) {
      dm <- lapply(dm_poisson(par$beta0, 0 * par$beta0, 0, Ni_test), as.vector)
      -forward_hsmm(dm, tpm_emb(par$eta), ap, ...)
    }
    RTMB::MakeADFun(nll, par, silent = TRUE)$report()
  }
  Gamma <- ref_Gamma(tpm_emb(par$eta), lapply(dm_hom, ref_hazard))
  rep <- report()
  expect_s4_class(rep$Gamma_sparse, "sparseMatrix")
  expect_close(as.matrix(rep$Gamma_sparse), Gamma)
  expect_s4_class(rep$delta_sparse, "sparseMatrix")
  expect_close(as.matrix(rep$delta_sparse), ref_delta_stat(Gamma), tolerance = 1e-6)
  # the reported matrices can be used for decoding
  expect_length(viterbi(as.vector(rep$delta_sparse), as.matrix(rep$Gamma_sparse), rep$allprobs[, state_test]), 30)
  # dm as passed, not as matrices
  expect_null(dim(rep$dm[[1]]))
  # tracks, with their original labels
  ID <- rep(c("b", "a"), c(12, 18))
  rep <- report(trackID = ID, delta = rbind(delta_N, rev(delta_N)))
  expect_identical(rep$trackID, ID)
  expect_length(rep$Gamma_sparse, 2)
  expect_close(as.matrix(rep$Gamma_sparse[[2]]), Gamma)
  expect_close(as.matrix(rep$Delta_sparse), rbind(delta_M, ref_delta_fresh(rev(delta_N), Ni_test)))
  expect_null(report(report = FALSE)$Gamma_sparse)
})

test_that("forward_ihsmm: reported quantities, omega convention and message", {
  par <- list(beta0 = beta0_test, beta1 = beta1_test, eta = cbind(c(0.3, -0.2, 0.1), c(0.1, 0.2, 0.3)))
  nll <- function(par) {
    dm <- dm_poisson(par$beta0, par$beta1, z_test, Ni_test)
    omega <- tpm_emb_g(cbind(1, z_test), par$eta, report = FALSE)
    -forward_ihsmm(dm, omega, ap, startInd = 9)
  }
  rep <- RTMB::MakeADFun(nll, par, silent = TRUE)$report()
  omega <- tpm_emb_g(cbind(1, z_test), par$eta, report = FALSE)
  # as before: matrices from max(Ni) on, delta at startInd
  expect_equal(rep$Gamma_sparse, tpm_ihsmm(omega, dm_inhom))
  expect_close(as.matrix(rep$delta_sparse), ref_delta_stat(ref_Gamma_t(9, omega[, , 9], dm_inhom)), tolerance = 1e-6)
  expect_close(rep$omega, omega) # as passed
  expect_close(rep$dm[[3]], dm_inhom[[3]])
  # forward_hsmm() with the slices moved by one
  expect_equal(forward_ihsmm(dm_inhom, omega, ap, report = FALSE),
               forward_hsmm(dm_inhom, omega[, , c(1, 1:29)], ap, report = FALSE))
  expect_error(forward_ihsmm(dm_inhom, omega[, , 1:29], ap, report = FALSE), "number of observations")
  expect_null(RTMB::MakeADFun(function(par) -forward_ihsmm(dm_inhom, tpm_emb(par$eta), ap, report = FALSE),
                              list(eta = 1:3 / 10), silent = TRUE)$report()$Gamma_sparse)
})

test_that("report_sparse does not build the matrices with AD", {
  inp <- prepare_hsmm(dm_hom, omega_test, ap)
  expect_silent(RTMB::MakeTape(function(x) { report_sparse(inp, stop("built")); x }, 1))
  expect_error(report_sparse(inp, stop("built")), "built")
})

# ---- bugs of earlier versions ----

test_that("forward_ihsmm: fixed dm and given delta", {
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ref <- ref_hsmm(delta_M, function(t) Gamma)
  expect_hsmm(function(ap) forward_ihsmm(dm_hom, omega_test, ap, delta = delta_N, report = FALSE), ref)
  Gamma_t <- function(t) ref_Gamma(omega_t_test[, , t], lapply(dm_hom, ref_hazard))
  expect_hsmm(function(ap) forward_ihsmm(dm_hom, omega_t_test, ap, delta = delta_N, report = FALSE),
              ref_hsmm(delta_M, Gamma_t))
})

test_that("forward_phsmm: fixed dm", {
  hazard <- lapply(dm_hom, ref_hazard)
  Gamma <- ref_Gamma(omega_test, hazard)
  for(delta in list(NULL, delta_N)) {
    start <- if(is.null(delta)) ref_delta_stat(Gamma) else delta_M
    expect_hsmm(function(ap) forward_phsmm(dm_hom, omega_test, ap, tod_test, delta = delta, report = FALSE),
                ref_hsmm(start, function(t) Gamma))
  }
  # periodic omega, also if tod does not advance by one
  Gamma_s <- function(s) ref_Gamma(omega_p_test[, , s], hazard)
  for(tod in list(tod_test, rep(c(2, 5), 15), rep(c(4, 4, 1), 10))) {
    start <- ref_delta_stat(Reduce(`%*%`, lapply((tod[1] + 0:5 - 1) %% 6 + 1, Gamma_s)))
    for(delta in list(NULL, delta_N)) {
      expect_hsmm(function(ap) forward_phsmm(dm_hom, omega_p_test, ap, tod, delta = delta, report = FALSE),
                  ref_hsmm(if(is.null(delta)) start else delta_M, function(t) Gamma_s(tod[t])))
    }
  }
})

test_that("forward_hsmm: delta can be estimated", {
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  delta <- c(0.5, 0.3, 0.2)
  F <- RTMB::MakeTape(function(delta) forward_hsmm(dm_hom, omega_test, ap, delta = delta, report = FALSE), delta)
  expect_close(F(delta), ref_hsmm(ref_delta_fresh(delta, Ni_test), function(t) Gamma))
})

test_that("forward_hsmm: several tracks with stationary delta under AD", {
  ID <- rep(1:2, c(12, 18))
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ref <- ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID)
  expect_hsmm(function(ap) forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE), ref)
})

test_that("forward_hsmm: as many tracks as states", {
  ID <- rep(1:3, c(8, 10, 12))
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  expect_close(forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE),
               ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID))
})

test_that("forward_ihsmm: several tracks", {
  ID <- rep(1:2, c(12, 18))
  first <- max(Ni_test)
  # time-varying dm
  for(omega in list(omega_test, omega_t_test)) {
    ref <- 0
    for(i in 1:2) {
      ind <- which(ID == i)
      rows <- ind[first:length(ind)]
      dm_i <- lapply(dm_inhom, function(d) d[ind, ])
      Gamma_t <- function(t) ref_Gamma_t(t - ind[1] + 1, if(is.matrix(omega)) omega else omega[, , t], dm_i)
      ref <- ref + ref_hsmm(ref_delta_stat(Gamma_t(rows[1])), Gamma_t, rows)
    }
    expect_hsmm(function(ap) forward_ihsmm(dm_inhom, omega, ap, trackID = ID, report = FALSE), ref)
  }
  # fixed dm
  Gamma_t <- function(t) ref_Gamma(omega_t_test[, , t], lapply(dm_hom, ref_hazard))
  ref <- ref_hsmm(ref_delta_stat(Gamma_t(1)), Gamma_t, 1:12) + ref_hsmm(ref_delta_stat(Gamma_t(13)), Gamma_t, 13:30)
  expect_hsmm(function(ap) forward_ihsmm(dm_hom, omega_t_test, ap, trackID = ID, report = FALSE), ref)
})

test_that("tpm_ihsmm and tpm_phsmm: aggregates of size 1 and 2", {
  for(Ni in list(c(5, 2, 4), c(5, 1, 4))) {
    dm <- dm_poisson(beta0_test, beta1_test, z_test, Ni)
    Gamma <- tpm_ihsmm(omega_test, dm)
    for(t in 5:30) expect_close(as.matrix(Gamma[[t - 4]]), ref_Gamma_t(t, omega_test, dm))
    dm <- lapply(dm, function(d) d[1:L_test, , drop = FALSE])
    Gamma <- tpm_phsmm(omega_test, dm)
    for(t in 1:L_test) expect_close(as.matrix(Gamma[[t]]), ref_Gamma_t(t, omega_test, dm, L_test))
  }
})

test_that("forward_phsmm: tod does not cover the full cycle", {
  Gamma_t <- function(t) ref_Gamma_t(t, omega_test, dm_per, L_test)
  # short series, times of day 5 and 6 are never observed
  tod <- 1:4
  delta <- ref_delta_stat(Reduce(`%*%`, lapply(1:L_test, Gamma_t)))
  expect_close(forward_phsmm(dm_per, omega_test, ap[1:4, ], tod, report = FALSE),
               ref_hsmm(delta, function(t) Gamma_t(tod[t]), 1:4))
})

test_that("forward_phsmm: one and several tracks", {
  Gamma_s <- function(s) ref_Gamma_t(s, omega_p_test[, , s], dm_per, L_test)
  delta <- ref_delta_stat(Reduce(`%*%`, lapply(1:L_test, Gamma_s)))
  expect_hsmm(function(ap) forward_phsmm(dm_per, omega_p_test, ap, tod_test, trackID = rep(1, 30), report = FALSE),
              ref_hsmm(delta, function(t) Gamma_s(tod_test[t])))
  # second track starts at time of day 3
  ID <- rep(1:2, c(12, 18))
  tod <- c(rep(1:6, 2), 3:6, rep(1:6, 2), 1:2)
  delta2 <- ref_delta_stat(Reduce(`%*%`, lapply(c(3:6, 1:2), Gamma_s)))
  expect_hsmm(function(ap) forward_phsmm(dm_per, omega_p_test, ap, tod, trackID = ID, report = FALSE),
              ref_hsmm(delta, function(t) Gamma_s(tod[t]), 1:12) + ref_hsmm(delta2, function(t) Gamma_s(tod[t]), 13:30))
})

test_that("forward_phsmm: reported quantities, omega convention and tod", {
  par <- list(beta0 = beta0_test, beta1 = beta1_test, eta = cbind(c(0.3, -0.2, 0.1), c(0.1, 0.2, 0.3)))
  zp <- sin(2 * pi * 1:L_test / L_test)
  omega <- tpm_emb_g(cbind(1, zp), par$eta, report = FALSE)
  nll <- function(par) {
    dm <- dm_poisson(par$beta0, par$beta1, zp, Ni_test)
    -forward_phsmm(dm, tpm_emb_g(cbind(1, zp), par$eta, report = FALSE), ap, tod_test)
  }
  rep <- RTMB::MakeADFun(nll, par, silent = TRUE)$report()
  expect_equal(rep$Gamma_sparse, tpm_phsmm(omega, dm_per))
  delta <- ref_delta_stat(Reduce(`%*%`, lapply(1:L_test, function(s) ref_Gamma_t(s, omega[, , s], dm_per, L_test))))
  expect_close(as.matrix(rep$delta_sparse), delta, tolerance = 1e-6)
  expect_close(rep$omega, omega) # as passed
  expect_equal(forward_phsmm(dm_per, omega, ap, tod_test, report = FALSE),
               forward_hsmm(dm_per, omega[, , c(6, 1:5)], ap, tod = tod_test, report = FALSE))
  # periodic dm needs regular tod
  expect_error(forward_phsmm(dm_per, omega, ap, rep(c(2, 5), 15), report = FALSE), "advance by one")
})

test_that("forward_hsmm: tracks of length one", {
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ID <- rep(1:2, c(1, 29))
  ref <- ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID)
  expect_close(forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE), ref)
  expect_close(forward_hsmm(dm_hom, omega_test, ap[1, , drop = FALSE], report = FALSE),
               ref_hsmm(ref_delta_stat(Gamma), function(t) Gamma, 1))
})
