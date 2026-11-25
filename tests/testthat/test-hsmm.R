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

test_that("forward algorithm with hsmm_step: value and gradient with AD", {
  delta <- ref_delta_fresh(c(0.5, 0.3, 0.2), Ni_test)
  par <- c(beta0_test, beta1_test, 0.3, -0.2, 0.1)

  # with the direct step, all observations, a fresh stay starts at t = 1
  f_step <- function(par) {
    stay <- lapply(dm_poisson(par[1:3], par[4:6], z_test, Ni_test), hsmm_stay)
    omega <- tpm_emb(par[7:9])
    # a stay entered before t = 1 is impossible here, row 1 is never used for it
    stay_into <- function(t) {
      padded <- lapply(stay, function(s) s[c(rep(1, max(Ni_test)), 1:30), ])
      hsmm_lag(padded, t + max(Ni_test))[1, ]
    }
    hsmm_forward(delta, stay_into, function(t) omega, ap, Ni_test)
  }
  # exact reference, entries before t = 1 get probability 0 of ending
  f_ref <- function(par) {
    dm <- dm_poisson(par[1:3], par[4:6], z_test, Ni_test)
    omega <- tpm_emb(par[7:9])
    pad <- max(Ni_test)
    dm_pad <- lapply(dm, function(d) rbind(matrix(0, pad, ncol(d)), d))
    ref_forward(delta, function(t, i) ref_Gamma_t(t - 1 + pad, omega, dm_pad), lp_M)
  }
  F <- RTMB::MakeTape(f_step, par)
  expect_close(F(par), f_ref(par))
  h <- 1e-5
  num <- sapply(seq_along(par), function(k) {
    e <- replace(0 * par, k, h)
    (f_ref(par + e) - f_ref(par - e)) / (2 * h)
  })
  expect_close(F$jacobian(par), num, tolerance = 1e-5)
})

test_that("hsmm_step: hazards of stays entered before t = 1 do not matter for a given delta", {
  delta <- ref_delta_fresh(c(0.5, 0.3, 0.2), Ni_test)
  stay <- lapply(dm_inhom, hsmm_stay)
  pad <- max(Ni_test)
  loglik <- function(before) {
    padded <- lapply(stay, function(s) rbind(matrix(before, pad, ncol(s)), s))
    hsmm_forward(delta, function(t) hsmm_lag(padded, t + pad)[1, ],
                 function(t) omega_test, ap, Ni_test)
  }
  expect_close(loglik(0.1), loglik(0.9), tolerance = 1e-12)
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
    inp <- prepare_hsmm(dm_inhom, omega, ap, startInd = 10)
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

test_that("prepare_hsmm: use_all_obs is exact, stays before the first observation have no effect", {
  # same model with 10 rows of nonsense dm before the first observation
  before <- function(d, value) rbind(matrix(value, 10, ncol(d)), d)
  Gamma_into <- function(t, value) ref_Gamma_t(t - 1 + 10, omega_test, lapply(dm_inhom, before, value))
  inp <- prepare_hsmm(dm_inhom, omega_test, ap, delta = delta_N, use_all_obs = TRUE)
  l <- hsmm_forward(inp$delta[1, ], function(t) inp$stay[t, ], function(t) omega_test, ap, Ni_test)
  for(value in c(0.01, 0.1)) {
    expect_close(l, ref_forward(delta_M, function(t, i) Gamma_into(t, value), lp_M))
  }
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
  expect_error(prep(dm = dm_inhom, startInd = 6), "startInd")
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

test_that("gradients of the forward functions match the direct step", {
  par <- c(beta0_test, beta1_test, 0.3, -0.2, 0.1)
  first <- max(Ni_test)
  f_ihsmm <- function(par) {
    dm <- dm_poisson(par[1:3], par[4:6], z_test, Ni_test)
    forward_ihsmm(dm, tpm_emb(par[7:9]), ap, report = FALSE)
  }
  f_step <- function(par) {
    dm <- dm_poisson(par[1:3], par[4:6], z_test, Ni_test)
    omega <- tpm_emb(par[7:9])
    delta <- stationary(as.matrix(tpm_ihsmm(omega, dm)[[1]]))
    stay <- lapply(dm, hsmm_stay)
    hsmm_forward(as.vector(delta), function(t) hsmm_lag(stay, t)[1, ],
                 function(t) omega, ap, Ni_test, start = first)
  }
  F_ihsmm <- RTMB::MakeTape(f_ihsmm, par)
  F_step <- RTMB::MakeTape(f_step, par)
  expect_close(F_step(par), F_ihsmm(par), tolerance = 1e-10)
  expect_close(F_step$jacobian(par), F_ihsmm$jacobian(par), tolerance = 1e-8)
})

# ---- known bugs of the current functions, run with LAMA_KNOWN_BUGS=true ----

test_that("forward_ihsmm and forward_phsmm: fixed dm and given delta", {
  skip_known_bug("fixed dm with delta fails, 'k' not found")
  delta <- c(0.5, 0.3, 0.2)
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ref <- ref_hsmm(ref_delta_fresh(delta, Ni_test), function(t) Gamma)
  expect_close(forward_ihsmm(dm_hom, omega_test, ap, delta = delta, report = FALSE), ref)
  expect_close(forward_phsmm(dm_hom, omega_test, ap, tod_test, delta = delta, report = FALSE), ref)
})

test_that("forward_hsmm: delta can be estimated", {
  skip_known_bug("delta as an AD parameter fails")
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  delta <- c(0.5, 0.3, 0.2)
  F <- RTMB::MakeTape(function(delta) forward_hsmm(dm_hom, omega_test, ap, delta = delta, report = FALSE), delta)
  expect_close(F(delta), ref_hsmm(ref_delta_fresh(delta, Ni_test), function(t) Gamma))
})

test_that("forward_hsmm: several tracks with stationary delta under AD", {
  skip_known_bug("Matrix() does not accept the AD stationary distribution")
  ID <- rep(1:2, c(12, 18))
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ref <- ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID)
  expect_hsmm(function(ap) forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE), ref)
})

test_that("forward_hsmm: as many tracks as states", {
  skip_known_bug("dm is read as one list per track if nTracks = nStates")
  ID <- rep(1:3, c(8, 10, 12))
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  expect_close(forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE),
               ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID))
})

test_that("forward_ihsmm: several tracks", {
  skip_known_bug("AD() does not know lists")
  ID <- rep(1:2, c(12, 18))
  first <- max(Ni_test)
  ref <- 0
  for(i in 1:2) {
    ind <- which(ID == i)
    rows <- ind[first:length(ind)]
    dm_i <- lapply(dm_inhom, function(d) d[ind, ])
    Gamma_t <- function(t) ref_Gamma_t(t - ind[1] + 1, omega_test, dm_i)
    ref <- ref + ref_hsmm(ref_delta_stat(Gamma_t(rows[1])), Gamma_t, rows)
  }
  expect_hsmm(function(ap) forward_ihsmm(dm_inhom, omega_test, ap, trackID = ID, report = FALSE), ref)
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
  skip_known_bug("cycle length is taken as the number of unique tod values")
  Gamma_t <- function(t) ref_Gamma_t(t, omega_test, dm_per, L_test)
  # short series, times of day 5 and 6 are never observed
  tod <- 1:4
  delta <- ref_delta_stat(Reduce(`%*%`, lapply(1:L_test, Gamma_t)))
  expect_close(forward_phsmm(dm_per, omega_test, ap[1:4, ], tod, report = FALSE),
               ref_hsmm(delta, function(t) Gamma_t(tod[t]), 1:4))
})

test_that("forward_phsmm: a single track", {
  skip_known_bug("delta matrix drops to a vector")
  Gamma_t <- function(t) ref_Gamma_t(t, omega_test, dm_per, L_test)
  delta <- ref_delta_stat(Reduce(`%*%`, lapply(1:L_test, Gamma_t)))
  expect_close(forward_phsmm(dm_per, omega_test, ap, tod_test, trackID = rep(1, 30), report = FALSE),
               ref_hsmm(delta, function(t) Gamma_t(tod_test[t])))
})

test_that("forward_hsmm: tracks of length one", {
  skip_known_bug("loops 2:n run backwards for nObs = 1")
  Gamma <- ref_Gamma(omega_test, lapply(dm_hom, ref_hazard))
  ID <- rep(1:2, c(1, 29))
  ref <- ref_forward(ref_delta_stat(Gamma), function(t, i) Gamma, lp_M, ID)
  expect_close(forward_hsmm(dm_hom, omega_test, ap, trackID = ID, report = FALSE), ref)
  expect_close(forward_hsmm(dm_hom, omega_test, ap[1, , drop = FALSE], report = FALSE),
               ref_hsmm(ref_delta_stat(Gamma), function(t) Gamma, 1))
})
