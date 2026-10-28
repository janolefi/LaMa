test_that("sim_states(): initial distribution, transitions and tracks", {
  set.seed(1)
  Gamma <- tpm(c(-2, -1.5)); delta <- c(0.9, 0.1)
  inp <- prepare_decoding(delta, Gamma, matrix(1, 20000, 2), rep(1:200, each = 100))
  Z <- sim_states(inp)
  expect_equal(mean(Z[inp$start] == 1), 0.9, tolerance = 0.05)
  inner <- setdiff(2:20000, inp$start)
  tr <- table(Z[inner - 1], Z[inner]) / as.vector(table(Z[inner - 1]))
  expect_equal(unname(tr), unname(Gamma), tolerance = 0.02)
})

test_that("sim_states(): time-varying Gamma uses the slice into t", {
  set.seed(2)
  n <- 40000
  G <- array(0, c(2, 2, n))
  for(t in 1:n) G[, , t] <- if(t %% 2) tpm(c(-4, -4)) else tpm(c(4, 4)) # odd: stay, even: switch
  Z <- sim_states(prepare_decoding(c(0.5, 0.5), G, matrix(1, n, 2)))
  expect_gt(mean(Z[seq(3, n, 2)] == Z[seq(2, n - 1, 2)]), 0.95)
  expect_gt(mean(Z[seq(2, n, 2)] != Z[seq(1, n - 1, 2)]), 0.95)
})

test_that("sim_states(): track-specific, periodic, continuous-time Gamma and one state", {
  set.seed(4)
  n <- 20000
  # share of transitions 1 -> 2 among the transitions into the given rows
  p12 <- function(Z, rows) mean(Z[rows][Z[rows - 1] == 1] == 2)
  G2 <- array(c(tpm(c(-2, -2)), tpm(c(1, 1))), c(2, 2, 2)) # two different matrices

  # one Gamma per track
  Z <- sim_states(prepare_decoding(c(0.5, 0.5), G2, matrix(1, n, 2), rep(1:2, each = n / 2)))
  expect_equal(c(p12(Z, 2:(n / 2)), p12(Z, (n / 2 + 2):n)), c(G2[1, 2, 1], G2[1, 2, 2]), tolerance = 0.1)

  # periodic: tod indexes the slices
  tod <- rep(1:2, length.out = n)
  Z <- sim_states(prepare_decoding(c(0.5, 0.5), G2, matrix(1, n, 2), tod = tod))
  expect_equal(c(p12(Z, which(tod == 1)[-1]), p12(Z, which(tod == 2))), c(G2[1, 2, 1], G2[1, 2, 2]), tolerance = 0.1)

  # continuous time: nObs - 1 matrices, slice t - 1 is the transition into t
  timediff <- rexp(n - 1, 2)
  Gt <- tpm_ct(generator(c(-1, -1.5)), timediff)
  Z <- sim_states(prepare_decoding(c(0.5, 0.5), Gt, matrix(1, n, 2)))
  from1 <- which(Z[-n] == 1)
  expect_equal(mean(Z[from1 + 1] == 2), mean(Gt[1, 2, from1]), tolerance = 0.1)

  # one state
  expect_equal(sim_states(prepare_decoding(1, matrix(1), matrix(1, 5, 1))), rep(1L, 5))
})

test_that("sim_stream(): per state, per observation, missing values, r-function lookup", {
  set.seed(3)
  n <- 20000
  Z <- rep(1:2, length.out = n)
  x <- rep(1, n); x[7] <- NA
  st <- list(x = x, dens = dnorm, name = "dnorm", pars = list(c(0, 10), 1), nStates = 2)
  y <- sim_stream(st, Z)
  expect_true(is.na(y[7]) && sum(is.na(y)) == 1)
  expect_equal(c(mean(y[Z == 1], na.rm = TRUE), mean(y[Z == 2], na.rm = TRUE)), c(0, 10), tolerance = 0.05)
  # parameters per observation and state, and a different namespace
  st$pars[[1]] <- cbind(seq(0, 1, length = n), seq(10, 11, length = n))
  y <- sim_stream(st, Z)
  expect_equal(mean(y[Z == 2 & seq_len(n) > n / 2], na.rm = TRUE), 10.75, tolerance = 0.05)
  st$dens <- dgamma2; st$name <- "LaMa::dgamma2" # from another namespace, name with prefix
  st$pars <- list(c(1, 5), c(0.5, 0.5))
  expect_equal(mean(sim_stream(st, Z)[Z == 2], na.rm = TRUE), 5, tolerance = 0.05)
  # no counterpart: error, unless given
  zig <- zero_inflate(dnorm)
  st$dens <- zig; st$name <- "zig"; st$pars <- list(0, c(0, 10), 1)
  expect_error(sim_stream(st, Z), "rdens")
  rzig <- function(n, zeroprob, mean, sd) rnorm(n, mean, sd)
  expect_equal(mean(sim_stream(st, Z, rzig)[Z == 2], na.rm = TRUE), 10, tolerance = 0.05)
})

# two tracks with a missing value, built with statedens()
d_sim <- list(step = trex$step[1:400], angle = trex$angle[1:400], id = rep(c("b", "a"), each = 200))
d_sim$step[10] <- NA
nll_sim <- function(par) {
  getAll(par, d_sim)
  Gamma <- tpm(eta)
  allprobs <- statedens(step, dgamma2, exp(log_mu), exp(log_sigma)) * statedens(angle, dvm, 0, exp(log_kappa))
  -forward(stationary(Gamma), Gamma, allprobs, trackID = id)
}
par_sim <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), log_kappa = log(c(0.2, 1.5)))
obj_sim <- RTMB::MakeADFun(nll_sim, par_sim, silent = TRUE)
nlminb(obj_sim$par, obj_sim$fn, obj_sim$gr)
mod_sim <- lama(obj_sim)

test_that("simulate(mod): observed design, missing values, seed", {
  sim <- simulate(mod_sim, seed = 1)
  expect_named(sim, c("trackID", "state", "step", "angle"))
  expect_equal(sim$trackID, rep(1:2, each = 200))
  expect_equal(is.na(sim$step), is.na(d_sim$step))
  expect_equal(is.na(sim$angle), is.na(d_sim$angle))
  expect_true(all(sim$step[!is.na(sim$step)] > 0) && all(abs(sim$angle[!is.na(sim$angle)]) <= pi))
  # seed: reproducible, and the random number stream of the session is left alone
  set.seed(5); a <- runif(1); set.seed(5)
  expect_equal(simulate(mod_sim, seed = 1), sim)
  expect_equal(runif(1), a)
  # several simulations
  expect_length(simulate(mod_sim, nsim = 3), 3)
})

test_that("simulate(mod): data are simulated from the fitted model", {
  mod <- mod_sim
  mod$Gamma <- tpm(c(-4, -4)); mod$delta <- c(1, 0)
  mod$streams$step$pars <- list(c(1, 10), c(0.1, 0.1))
  mod$allprobs <- Reduce(`*`, lapply(mod$streams, function(s) dens_matrix(s$x, s$dens, s$pars, s$nStates)))
  sim <- simulate(mod, nsim = 1, seed = 2)
  expect_equal(as.vector(tapply(sim$step, sim$state, mean, na.rm = TRUE)), c(1, 10), tolerance = 0.1)
})

test_that("simulate(mod): errors without statedens() or for densities without generator", {
  mod <- mod_sim
  mod$streams <- NULL
  expect_error(simulate(mod), "statedens")
  mod <- mod_sim
  mod$streams$step$name <- "dmy"
  expect_error(simulate(mod), "rdens")
  expect_equal(is.na(simulate(mod, rdens = list(step = rgamma2))$step), is.na(d_sim$step))
})

test_that("simulate(mod): all-NA streams as placeholders, simulate then fit", {
  d <- d_sim; d$step[] <- NA; d$angle[] <- NA
  nll <- nll_sim; environment(nll) <- list2env(list(d_sim = d), parent = environment(nll_sim))
  par <- list(eta = c(-2, -2), log_mu = log(c(1, 10)), log_sigma = log(c(0.5, 2)), log_kappa = log(c(1, 5)))
  mod0 <- lama(RTMB::MakeADFun(nll, par, silent = TRUE))
  sim <- simulate(mod0, seed = 1)
  expect_false(anyNA(sim))
  expect_equal(nrow(sim), 400)
  # data from the chosen parameters
  expect_equal(as.vector(tapply(sim$step, sim$state, mean)), c(1, 10), tolerance = 0.3)
  # and the likelihood can be fitted to them
  d$step <- sim$step; d$angle <- sim$angle
  obj <- RTMB::MakeADFun(nll, par, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  expect_equal(opt$convergence, 0)
  expect_lt(max(abs(opt$par[3:4] - log(c(1, 10)))), 0.5)
})

test_that("simulate(mod): no trackID column for a single track, parameters per observation", {
  x <- rep(NA, 100); z <- seq(-1, 1, length = 100)
  nll <- function(par) {
    getAll(par)
    Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, statedens(x, dnorm, cbind(a[1] + b * z, a[2]), 1))
  }
  par <- list(eta = c(-3, -3), a = c(0, 10), b = 5)
  sim <- simulate(lama(RTMB::MakeADFun(nll, par, silent = TRUE)), seed = 1)
  expect_named(sim, c("state", "x"))
  expect_equal(nrow(sim), 100)
  s1 <- sim$state == 1 # mean 5 * z in state 1
  expect_gt(cor(sim$x[s1], z[s1]), 0.5)
})

# random effect on the mean of state 1 per track, for obj$simulate()
n_re <- 2000; id_re <- rep(1:20, each = 100)
step_re <- rep(1, n_re); step_re[5] <- NA
nll_re <- function(par) {
  getAll(par)
  u %~% dnorm(0, exp(log_su))
  mu <- cbind(exp(log_mu[1] + u[id_re]), exp(log_mu[2]))
  Gamma <- tpm(eta)
  .nll <- .nll - forward(stationary(Gamma), Gamma, statedens(step_re, dgamma2, mu, exp(log_sigma)), id_re)
  .nll
}
par_re <- list(eta = c(-3, -3), log_mu = log(c(1, 10)), log_sigma = log(c(0.2, 1)), log_su = log(0.5), u = rep(0, 20))
obj_re <- RTMB::MakeADFun(nll_re, par_re, random = "u", silent = TRUE)

test_that("obj$simulate(): random effects first, then states and data given them", {
  set.seed(1)
  sim <- obj_re$simulate(obj_re$env$last.par)
  expect_named(sim, c("u", "step_re", "states"))
  expect_equal(is.na(sim$step_re), is.na(step_re))
  expect_equal(as.vector(tapply(sim$step_re, sim$states, mean, na.rm = TRUE))[2], 10, tolerance = 0.1)
  # data in state 1 follow the simulated random effect of their track
  m1 <- tapply(log(sim$step_re[sim$states == 1]), id_re[sim$states == 1], mean, na.rm = TRUE)
  expect_gt(cor(m1, sim$u[as.integer(names(m1))]), 0.9)
  # new random effects in every simulation
  expect_false(isTRUE(all.equal(sim$u, obj_re$simulate(obj_re$env$last.par)$u)))
  # the likelihood itself is unchanged
  expect_true(is.finite(obj_re$fn(obj_re$par)))
})

test_that("obj$simulate(): data that are NA throughout are simulated completely", {
  nll <- nll_re; environment(nll) <- list2env(list(step_re = rep(NA_real_, n_re)), parent = environment(nll_re))
  obj <- RTMB::MakeADFun(nll, par_re, random = "u", silent = TRUE)
  sim <- obj$simulate(obj$env$last.par)
  expect_false(anyNA(sim$step_re))
})

test_that("simulate(mod, random = TRUE) uses obj$simulate()", {
  mod <- lama(obj_re)
  sim <- simulate(mod, random = TRUE, seed = 1)
  expect_named(sim, c("trackID", "state", "step_re"))
  expect_equal(is.na(sim$step_re), is.na(step_re))
  expect_equal(simulate(mod, random = TRUE, seed = 1), sim)
  expect_length(simulate(mod, random = TRUE, nsim = 2), 2)
  # fixed random effects (default) differ from freshly drawn ones: spread of the track means in state 1
  spread <- function(s) sd(tapply(log(s$step_re[s$state == 1]), s$trackID[s$state == 1], mean, na.rm = TRUE))
  set.seed(2)
  expect_gt(mean(sapply(simulate(mod, random = TRUE, nsim = 5), spread)), 0.3) # sd of u is 0.5
  mod <- mod_sim
  mod$obj <- RTMB::MakeADFun(function(par) -forward(c(0.5, 0.5), tpm(par), matrix(1, 5, 2)), c(-2, -2), silent = TRUE)
  expect_error(simulate(mod, random = TRUE), "statedens")
})

test_that("random = TRUE: density of the random effect as a sum, before its use", {
  nll_sum <- function(par) {
    getAll(par)
      nll <- -sum(dnorm(u, 0, exp(log_su), log = TRUE))
    mu <- cbind(exp(log_mu[1] + u[id_re]), exp(log_mu[2]))
    Gamma <- tpm(eta)
    nll - forward(stationary(Gamma), Gamma, statedens(step_re, dgamma2, mu, exp(log_sigma)), id_re)
  }
  obj <- RTMB::MakeADFun(nll_sum, par_re, random = "u", silent = TRUE)
  sim <- simulate(lama(obj), random = TRUE, seed = 1)
  expect_named(sim, c("trackID", "state", "step_re"))
  expect_equal(sim, simulate(lama(obj_re), random = TRUE, seed = 1)) # same as with %~%

  # density only after the use: random effect not simulated yet
  nll_late <- function(par) {
    getAll(par)
      mu <- cbind(exp(log_mu[1] + u[id_re]), exp(log_mu[2]))
    Gamma <- tpm(eta)
    nll <- -forward(stationary(Gamma), Gamma, statedens(step_re, dgamma2, mu, exp(log_sigma)), id_re)
    nll - sum(dnorm(u, 0, exp(log_su), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll_late, par_re, random = "u", silent = TRUE)
  expect_error(simulate(lama(obj), random = TRUE), "model hierarchy")
})

test_that("statedens(): informative error for random effects that are not simulated yet", {
  # random effect passed directly, and in a vector of parameters per state
  for(mean_par in list(quote(u), quote(exp(1 + u)))) {
    nll <- eval(bquote(function(par) {
      getAll(par)
          Gamma <- tpm(eta)
      nll <- -forward(stationary(Gamma), Gamma, statedens(step_re, dnorm, .(mean_par), 1))
      nll - sum(dnorm(u, 0, 1, log = TRUE))
    }))
    obj <- RTMB::MakeADFun(nll, list(eta = c(-3, -3), u = c(0, 5)), random = "u", silent = TRUE)
    expect_error(obj$simulate(), "model hierarchy")
  }
})

test_that("obj$simulate(): data in a list, explicit OBS(), no stale data after MakeTape", {
  dat <- list(step = replace(rep(1, 100), 7, NA))
  nll <- function(par) {
    getAll(par, dat)
    Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, statedens(step, dgamma2, exp(log_mu), 1))
  }
  par <- list(eta = c(-3, -3), log_mu = log(c(1, 10)))
  sim <- RTMB::MakeADFun(nll, par, silent = TRUE)$simulate()
  expect_named(sim, c("step", "states"))
  expect_equal(is.na(sim$step), is.na(dat$step))
  # same with the data as dat$step, and with OBS() called by the user
  nll2 <- function(par) {
    getAll(par)
    step <- OBS(dat$step)
    Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, statedens(dat$step, dgamma2, exp(log_mu), 1))
  }
  expect_named(RTMB::MakeADFun(nll2, par, silent = TRUE)$simulate(), c("step", "states"))
  # the registration of the data in a tape does not change later calls of statedens()
  step <- dat$step
  RTMB::MakeTape(function(p) sum(statedens(step, dgamma2, exp(p), 1)), log(c(1, 10)))
  step[3] <- NA
  expect_equal(statedens(step, dgamma2, c(1, 10), 1)[3, ], c(1, 1))
})

test_that("simulate(fill_na = TRUE) fills in missing observations", {
  expect_true(anyNA(simulate(mod_sim, seed = 1)$step))
  sim <- simulate(mod_sim, fill_na = TRUE, seed = 1)
  expect_false(anyNA(sim))
  # the other draws are the same: only the missing rows are added
  expect_equal(sim$step[!is.na(d_sim$step)], simulate(mod_sim, seed = 1)$step[!is.na(d_sim$step)])
  # with random = TRUE
  mod <- lama(obj_re)
  expect_true(anyNA(simulate(mod, random = TRUE, seed = 1)$step_re))
  expect_false(anyNA(simulate(mod, random = TRUE, fill_na = TRUE, seed = 1)))
})
