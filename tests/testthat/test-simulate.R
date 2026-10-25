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
mod_sim <- report(obj_sim)

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
  mod0 <- report(RTMB::MakeADFun(nll, par, silent = TRUE))
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
  sim <- simulate(report(RTMB::MakeADFun(nll, par, silent = TRUE)), seed = 1)
  expect_named(sim, c("state", "x"))
  expect_equal(nrow(sim), 100)
  s1 <- sim$state == 1 # mean 5 * z in state 1
  expect_gt(cor(sim$x[s1], z[s1]), 0.5)
})
