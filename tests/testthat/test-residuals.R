test_that("forward: record mode stores allprobs on the natural scale", {
  env <- getFromNamespace(".lama_env", "LaMa")
  env$record <- list()
  on.exit(env$record <- NULL)
  ap <- exp(logP_test)
  expect_equal(forward(d_test, G1, ap), 0)
  forward(d_test, G1, logP_test, logspace = TRUE)
  expect_length(env$record, 2)
  expect_equal(env$record[[1]], ap)
  expect_equal(env$record[[2]], ap)
})

# gamma HMM: loop over states, covariate-dependent mean, missing values, prior on eta
dat_res <- list(step = trex$step[1:300], z = seq(-1, 1, length.out = 300))
dat_res$step[c(10, 50)] <- NA
par_res <- list(log_mu = log(c(0.3, 2.5)), b = c(0.2, -0.1), log_sigma = log(c(0.2, 1.5)), eta = rep(-2, 2))
nll_res <- function(par) {
  getAll(par, dat_res)
  Gamma <- tpm(eta)
  sigma <- exp(log_sigma)
  allprobs <- matrix(1, length(step), 2)
  ind <- which(!is.na(step))
  for(j in 1:2) {
    mu <- exp(log_mu[j] + b[j] * z[ind])
    allprobs[ind, j] <- dgamma2(step[ind], mu, sigma[j])
  }
  -forward(stationary(Gamma), Gamma, allprobs) - sum(dnorm(eta, -2, 2, log = TRUE)) # prior
}
obj_res <- RTMB::MakeADFun(nll_res, par_res, silent = TRUE)

test_that("stream_cdfs: CDFs in the cells of allprobs", {
  sc <- stream_cdfs(obj_res)
  expect_named(sc, "step") # the prior on eta does not enter allprobs
  mu <- sapply(1:2, function(j) exp(par_res$log_mu[j] + par_res$b[j] * dat_res$z))
  F <- sapply(1:2, function(j) pgamma2(dat_res$step, mu[, j], exp(par_res$log_sigma[j])))
  obs <- !is.na(dat_res$step)
  expect_equal(sc$step$upper[obs, ], F[obs, ])
  expect_equal(sc$step$observed, obs)
  expect_null(sc$step$lower)
})

test_that("stream_cdfs: F(x - 1) for discrete families", {
  cnt <- rpois(50, 3)
  obj <- RTMB::MakeADFun(function(par) {
    getAll(par); Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, cbind(dpois(cnt, exp(l[1])), dpois(cnt, exp(l[2]))))
  }, list(l = log(c(1, 5)), eta = rep(-2, 2)), silent = TRUE)
  sc <- stream_cdfs(obj)
  expect_equal(sc$cnt$lower, cbind(ppois(cnt - 1, 1), ppois(cnt - 1, 5)))
})

test_that("stream_cdfs: stops if densities do not enter allprobs as product or mixture", {
  pos <- dat_res$step[!is.na(dat_res$step)]
  zc <- c(0, 0, 3, 5, 0, 2, 7, 1)
  g <- seq(-3, 6, length.out = 20)
  obj_of <- function(nll, par) RTMB::MakeADFun(nll, c(par, list(eta = rep(-2, 2))), silent = TRUE)
  # truncation by hand
  expect_error(stream_cdfs(obj_of(function(par) {
    getAll(par); Gamma <- tpm(eta)
    allprobs <- sapply(1:2, function(j) dnorm(pos, mu[j], 1) / (1 - pnorm(0, mu[j], 1)))
    -forward(stationary(Gamma), Gamma, allprobs)
  }, list(mu = c(0.3, 2)))), "not 1")
  # Jacobian by hand
  expect_error(stream_cdfs(obj_of(function(par) {
    getAll(par); Gamma <- tpm(eta)
    allprobs <- sapply(1:2, function(j) dnorm(log(pos), mu[j], 1) / pos)
    -forward(stationary(Gamma), Gamma, allprobs)
  }, list(mu = c(-1, 0.5)))), "not 1")
  # zero inflation via an indicator
  expect_error(stream_cdfs(obj_of(function(par) {
    getAll(par); Gamma <- tpm(eta)
    allprobs <- sapply(1:2, function(j) 0.3 * (zc == 0) + 0.7 * dpois(zc, exp(l[j])))
    -forward(stationary(Gamma), Gamma, allprobs)
  }, list(l = log(c(1, 5))))), "not 1")
  # quadrature over a latent variable
  expect_error(stream_cdfs(obj_of(function(par) {
    getAll(par); Gamma <- tpm(eta)
    allprobs <- sapply(1:2, function(j)
      sapply(pos[1:30], function(y) sum(dnorm(y, g, 0.5) * dnorm(g, mu[j], 1) * 0.47)))
    -forward(stationary(Gamma), Gamma, allprobs)
  }, list(mu = c(0, 3)))), "not 1")
  # density entering twice
  expect_error(stream_cdfs(obj_of(function(par) {
    getAll(par); Gamma <- tpm(eta)
    allprobs <- sapply(1:2, function(j) dnorm(pos, mu[j], 1)^2)
    -forward(stationary(Gamma), Gamma, allprobs)
  }, list(mu = c(0.3, 2)))), "non-linearly")
})

test_that("stream_cdfs: densities in helper functions are found", {
  pos <- dat_res$step[!is.na(dat_res$step)]
  helper <- function(mu) sapply(1:2, function(j) dnorm(pos, mu[j], 1))
  obj <- RTMB::MakeADFun(function(par) {
    getAll(par); Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, helper(mu))
  }, list(mu = c(0.3, 2), eta = rep(-2, 2)), silent = TRUE)
  expect_equal(stream_cdfs(obj)$pos$upper, cbind(pnorm(pos, 0.3, 1), pnorm(pos, 2, 1)))
})

test_that("stream_cdfs: only one forward() call", {
  obj <- RTMB::MakeADFun(function(par) {
    ap <- cbind(dnorm(dat_res$z, par$mu[1]), dnorm(dat_res$z, par$mu[2]))
    -forward(c(0.5, 0.5), G1, ap[1:100, ]) - forward(c(0.5, 0.5), G1, ap[101:300, ])
  }, list(mu = c(0, 1)), silent = TRUE)
  expect_error(stream_cdfs(obj), "one call of forward")
})

test_that("pseudo_res: explicit distribution and parameters", {
  sp <- stateprobs(d_test, G1, exp(logP_test), forecast = TRUE)
  x <- c(rnorm(19, 1), NA)
  r <- pseudo_res(x, "norm", list(mean = c(0, 3), sd = c(1, 1)), sp)
  expect_s3_class(r, "LaMaResiduals")
  expect_equal(unclass(r), qnorm(rowSums(sp * cbind(pnorm(x, 0, 1), pnorm(x, 3, 1)))))
  expect_true(is.na(r[20]))
  # discrete, randomised
  cnt <- c(rpois(19, 2), NA)
  set.seed(5)
  r <- suppressMessages(pseudo_res(cnt, "pois", list(lambda = c(1, 4)), sp))
  up <- rowSums(sp * cbind(ppois(cnt, 1), ppois(cnt, 4)))
  lo <- rowSums(sp * cbind(ppois(cnt - 1, 1), ppois(cnt - 1, 4)))
  set.seed(5)
  expect_equal(unclass(r), c(qnorm(runif(19, lo[1:19], up[1:19])), NA))
  # discrete, not randomised: NA where not observed
  r <- suppressMessages(pseudo_res(cnt, "pois", list(lambda = c(1, 4)), sp, randomise = FALSE, normal = FALSE))
  expect_equal(r$upper, replace(up, 20, NA))
  expect_equal(r$lower, replace(lo, 20, NA))
})

test_that("pseudo_res: parameters per observation with missing observations", {
  sp <- stateprobs(d_test, G1, exp(logP_test), forecast = TRUE)
  x <- replace(rnorm(20, 1), c(3, 8), NA)
  m <- cbind(seq(0, 1, length.out = 20), 3) # mean per observation and state
  ref <- qnorm(rowSums(sp * cbind(pnorm(x, m[, 1]), pnorm(x, m[, 2]))))
  expect_equal(unclass(pseudo_res(x, "norm", list(mean = m, sd = 1), sp)), ref)
  # vector per observation, common to all states
  sd <- seq(1, 2, length.out = 20)
  ref <- qnorm(rowSums(sp * cbind(pnorm(x, 0, sd), pnorm(x, 3, sd))))
  expect_equal(unclass(pseudo_res(x, "norm", list(mean = c(0, 3), sd = sd), sp)), ref)
})

fit_mod <- function(nll, par) {
  obj <- RTMB::MakeADFun(nll, par, silent = TRUE)
  nlminb(obj$par, obj$fn, obj$gr)
  lama(obj)
}

test_that("pseudo_res(mod): same as with family and parameters given", {
  mod <- fit_mod(nll_res, par_res) # loop over states, covariate, missing values, prior
  p <- mod$par
  mu <- sapply(1:2, function(j) exp(p$log_mu[j] + p$b[j] * dat_res$z))
  ref <- pseudo_res(dat_res$step, "gamma2", list(mean = mu, sd = exp(p$log_sigma)), mod = mod)
  expect_equal(pseudo_res(mod), ref)
  expect_equal(pseudo_res(mod = mod), ref)
  expect_equal(pseudo_res(mod$obj), ref) # RTMB object
  expect_true(all(is.na(pseudo_res(mod)[c(10, 50)])))
})

test_that("pseudo_res(mod): several data streams, one residual vector each", {
  d2 <- list(step = trex$step[1:300], angle = trex$angle[1:300])
  mod <- fit_mod(function(par) {
    getAll(par, d2)
    Gamma <- tpm(eta)
    mu <- exp(log_mu); sigma <- exp(log_sigma); kappa <- exp(log_kappa)
    allprobs <- cbind(dgamma2(step, mu[1], sigma[1]), dgamma2(step, mu[2], sigma[2]))
    ind <- which(!is.na(angle))
    for(j in 1:2) allprobs[ind, j] <- allprobs[ind, j] * dvm(angle[ind], 0, kappa[j])
    -forward(stationary(Gamma), Gamma, allprobs)
  }, list(log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), log_kappa = log(c(0.5, 2)), eta = rep(-2, 2)))
  p <- mod$par
  r <- pseudo_res(mod)
  expect_named(r, c("step", "angle"))
  expect_s3_class(r, "LaMaResiduals")
  expect_equal(r$step, pseudo_res(d2$step, "gamma2", list(mean = exp(p$log_mu), sd = exp(p$log_sigma)), mod = mod))
  expect_equal(r$angle, pseudo_res(d2$angle, "vm", list(mu = 0, kappa = exp(p$log_kappa)), mod = mod))
  # plot: one row per stream, mfrow restored afterwards
  pdf(NULL)
  expect_null(plot(r))
  expect_equal(par("mfrow"), c(1, 1))
  expect_equal(plot(mod), r) # plot.LaMaModel returns the residuals
  dev.off()
})

test_that("pseudo_res(mod): two streams with statedens(), also on the log scale", {
  d2 <- list(step = trex$step[1:300], angle = trex$angle[1:300])
  par <- list(log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), log_kappa = log(c(0.5, 2)), eta = rep(-2, 2))
  for(lg in c(FALSE, TRUE)) {
    mod <- fit_mod(function(par) {
      getAll(par, d2)
      Gamma <- tpm(eta)
      P1 <- statedens(step, dgamma2, exp(log_mu), exp(log_sigma), log = lg)
      P2 <- statedens(angle, dvm, 0, exp(log_kappa), log = lg)
      -forward(stationary(Gamma), Gamma, if(lg) P1 + P2 else P1 * P2, logspace = lg)
    }, par)
    p <- mod$par
    r <- pseudo_res(mod)
    expect_equal(r$step, pseudo_res(d2$step, "gamma2", list(mean = exp(p$log_mu), sd = exp(p$log_sigma)), mod = mod))
    expect_equal(r$angle, pseudo_res(d2$angle, "vm", list(mu = 0, kappa = exp(p$log_kappa)), mod = mod))
  }
})

test_that("pseudo_res(mod): loop over time, two streams summed on the log scale", {
  set.seed(3)
  s <- rep(1:2, each = 25, length.out = 300)
  d5 <- list(x = rnorm(300, c(0, 3)[s], 1), y = rgamma(300, shape = c(2, 6)[s], rate = 1))
  mod <- fit_mod(function(par) {
    getAll(par, d5)
    Gamma <- tpm(eta)
    lallprobs <- matrix(0, length(x), 2)
    for(t in seq_along(x)) lallprobs[t, ] <- dnorm(x[t], mu, exp(log_sd), log = TRUE) +
      dgamma(y[t], shape = exp(log_shape), rate = 1, log = TRUE)
    -forward(stationary(Gamma), Gamma, lallprobs, logspace = TRUE)
  }, list(mu = c(0, 3), log_sd = c(0, 0), log_shape = log(c(2, 6)), eta = rep(-2, 2)))
  p <- mod$par
  r <- pseudo_res(mod)
  expect_equal(r$x, pseudo_res(d5$x, "norm", list(mean = p$mu, sd = exp(p$log_sd)), mod = mod))
  expect_equal(r$y, pseudo_res(d5$y, "gamma", list(shape = exp(p$log_shape), rate = 1), mod = mod))
})

test_that("pseudo_res(mod): discrete, randomised as with family and parameters given", {
  set.seed(1)
  cnt <- rpois(200, rep(c(2, 8), each = 20))
  mod <- fit_mod(function(par) {
    getAll(par)
    Gamma <- tpm(eta)
    lambda <- exp(log_lambda)
    -forward(stationary(Gamma), Gamma, cbind(dpois(cnt, lambda[1]), dpois(cnt, lambda[2])))
  }, list(log_lambda = log(c(1, 5)), eta = rep(-2, 2)))
  set.seed(2); r <- suppressMessages(pseudo_res(mod))
  set.seed(2); ref <- suppressMessages(pseudo_res(cnt, "pois", list(lambda = exp(mod$par$log_lambda)), mod = mod))
  expect_equal(r, ref)
  r <- suppressMessages(pseudo_res(mod, randomise = FALSE))
  expect_named(r, c("lower", "upper", "mean"))
})

test_that("pseudo_res(mod) and scores(mod) for qreml() models use the final fit", {
  modmat <- modmat(~ s(tod, bs = "cp"), data = data.frame(tod = 1:24), knots = list(tod = c(0, 24)))
  dat <- list(step = trex$step[1:500], angle = trex$angle[1:500], tod = trex$tod[1:500],
              Z = modmat$Z, S = modmat$S, lambda = rep(100, 2))
  pnll <- function(par) {
    getAll(par, dat)
    Gamma <- tpm(cbind(beta0, betaSpline), Z)
    Delta <- stationary_p(Gamma)
    mu <- exp(logmu); REPORT(mu)
    sigma <- exp(logsigma); REPORT(sigma)
    kappa <- exp(logkappa); REPORT(kappa)
    allprobs <- matrix(1, length(step), 2)
    ind <- which(!is.na(step) & !is.na(angle))
    for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], mu[j], sigma[j]) * dvm(angle[ind], 0, kappa[j])
    -forward(Delta[tod[1], ], Gamma[, , tod], allprobs) + penalty(betaSpline, S, lambda)
  }
  par <- list(logmu = log(c(0.3, 2.5)), logsigma = log(c(0.3, 2)), logkappa = log(c(0.1, 2)),
              beta0 = c(-2, -2), betaSpline = matrix(0, 2, ncol(modmat$Z) - 1))
  mod <- suppressMessages(qreml(pnll, par, dat, random = "betaSpline", silent = 2))
  expect_equal(est_par(mod), mod$par_vec)
  r <- pseudo_res(mod)
  ok <- !is.na(dat$step) & !is.na(dat$angle)
  ref <- pseudo_res(dat$step, "gamma2", list(mean = mod$mu, sd = mod$sigma), mod = mod)
  expect_equal(unclass(r$step)[ok], unclass(ref)[ok])
  ref <- pseudo_res(dat$angle, "vm", list(mu = 0, kappa = mod$kappa), mod = mod)
  expect_equal(unclass(r$angle)[ok], unclass(ref)[ok])
  expect_equal(scores(mod), scores(mod$obj, mod$par_vec))
})
