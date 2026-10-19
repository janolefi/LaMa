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
  # density in a helper function outside the likelihood function
  helper <- function(mu) sapply(1:2, function(j) dnorm(pos, mu[j], 1))
  expect_error(stream_cdfs(obj_of(function(par) {
    getAll(par); Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, helper(mu))
  }, list(mu = c(0.3, 2)))), "No density")
})

test_that("stream_cdfs: only one forward() call", {
  obj <- RTMB::MakeADFun(function(par) {
    ap <- cbind(dnorm(dat_res$z, par$mu[1]), dnorm(dat_res$z, par$mu[2]))
    -forward(c(0.5, 0.5), G1, ap[1:100, ]) - forward(c(0.5, 0.5), G1, ap[101:300, ])
  }, list(mu = c(0, 1)), silent = TRUE)
  expect_error(stream_cdfs(obj), "one call of forward")
})
