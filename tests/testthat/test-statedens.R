step <- trex$step[1:100]
angle <- trex$angle[1:100]
mu <- c(0.3, 2.5); sigma <- c(0.2, 1.5); kappa <- c(0.2, 1.5)

# allprobs of one stream by hand
by_hand <- function(x, f) {
  P <- matrix(1, length(x), 2)
  ind <- which(!is.na(x))
  for(j in 1:2) P[ind, j] <- f(x[ind], j, ind)
  P
}

test_that("statedens() equals the loop by hand", {
  expect_equal(statedens(step, dgamma2, mu, sigma),
               by_hand(step, function(x, j, ind) dgamma2(x, mu[j], sigma[j])))
  # single value for all states, named arguments
  expect_equal(statedens(angle, dvm, mu = 0, kappa = kappa),
               by_hand(angle, function(x, j, ind) dvm(x, 0, kappa[j])))
})

test_that("statedens(): parameters per observation and state", {
  Mu <- cbind(seq(0.2, 0.4, length = 100), seq(2, 3, length = 100))
  expect_equal(statedens(step, dgamma2, Mu, sigma),
               by_hand(step, function(x, j, ind) dgamma2(x, Mu[ind, j], sigma[j])))
})

test_that("statedens(): missing values per stream", {
  step[3] <- NA; angle[5] <- NA
  P <- statedens(step, dgamma2, mu, sigma) * statedens(angle, dvm, 0, kappa)
  expect_equal(P[3, ], dvm(angle[3], 0, kappa))
  expect_equal(P[5, ], dgamma2(step[5], mu, sigma))
})

test_that("statedens(): nStates for single values, and with AD", {
  expect_equal(dim(statedens(step, dgamma2, 1, 1, nStates = 3)), c(100, 3))
  F <- RTMB::MakeTape(function(p) sum(statedens(step, dgamma2, exp(p[1:2]), exp(p[3:4]))), log(c(mu, sigma)))
  expect_equal(F(log(c(mu, sigma))), sum(statedens(step, dgamma2, mu, sigma)))
})

test_that("statedens(log = TRUE): log-densities, also without a log argument", {
  step[3] <- NA
  expect_equal(statedens(step, dgamma2, mu, sigma, log = TRUE),
               log(statedens(step, dgamma2, mu, sigma)))
  # zero_inflate() has no log argument, so log() is taken afterwards
  zig <- zero_inflate(dgamma2)
  step[1] <- 0
  expect_equal(statedens(step, zig, 0.1, mu, sigma, log = TRUE),
               log(statedens(step, zig, 0.1, mu, sigma)))
  # same likelihood as on the natural scale
  Gamma <- tpm(c(-2, -2)); delta <- stationary(Gamma)
  expect_equal(forward(delta, Gamma, statedens(step, zig, 0.1, mu, sigma, log = TRUE), logspace = TRUE),
               forward(delta, Gamma, statedens(step, zig, 0.1, mu, sigma)))
})

# two-stream model built with statedens(), for the tests below
d_sd <- list(step = trex$step[1:300], angle = trex$angle[1:300])
nll_sd <- function(par) {
  getAll(par, d_sd)
  Gamma <- tpm(eta)
  allprobs <- statedens(step, dgamma2, exp(log_mu), exp(log_sigma)) * statedens(angle, dvm, 0, exp(log_kappa))
  -forward(stationary(Gamma), Gamma, allprobs)
}
par_sd <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), log_kappa = log(c(0.2, 1.5)))
obj_sd <- RTMB::MakeADFun(nll_sd, par_sd, silent = TRUE)
opt_sd <- nlminb(obj_sd$par, obj_sd$fn, obj_sd$gr)
mod_sd <- lama(obj_sd)

test_that("forward() reports the streams of statedens() and clears them", {
  expect_named(mod_sd$streams, c("step", "angle"))
  expect_equal(mod_sd$streams$step$name, "dgamma2")
  expect_equal(mod_sd$streams$angle$pars[[2]], exp(mod_sd$par$log_kappa))
  expect_null(getFromNamespace(".lama_env", "LaMa")$streams)
})

test_that("stream_curves(): weighted densities add up to the marginal", {
  st <- get_streams(mod_sd)
  for(w in list("stateprobs", "viterbi", "stationary", c(0.5, 0.5))) {
    W <- state_weights(mod_sd, w)
    cv <- stream_curves(st$step, W)
    expect_equal(sum(cv$weights), 1)
    expect_equal(sum(rowSums(cv$density)) * diff(cv$x[1:2]), 1, tolerance = 0.01) # integrates to ~1
  }
  # parameters per observation: same as constant ones repeated
  st2 <- st$step
  st2$pars[[1]] <- matrix(st2$pars[[1]], 300, 2, byrow = TRUE)
  expect_equal(stream_curves(st2, W)$density, stream_curves(st$step, W)$density)
})

test_that("get_streams(): only if allprobs is the product of the streams", {
  mod <- mod_sd
  mod$allprobs[1, ] <- 0.5 * mod$allprobs[1, ]
  expect_error(get_streams(mod), "not the product")
  mod$streams <- NULL
  expect_error(get_streams(mod), "statedens")
})

test_that("plot(mod, type = \"statedens\") returns the curves", {
  pdf(NULL)
  on.exit(dev.off())
  cv <- plot(mod_sd, type = "statedens")
  expect_named(cv, c("step", "angle"))
  expect_equal(par("mfrow"), c(1, 1))
  cv <- plot(mod_sd, type = "statedens", stream = "angle", weights = "viterbi", breaks = 10, legend = FALSE)
  expect_equal(range(cv$x), c(-pi, pi))
  expect_error(plot(mod_sd, type = "statedens", stream = "speed"), "one of")
  expect_s3_class(plot(mod_sd, stream = "step"), "LaMaResiduals") # residuals of one stream
})

test_that("forward() reports the streams of statedens() called in its argument", {
  x <- trex$step[1:50]; z <- seq(0, 1, length = 50)
  nll <- function(par) {
    getAll(par)
    Gamma <- tpm(eta)
    -forward(stationary(Gamma), Gamma, statedens(x, dnorm, cbind(a[1] + b * z, a[2]), 1))
  }
  obj <- RTMB::MakeADFun(nll, list(eta = c(-2, -2), a = c(0, 3), b = 0.5), silent = TRUE)
  mod <- lama(obj) # first report: used to hold the streams of the taping instead
  expect_equal(mod$streams$x$pars[[1]], cbind(0 + 0.5 * z, 3))
})
