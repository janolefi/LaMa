# small 2-state model with several local optima
dat_ms <- list(step = trex$step[1:300])
par_ms <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)))
nll_ms <- function(par) {
  getAll(par, dat_ms)
  Gamma <- tpm(eta)
  allprobs <- cbind(dgamma2(step, exp(log_mu[1]), exp(log_sigma[1])), dgamma2(step, exp(log_mu[2]), exp(log_sigma[2])))
  -forward(stationary(Gamma), Gamma, allprobs)
}

test_that("multistart: best run, lama() builds the model from it", {
  obj <- RTMB::MakeADFun(nll_ms, par_ms, silent = TRUE)
  opt <- suppressWarnings(multistart(obj, nstarts = 6, stop = "none", seed = 1, silent = TRUE))
  expect_equal(nrow(opt$runs), 6)
  expect_equal(dim(opt$starts), c(6, length(obj$par)))
  expect_equal(opt$starts[1, ], unname(obj$par)) # first start is obj$par
  expect_equal(opt$objective, min(opt$runs$nll, na.rm = TRUE))
  expect_equal(lama(obj)$llk, -opt$objective)
})

test_that("multistart: seed reproduces starts and keeps the random number stream", {
  obj <- RTMB::MakeADFun(nll_ms, par_ms, silent = TRUE)
  set.seed(3); r1 <- runif(1)
  set.seed(3)
  o1 <- suppressWarnings(multistart(obj, nstarts = 3, stop = "none", seed = 1, silent = TRUE))
  expect_equal(runif(1), r1)
  o2 <- suppressWarnings(multistart(obj, nstarts = 3, stop = "none", seed = 1, silent = TRUE))
  expect_equal(o1$starts, o2$starts)
})

test_that("multistart: sd per parameter block and given starts", {
  obj <- RTMB::MakeADFun(nll_ms, par_ms, silent = TRUE)
  o <- suppressWarnings(multistart(obj, nstarts = 3, sd = list(eta = 0, log_mu = 0, log_sigma = 1),
                                   stop = "none", seed = 1, silent = TRUE))
  expect_equal(o$starts[, 1:4], matrix(obj$par[1:4], 3, 4, byrow = TRUE), ignore_attr = TRUE)
  expect_error(multistart(obj, sd = list(eta = 1)), "sd needs a value")
  S <- rbind(obj$par, obj$par + 0.1)
  o <- multistart(obj, starts = S, stop = "none", silent = TRUE)
  expect_equal(o$starts, unname(S))
  expect_error(multistart(obj, starts = S[, 1:2]), "one column per element")
})

test_that("multistart: stopping rules", {
  obj <- RTMB::MakeADFun(nll_ms, par_ms, silent = TRUE)
  # starts at the same point: one optimum, hits stops after 3
  S <- matrix(obj$par, 20, length(obj$par), byrow = TRUE)
  expect_equal(nrow(multistart(obj, starts = S, stop = "hits", silent = TRUE)$runs), 3)
  # Bayes: not for given starts, else after 8 random starts besides obj$par
  expect_equal(nrow(multistart(obj, starts = S, stop = "bayes", silent = TRUE)$runs), 20)
  expect_equal(nrow(multistart(obj, sd = 0, stop = "bayes", silent = TRUE)$runs), 9)
  expect_equal(n_optima(c(1, 1 + 1e-9, 2), 1e-6), 2)
  expect_output(multistart(obj, starts = S[1:2, ], stop = "none"), "2 starts, 2 converged, 1 distinct optima")
})

test_that("multistart: unit information scatter", {
  obj <- RTMB::MakeADFun(nll_ms, par_ms, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  Sigma <- unit_info(obj, opt$par, sd_max = 2)
  H <- obj$he(opt$par)
  s <- sqrt(diag(300 * solve(H)))
  expect_equal(sqrt(diag(Sigma)), pmin(s, 2))
  expect_equal(cov2cor(Sigma), cov2cor(solve(H)))
  o <- suppressWarnings(multistart(obj, nstarts = 4, scale = "unit_info", stop = "none", seed = 1, silent = TRUE))
  expect_equal(nrow(o$runs), 4)
})
