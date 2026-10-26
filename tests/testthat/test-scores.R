lp <- logP_test
w <- seq(0.5, 2, length.out = 20) # arbitrary weights

test_that("forward: weighted log-likelihood contributions are exact", {
  ID <- rep(1:3, c(5, 8, 7)); tod <- rep(1:4, 5)
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2))
  # f(allprobs) and Gam(t, i) as in ref_forward()
  expect_weighted <- function(f, delta, Gam, trackID = rep(1, 20), logspace = FALSE) {
    contrib <- ref_contrib(delta, Gam, lp, trackID)
    expect_equal(forward_weighted(f, lp, w, logspace), sum(w * contrib))
    # single contribution
    expect_equal(forward_weighted(f, lp, replace(numeric(20), 9, 1), logspace), contrib[9])
  }
  expect_weighted(function(ap) forward(d_test, G1, ap, report = FALSE), d_test, function(t, i) G1)
  expect_weighted(function(ap) forward(D_test, Gk_test, ap, trackID = ID, report = FALSE),
                  D_test, function(t, i) Gk_test[, , i], ID)
  expect_weighted(function(ap) forward(d_test, Gt_test, ap, trackID = ID, report = FALSE),
                  d_test, function(t, i) Gt_test[, , t], ID)
  expect_weighted(function(ap) forward(d_test, Gp, ap, report = FALSE, tod = tod),
                  d_test, function(t, i) Gp[, , tod[t]])
  expect_weighted(function(ap) forward(d_test, G1, ap, logspace = TRUE, report = FALSE),
                  d_test, function(t, i) G1, logspace = TRUE)
  # banded with bandwidth larger than the tracks: exact
  expect_weighted(function(ap) forward(d_test, G1, ap, trackID = ID, bw = 10, report = FALSE),
                  d_test, function(t, i) G1, ID)
})

test_that("forward: weights for one state and with weights 1", {
  lp1 <- lp[, 1, drop = FALSE]
  expect_equal(forward_weighted(function(ap) forward(1, matrix(1), ap, logspace = TRUE, report = FALSE),
                                lp1, w, logspace = TRUE), sum(w * lp1))
  f <- function(ap) forward(d_test, Gt_test, ap, report = FALSE)
  expect_equal(forward_weighted(f, lp, rep(1, 20)), ref_forward(d_test, function(t, i) Gt_test[, , t], lp))
})

test_that("forward: several forward() calls take consecutive blocks of weights", {
  ID <- rep(1:2, c(5, 15))
  contrib <- ref_contrib(D_test[1:2, ], function(t, i) G1, lp, ID)
  f <- function(ap) forward(D_test[1, ], G1, ap[1:5, ], report = FALSE) +
    forward(D_test[2, ], G1, ap[6:20, ], report = FALSE)
  expect_equal(forward_weighted(f, lp, w), sum(w * contrib))
  # more observations than weights
  expect_error(forward_weighted(f, lp, w[1:10]), "More observations")
})

test_that("count_obs: observations of all forward() calls", {
  nll <- function(par) {
    ap <- exp(lp + par$x)
    -forward(D_test[1, ], G1, ap[1:5, ], report = FALSE) - forward(D_test[2, ], G1, ap[6:20, ], report = FALSE)
  }
  expect_equal(count_obs(nll, list(x = 0)), 20)
  nll <- function(par) -forward(d_test, G1, exp(lp + par$x), trackID = rep(1:2, each = 10), report = FALSE)
  expect_equal(count_obs(nll, list(x = 0)), 20)
  expect_null(getFromNamespace(".lama_env", "LaMa")$count) # counting mode ends
})

test_that("forward: weights are only used with AD", {
  env <- getFromNamespace(".lama_env", "LaMa")
  env$weights <- w
  env$pos <- 0
  on.exit(env$weights <- NULL)
  expect_error(forward(d_test, G1, exp(lp), ad = FALSE), "only used with AD")
})

# small model with time-varying Gamma (matrix parameter beta) and two tracks
set.seed(4)
dat_sc <- list(x = rnorm(30, rep(c(0, 3, 0), each = 10)), z = runif(30), ID = rep(1:2, each = 15))
par_sc <- list(mu = c(0, 3), log_sigma = c(0, 0), beta = matrix(c(-2, -2, 0.5, -0.5), 2))
nll_sc <- function(par) {
  getAll(par, dat_sc)
  Gamma <- tpm(beta, z, report = FALSE)
  allprobs <- cbind(dnorm(x, mu[1], exp(log_sigma[1])), dnorm(x, mu[2], exp(log_sigma[2])))
  -forward(c(0.5, 0.5), Gamma, allprobs, trackID = ID)
}
est_sc <- c(-0.2, 2.8, 0.1, -0.1, -1.5, -2.5, 0.3, -0.2)

obj_sc <- RTMB::MakeADFun(nll_sc, par_sc, silent = TRUE)

test_that("scores: column sums equal the gradient", {
  U <- scores(obj_sc, est_sc)
  expect_equal(dim(U), c(30, 8))
  expect_equal(colSums(U), -obj_sc$gr(est_sc)[1, ], ignore_attr = TRUE)
})

test_that("scores: rows are derivatives of the log-likelihood contributions", {
  # exact contributions as a function of the parameters
  contrib <- function(theta) {
    par <- utils::relist(theta, par_sc)
    Gamma <- tpm(par$beta, dat_sc$z, report = FALSE)
    logP <- cbind(dnorm(dat_sc$x, par$mu[1], exp(par$log_sigma[1]), log = TRUE),
                  dnorm(dat_sc$x, par$mu[2], exp(par$log_sigma[2]), log = TRUE))
    ref_contrib(c(0.5, 0.5), function(t, i) Gamma[, , t], logP, dat_sc$ID)
  }
  h <- 1e-5
  num <- sapply(seq_along(est_sc), function(k) {
    e <- replace(0 * est_sc, k, h)
    (contrib(est_sc + e) - contrib(est_sc - e)) / (2 * h)
  })
  expect_equal(unname(scores(obj_sc, est_sc)), num, tolerance = 1e-6)
})

test_that("scores: respects map", {
  obj <- RTMB::MakeADFun(nll_sc, par_sc, map = list(log_sigma = factor(c(1, 1))), silent = TRUE)
  est <- est_sc[-4] # one common log_sigma
  U <- scores(obj, est)
  expect_equal(ncol(U), 7)
  expect_equal(colSums(U), -obj$gr(est)[1, ], ignore_attr = TRUE)
})

test_that("scores: model object from lama() and best parameters", {
  obj <- RTMB::MakeADFun(nll_sc, par_sc, silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  U <- scores(lama(obj))
  expect_equal(U, scores(obj, opt$par))
  expect_equal(unname(colSums(U)), rep(0, 8), tolerance = 1e-4) # MLE
})

test_that("scores: one forward() call per track gives the same", {
  nll2 <- function(par) {
    getAll(par, dat_sc)
    Gamma <- tpm(beta, z, report = FALSE)
    allprobs <- cbind(dnorm(x, mu[1], exp(log_sigma[1])), dnorm(x, mu[2], exp(log_sigma[2])))
    -forward(c(0.5, 0.5), Gamma[, , 1:15], allprobs[1:15, ], report = FALSE) -
      forward(c(0.5, 0.5), Gamma[, , 16:30], allprobs[16:30, ], report = FALSE)
  }
  obj2 <- RTMB::MakeADFun(nll2, par_sc, silent = TRUE)
  expect_equal(scores(obj2, est_sc), scores(obj_sc, est_sc))
})

test_that("scores: informative errors", {
  expect_error(scores(obj_sc, est_sc[-1]), "same length")
  obj0 <- RTMB::MakeADFun(function(par) sum(par$x^2), list(x = 1), silent = TRUE)
  expect_error(scores(obj0), "no call of forward")
  obj_re <- RTMB::MakeADFun(nll_sc, par_sc, random = "mu", silent = TRUE)
  expect_error(scores(obj_re), "random effects")
  expect_error(scores(list(a = 1)), "Expected an RTMB object")
})
