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

test_that("forward: weights are only used with AD", {
  env <- getFromNamespace(".lama_env", "LaMa")
  env$weights <- w
  on.exit(env$weights <- NULL)
  expect_error(forward(d_test, G1, exp(lp), ad = FALSE), "only used with AD")
})
