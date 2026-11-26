test_that("old HSMM functions are deprecated and give the same results", {
  omega <- tpm_switch()
  dm <- list(dpois(0:9, 3), dpois(0:9, 5))
  G <- tpm_hsmm(omega, dm, sparse = FALSE)
  ap <- matrix(runif(40), 20, 2)
  d <- stationary(G)
  expect_message(G2 <- tpm_hsmm2(omega, dm), "deprecated")
  expect_message(l <- forward_s(d, G, ap, c(10, 10)), "deprecated")
  expect_equal(l, forward(d, G, ap[, rep(1:2, c(10, 10))]))

  dmp <- lapply(dm, function(p) matrix(p, 4, 10, byrow = TRUE)) # same in all 4 slices of the cycle
  expect_message(Gp <- tpm_phsmm2(omega, dmp), "deprecated")
  expect_message(l <- forward_sp(stationary_p(Gp)[1, ], Gp, ap, c(10, 10), rep(1:4, 5)), "deprecated")
  expect_equal(l, forward(d, G, ap[, rep(1:2, c(10, 10))]))
})

test_that("forward_ihsmm is deprecated and gives the same result as forward_hsmm", {
  omega <- tpm_switch()
  dm <- lapply(c(3, 5), function(lambda) sapply(0:9, dpois, lambda = rep(lambda, 20)))
  ap <- matrix(runif(40), 20, 2)
  expect_message(l <- forward_ihsmm(dm, omega, ap), "deprecated")
  expect_equal(l, forward_hsmm(dm, omega, ap))
})

test_that("forward_phsmm is deprecated and gives the same result as forward_hsmm", {
  omega <- tpm_switch()
  dm <- lapply(c(3, 5), function(lambda) sapply(0:9, dpois, lambda = lambda * 1:4))
  ap <- matrix(runif(40), 20, 2)
  tod <- rep(1:4, 5)
  expect_message(l <- forward_phsmm(dm, omega, ap, tod), "deprecated")
  expect_equal(l, forward_hsmm(dm, omega, ap, tod = tod))
})

test_that("other old functions are deprecated and give the same results", {
  Z <- cbind(1, seq(-1, 1, length = 10))
  beta <- matrix(c(-2, -2, 1, -1), 2)
  expect_message(G <- tpm_g2(Z, beta, report = FALSE), "deprecated")
  expect_equal(G, tpm(beta, Z, report = FALSE), ignore_attr = TRUE)

  Q <- generator(c(-1, -2))
  expect_message(P <- stationary_cont(Q), "deprecated")
  expect_equal(P, stationary_ct(Q))
  expect_message(G <- tpm_cont(Q, c(1, 2), report = FALSE), "deprecated")
  expect_equal(G, tpm_ct(Q, c(1, 2), report = FALSE))

  expect_message(expect_error(qreml_old("not a function", list(), list(), "b")), "deprecated")
})

test_that("tpm_emb and tpm_emb_g are deprecated and give the same results as tpm_switch", {
  expect_message(omega <- tpm_emb(c(0.3, -0.2, 0.1)), "deprecated")
  expect_equal(omega, tpm_switch(c(0.3, -0.2, 0.1)))
  beta <- cbind(c(0.3, -0.2, 0.1), c(0.4, -0.5, 0.2))
  expect_message(omega <- tpm_emb_g(cbind(1, 1:5), beta), "deprecated")
  expect_equal(omega, tpm_switch(beta, 1:5))
  expect_equal(tpm_emb_g(rep(1, 5), beta[, 1]), tpm_switch(beta[, 1, drop = FALSE], rep(1, 5))) # vector beta: no covariates
})
