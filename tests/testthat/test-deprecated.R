test_that("old HSMM functions are deprecated and give the same results", {
  omega <- tpm_emb()
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
