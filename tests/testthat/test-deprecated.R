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

test_that("stationary_sparse and stationary_p_sparse are deprecated and give the same results", {
  dm <- list(dpois(0:9, 3), dpois(0:7, 2))
  Gamma <- tpm_hsmm(tpm_switch(), dm)
  expect_message(delta <- stationary_sparse(Gamma), "deprecated")
  expect_equal(delta, stationary(as.matrix(Gamma)))
  expect_equal(as.vector(delta %*% as.matrix(Gamma)), as.vector(delta))
  dmp <- lapply(c(3, 2), function(l) sapply(0:9, dpois, lambda = l * (1 + 0.3 * sin(2 * pi * 1:4 / 4))))
  Gp <- tpm_phsmm(tpm_switch(), dmp)
  Gp_array <- simplify2array(lapply(Gp, as.matrix))
  expect_message(Delta <- stationary_p_sparse(Gp), "deprecated")
  expect_equal(Delta, stationary_p(Gp_array))
  expect_equal(stationary_p_sparse(Gp, 3), stationary_p(Gp_array, 3))
  # with AD
  f <- function(x) sum(stationary_p_sparse(lapply(Gp, function(G) as.matrix(G) * (1 + 0 * x))) * 1:40)
  expect_equal(RTMB::MakeADFun(function(par) f(par$x), list(x = 1), silent = TRUE)$fn(), f(1))
})
