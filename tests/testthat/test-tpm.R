test_that("row sums are equal to one (for two states)", {
  param = rep(-1, 2)
  expect_equal(rowSums(tpm(param)), c(S1 = 1, S2 = 1))
})

test_that("row sums are equal to one (for three states)", {
  param = rep(-2, 6)
  expect_equal(rowSums(tpm(param)), c(S1 = 1, S2 = 1, S3 = 1))
})

test_that("tpm_g() and generator_g() are deprecated and give the same results", {
  z <- seq(-1, 1, length.out = 10); Z <- cbind(z, z^2)
  beta <- matrix(c(-2, -2, 1, -1, 0.5, 0.5), 2)
  expect_message(G <- tpm_g(Z, beta, report = FALSE), "deprecated")
  expect_equal(G, tpm(beta, Z, report = FALSE))
  expect_equal(dim(G), c(2, 2, 10))
  expect_message(Q <- generator_g(Z, beta, report = FALSE), "deprecated")
  expect_equal(Q, generator(beta, Z, report = FALSE))
  # a single covariate as a vector
  expect_equal(tpm(beta[, 1:2], z, report = FALSE), tpm(beta[, 1:2], cbind(1, z), report = FALSE))
})

test_that("deprecated(): message once per session", {
  expect_message(deprecated("old() is deprecated, test."), "old")
  expect_silent(deprecated("old() is deprecated, test."))
})

test_that("tpm_ct(): 2-state closed form (AD) agrees with Matrix::expm", {
  set.seed(1)
  n <- 30; dt <- rexp(n, 0.5); z <- rnorm(n)
  beta <- c(-1, -2, 0.5, -0.3)
  # old AD path: Matrix::expm per interval
  expm_loop <- function(Q, rates) {
    "[<-" <- RTMB::ADoverload("[<-")
    G <- array(0, c(2, 2, n))
    for(t in 1:n) G[, , t] <- as.matrix(Matrix::expm((Q[, , t] - diag(rates)) * dt[t]))
    G
  }
  for(rates in list(NULL, c(0.5, 3))) {
    r <- if(is.null(rates)) c(0, 0) else rates
    # constant Q (matrix) and time-varying Q (array)
    fQ <- function(b) generator(b[1:2], report = FALSE)
    gQ <- function(b) generator(matrix(b, 2), Z = z, report = FALSE)
    for(Qfun in list(fQ, gQ)) {
      f_new <- function(b) c(tpm_ct(Qfun(b), dt, rates, report = FALSE))
      f_old <- function(b) {
        Q <- Qfun(b)
        if(is.matrix(Q)) Q <- array(Q, c(2, 2, n))
        c(expm_loop(Q, r))
      }
      F_new <- RTMB::MakeTape(f_new, beta)
      F_old <- RTMB::MakeTape(f_old, beta)
      expect_equal(F_new(beta), c(expm_loop(array(Qfun(beta), c(2, 2, n)), r)), tolerance = 1e-12)
      expect_equal(F_new(beta), f_new(beta), tolerance = 1e-12) # AD vs C++
      expect_equal(F_new$jacobian(beta), F_old$jacobian(beta), tolerance = 1e-10)
    }
  }
})

test_that("tpm_ct(): 2-state closed form handles extreme rates and zero Q", {
  dt <- c(1e-8, 1, 1e3)
  F <- RTMB::MakeTape(function(b) c(tpm_ct(generator(b, report = FALSE), dt, ad = TRUE, report = FALSE)), c(0, 0))
  # exact: Gamma_12 = q12 (1 - exp(-s dt)) / s (expm loses accuracy for huge s dt)
  for(b in list(c(-15, -10), c(8, 9), c(-10, 8))) {
    q21 <- exp(b[1]); q12 <- exp(b[2]); s <- q12 + q21
    p12 <- -q12 * expm1(-s * dt) / s; p21 <- -q21 * expm1(-s * dt) / s
    expect_equal(F(b), c(rbind(1 - p12, p21, p12, 1 - p21)), tolerance = 1e-14)
  }
  # no switching at all: identity
  G <- RTMB::MakeTape(function(x) c(tpm_ct(x * matrix(0, 2, 2), dt, ad = TRUE, report = FALSE)), 1)
  expect_equal(G(1), rep(c(1, 0, 0, 1), 3))
})

test_that("tpm_ct(): constant Q with repeated time differences", {
  set.seed(2)
  n <- 40; dt <- sample(c(1, 2, 7, 0.5), n, replace = TRUE)
  for(N in 2:3) {
    b <- rnorm(N * (N - 1), -1)
    # reference: same Q as array (one slice per interval), no deduplication there
    f_ref <- function(b) c(tpm_ct(generator(matrix(b, ncol = 1), Z = matrix(1, n, 1), report = FALSE), dt, report = FALSE))
    f_new <- function(b) c(tpm_ct(generator(b, report = FALSE), dt, report = FALSE))
    G <- tpm_ct(generator(b, report = FALSE), dt, report = FALSE)
    expect_equal(dim(G), c(N, N, n))
    expect_equal(attr(G, "time"), "continuous")
    expect_equal(c(G), f_ref(b), tolerance = 1e-12) # C++
    F_new <- RTMB::MakeTape(f_new, b); F_ref <- RTMB::MakeTape(f_ref, b)
    expect_equal(F_new(b), f_ref(b), tolerance = 1e-12) # AD
    expect_equal(F_new$jacobian(b), F_ref$jacobian(b), tolerance = 1e-10)
  }
  # with rates
  Q <- generator(c(-1, -2), report = FALSE)
  expect_equal(tpm_ct(Q, dt, rates = c(0.5, 2), report = FALSE)[, , 5], 
               as.matrix(Matrix::expm((Q - diag(c(0.5, 2))) * dt[5])), tolerance = 1e-12, ignore_attr = TRUE)
})
