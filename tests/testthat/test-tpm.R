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
