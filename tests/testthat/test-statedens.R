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
