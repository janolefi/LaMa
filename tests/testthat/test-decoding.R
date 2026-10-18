# small periodic model with cycle length 6, fitted with RTMB
set.seed(3)
x_dec <- rnorm(60, rep(c(0, 3), each = 6), 1)
tod_dec <- rep(c(3:6, 1:2), 10)
fit_p <- function(trackID = NULL) {
  nll <- function(par) {
    getAll(par)
    Gamma <- tpm_p(1:6, L = 6, beta = beta)
    delta <- stationary_p(Gamma, tod_dec[1])
    allprobs <- cbind(dnorm(x_dec, mu[1], 1), dnorm(x_dec, mu[2], 1))
    REPORT(mu)
    -forward_p(delta, Gamma, allprobs, tod_dec, trackID = trackID)
  }
  obj <- RTMB::MakeADFun(nll, list(beta = matrix(c(-2, -2, 0, 0, 0, 0), 2), mu = c(0, 3)), silent = TRUE)
  nlminb(obj$par, obj$fn, obj$gr)
  obj$report()
}

test_that("decoding: periodic model reports L slices and decodes as with one slice per observation", {
  for(ID in list(NULL, rep(c("b", "a", "c"), each = 20))) {
    mod <- fit_p(ID)
    expect_equal(mod$type, "periodic")
    expect_equal(dim(mod$Gamma), c(2, 2, 6))
    # old reporting: expanded Gamma, type inhomogeneous
    G <- mod$Gamma[, , mod$tod]
    old <- mod; old$Gamma <- G; old$type <- "inhomogeneous"
    expect_equal(viterbi(mod = mod), viterbi_g(mod$delta, G, mod$allprobs, mod$trackID))
    expect_equal(viterbi(mod = mod), viterbi(mod = old))
    expect_equal(stateprobs(mod = mod), stateprobs_g(mod$delta, G, mod$allprobs, mod$trackID))
    expect_equal(stateprobs(mod = mod, forecast = TRUE),
                 stateprobs_g(mod$delta, G, mod$allprobs, mod$trackID, forecast = TRUE))
    expect_equal(pseudo_res(x_dec, "norm", list(mean = mod$mu, sd = 1), mod = mod),
                 pseudo_res(x_dec, "norm", list(mean = mod$mu, sd = 1), mod = old))
    # unknown type: dispatch on dim(Gamma)[3]
    old$type <- NULL
    expect_equal(viterbi(mod = old), viterbi(mod = mod))
    expect_equal(stateprobs(mod = old), stateprobs(mod = mod))
  }
})

test_that("decoding: forecast probabilities with trackID use the transition into t", {
  ap <- exp(logP_test)
  ID <- rep(1:2, each = 10)
  ref <- rbind(stateprobs_g(d_test, Gt_test[, , 2:10], ap[1:10, ], forecast = TRUE),
               stateprobs_g(d_test, Gt_test[, , 12:20], ap[11:20, ], forecast = TRUE))
  expect_equal(stateprobs_g(d_test, Gt_test, ap, trackID = ID, forecast = TRUE), ref)
})
