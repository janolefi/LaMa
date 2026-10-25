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
    -forward(delta, Gamma, allprobs, tod = tod_dec, trackID = trackID)
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
    # reporting of LaMa <= 2.1.1: expanded Gamma, with tod, type inhomogeneous
    G <- mod$Gamma[, , mod$tod]
    old <- mod; old$Gamma <- G; old$type <- "inhomogeneous"
    expect_equal(viterbi(mod = mod), viterbi(mod$delta, G, mod$allprobs, mod$trackID))
    expect_equal(viterbi(mod = mod), viterbi(mod = old))
    expect_equal(stateprobs(mod = mod), stateprobs(mod$delta, G, mod$allprobs, mod$trackID))
    expect_equal(stateprobs(mod = mod, forecast = TRUE),
                 stateprobs(mod$delta, G, mod$allprobs, mod$trackID, forecast = TRUE))
    expect_equal(pseudo_res(x_dec, "norm", list(mean = mod$mu, sd = 1), mod = mod),
                 pseudo_res(x_dec, "norm", list(mean = mod$mu, sd = 1), mod = old))
    # mod without type
    mod$type <- NULL
    expect_equal(viterbi(mod = mod), viterbi(mod = old))
    expect_equal(stateprobs(mod = mod), stateprobs(mod = old))
  }
})

test_that("decoding: stationary_p(Gamma, tod[1] %% L + 1) is periodically stationary under forward()'s tod", {
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -2, 1, -1, 0, 1), 2)); tod <- rep(c(3, 4, 1, 2), 5)
  # no information in allprobs: stateprobs are the marginal state distributions
  P <- stateprobs(stationary_p(Gp, tod[1] %% 4 + 1), Gp, matrix(1, 20, 2), tod = tod)
  expect_equal(unname(P), unname(stationary_p(Gp)[tod %% 4 + 1, ]))
})

test_that("decoding: forecast probabilities reproduce the log-likelihood", {
  ap <- exp(logP_test)
  llk_fc <- function(fc) sum(log(rowSums(fc * ap))) # sum of log f(x_t | x_1, ..., x_t-1)
  ID <- rep(1:3, c(2, 8, 10)) # with a track of length 2
  expect_equal(llk_fc(stateprobs(d_test, G1, ap, forecast = TRUE)),
               forward(d_test, G1, ap, ad = FALSE))
  expect_equal(llk_fc(stateprobs(d_test, Gt_test, ap, forecast = TRUE)),
               forward(d_test, Gt_test, ap, ad = FALSE))
  expect_equal(llk_fc(stateprobs(D_test, Gt_test, ap, trackID = ID, forecast = TRUE)),
               forward(D_test, Gt_test, ap, trackID = ID, ad = FALSE))
})

test_that("viterbi: exact for all kinds of Gamma, also for tracks of length 1 and 2", {
  ap <- exp(logP_test)
  ID <- rep(c("d", "a", "c", "b"), c(1, 2, 5, 12))
  D <- rbind(D_test, c(0.4, 0.6))
  Gk <- array(c(G1, G2, G3, G2), dim = c(2, 2, 4)) # track-specific
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2)); tod <- rep(1:4, 5)
  lp <- logP_test[1:12, ] # single track
  expect_equal(viterbi(d_test, G1, ap[1:12, ]), ref_viterbi(d_test, function(t, i) G1, lp))
  expect_equal(viterbi(d_test, Gt_test[, , 2:12], ap[1:12, ]),
               ref_viterbi(d_test, function(t, i) Gt_test[, , t], lp))
  expect_equal(viterbi(D, G1, ap, ID), ref_viterbi(D, function(t, i) G1, logP_test, ID))
  expect_equal(viterbi(D, Gk, ap, ID), ref_viterbi(D, function(t, i) Gk[, , i], logP_test, ID))
  expect_equal(viterbi(D, Gt_test, ap, ID), ref_viterbi(D, function(t, i) Gt_test[, , t], logP_test, ID))
  expect_equal(viterbi(D, Gp, ap, ID, tod = tod), ref_viterbi(D, function(t, i) Gp[, , tod[t]], logP_test, ID))
})

test_that("stateprobs: exact for all kinds of Gamma, also for tracks of length 1 and 2", {
  ap <- exp(logP_test)
  ID <- rep(c("d", "a", "c", "b"), c(1, 2, 5, 12))
  D <- rbind(D_test, c(0.4, 0.6))
  Gk <- array(c(G1, G2, G3, G2), dim = c(2, 2, 4)) # track-specific
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2)); tod <- rep(1:4, 5)
  # sp(forecast) calls a stateprobs function, Gam as in ref_forward()
  expect_stateprobs <- function(sp, delta, Gam, trackID = rep(1, 20)) {
    expect_equal(unname(sp(FALSE)), ref_stateprobs(delta, Gam, logP_test, trackID))
    expect_equal(unname(sp(TRUE)), ref_forecast(delta, Gam, logP_test, trackID))
  }
  expect_stateprobs(function(f) stateprobs(d_test, G1, ap, forecast = f), d_test, function(t, i) G1)
  expect_stateprobs(function(f) stateprobs(d_test, Gt_test[, , -1], ap, forecast = f),
                    d_test, function(t, i) Gt_test[, , t])
  expect_stateprobs(function(f) stateprobs(D, Gk, ap, ID, forecast = f), D, function(t, i) Gk[, , i], ID)
  expect_stateprobs(function(f) stateprobs(D, Gt_test, ap, ID, forecast = f),
                    D, function(t, i) Gt_test[, , t], ID)
  expect_stateprobs(function(f) stateprobs(D, Gp, ap, ID, forecast = f, tod = tod),
                    D, function(t, i) Gp[, , tod[t]], ID)
})

test_that("deprecated _g and _p functions give a message and the same results", {
  ap <- exp(logP_test); ID <- rep(1:3, c(5, 8, 7)); tod <- rep(1:4, 5)
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2))
  expect_message(l <- forward_g(d_test, Gt_test, ap, ad = FALSE), "deprecated")
  expect_equal(l, forward(d_test, Gt_test, ap, ad = FALSE))
  expect_message(l <- forward_p(d_test, Gp, ap, tod, ad = FALSE), "deprecated")
  expect_equal(l, forward(d_test, Gp, ap, ad = FALSE, tod = tod))
  expect_message(s <- viterbi_g(d_test, Gt_test, ap, ID), "deprecated")
  expect_equal(s, viterbi(d_test, Gt_test, ap, ID))
  expect_message(s <- viterbi_p(d_test, Gp, ap, tod, ID), "deprecated")
  expect_equal(s, viterbi(d_test, Gp, ap, ID, tod = tod))
  expect_message(p <- stateprobs_g(d_test, Gt_test, ap, ID, forecast = TRUE), "deprecated")
  expect_equal(p, stateprobs(d_test, Gt_test, ap, ID, forecast = TRUE))
  expect_message(p <- stateprobs_p(d_test, Gp, ap, tod, ID), "deprecated")
  expect_equal(p, stateprobs(d_test, Gp, ap, ID, tod = tod))
})

test_that("decoding: model object as the first argument", {
  mod <- fit_p(rep(1:3, each = 20))
  expect_equal(viterbi(mod), viterbi(mod = mod))
  expect_equal(stateprobs(mod, forecast = TRUE), stateprobs(mod = mod, forecast = TRUE))
  qmod <- structure(mod, class = "qremlModel") # as returned by qreml()
  expect_equal(viterbi(qmod), viterbi(mod = mod))
  expect_equal(viterbi(mod$delta, mod$Gamma, mod$allprobs, mod$trackID, tod = mod$tod), viterbi(mod))
})
