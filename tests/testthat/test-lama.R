# small fitted model
dat_rep <- list(step = replace(trex$step[1:200], c(5, 50), NA)) # two missing
par_rep <- list(log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), eta = rep(-2, 2))
nll_rep <- function(par) {
  getAll(par, dat_rep)
  Gamma <- tpm(eta)
  mu <- exp(log_mu); REPORT(mu); sigma <- exp(log_sigma); REPORT(sigma)
  allprobs <- matrix(1, length(step), 2)
  ind <- which(!is.na(step))
  for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], mu[j], sigma[j])
  -forward(stationary(Gamma), Gamma, allprobs)
}
obj_rep <- RTMB::MakeADFun(nll_rep, par_rep, silent = TRUE)
opt_rep <- nlminb(obj_rep$par, obj_rep$fn, obj_rep$gr)

test_that("lama: stores obj itself, get_obj() finds it", {
  mod <- lama(obj_rep)
  expect_s3_class(mod, "LaMaModel")
  expect_identical(mod$obj$env, obj_rep$env) # same environment, not a copy
  expect_identical(get_obj(mod), get_obj(obj_rep))
  expect_identical(get_obj(structure(list(obj = obj_rep), class = "qremlModel"))$env, obj_rep$env)
  expect_error(get_obj(list(a = 1)), "Expected an RTMB object")
})

test_that("lama: log-likelihood, missing observations not counted, gradient", {
  mod <- lama(obj_rep)
  expect_equal(mod$llk, -opt_rep$objective)
  expect_equal(mod$nobs, sum(!is.na(dat_rep$step)))
  expect_equal(attr(logLik(mod), "nobs"), mod$nobs)
  expect_lt(mod$mgc, 1e-3)
})

test_that("print: state process, missing observations and reported quantities", {
  mod <- lama(obj_rep)
  expect_output(print(mod), "2 states \\(discrete time, homogeneous\\)")
  expect_output(print(mod), "198 observations \\(2 missing\\), 6 parameters")
  expect_output(print(mod), "Reported: mu, sigma, delta, Gamma, allprobs")
})

test_that("summary: state process and reported quantities by state", {
  mod <- lama(obj_rep)
  s <- summary(mod)
  expect_s3_class(s, "summary.LaMaModel")
  expect_equal(unname(s$Gamma), unname(mod$Gamma))
  expect_equal(unname(s$stationary), unname(stationary(mod$Gamma)))
  expect_equal(unname(s$dwell), unname(1 / (1 - diag(mod$Gamma))))
  expect_equal(s$bystate, rbind(mu = mod$mu, sigma = mod$sigma), ignore_attr = TRUE)
  expect_output(print(s), "Mean dwell times")
})

test_that("summary: continuous time, averages over tracks, periodic", {
  Q <- matrix(c(-1, 1, 2, -2), 2, byrow = TRUE)
  mod <- structure(list(allprobs = matrix(0.5, 10, 2), type = "inhomogeneous", Q = Q,
                        llk = -10, df = 4, nobs = 10), class = "LaMaModel")
  s <- summary(mod)
  expect_equal(unname(s$Q), Q)
  expect_equal(unname(s$stationary), c(2, 1) / 3)
  expect_equal(unname(s$dwell), c(1, 0.5))

  G <- array(c(0.9, 0.2, 0.1, 0.8, 0.5, 0.5, 0.5, 0.5), c(2, 2, 2))
  mod <- structure(list(allprobs = matrix(0.5, 10, 2), type = "track-specific", Gamma = G,
                        trackID = rep(1:2, each = 5), llk = -10, df = 4, nobs = 10), class = "LaMaModel")
  s <- summary(mod)
  expect_equal(s$nTracks, 2)
  expect_equal(unname(s$Gamma), (G[,,1] + G[,,2]) / 2)
  expect_equal(unname(s$stationary), (c(2, 1) / 3 + c(0.5, 0.5)) / 2)
  expect_equal(unname(s$dwell), (c(10, 5) + c(2, 2)) / 2)

  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -2, 0.5, 0.5, 0, 0), 2))
  mod <- structure(list(allprobs = matrix(0.5, 8, 2), type = "periodic", Gamma = Gp, tod = rep(1:4, 2),
                        llk = -10, df = 4, nobs = 8), class = "LaMaModel")
  expect_equal(unname(summary(mod)$occupancy), unname(colMeans(stationary_p(Gp))))
})

test_that("summary: qreml models add smooth terms", {
  mod <- structure(list(allprobs = matrix(0.5, 10, 2), type = "inhomogeneous", llk = -10, df = 8.5, nobs = 10,
                        n_fixpar = 4, edf = list("s(x)" = c(2, 2.5)), lambda = c("s(x).1" = 1, "s(x).2" = 3),
                        spname = "lambda", parname = "par", converged = TRUE, iter = 3),
                   class = c("qremlModel", "LaMaModel"))
  s <- summary(mod)
  expect_s3_class(s, "summary.qremlModel")
  expect_equal(s$edf, c("s(x).1" = 2, "s(x).2" = 2.5))
  expect_null(s$other$lambda) # set by qreml(), not reported
  expect_output(print(s), "4 fixed parameters \\+ smooths with 4.5 edf")
  expect_output(print(s), "Converged: TRUE \\(3 outer iterations\\)")
})

test_that("MCreport: model object gives the same as obj", {
  mod <- lama(obj_rep)
  set.seed(1); a <- suppressMessages(MCreport(obj_rep, nSamples = 20))
  set.seed(1); b <- suppressMessages(MCreport(mod, nSamples = 20))
  expect_equal(a, b)
})

test_that("report: deprecated, same as lama()", {
  expect_warning(mod <- report(obj_rep), "deprecated")
  expect_equal(mod, lama(obj_rep))
})
