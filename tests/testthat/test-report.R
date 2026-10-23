# small fitted model
dat_rep <- list(step = replace(trex$step[1:200], c(5, 50), NA)) # two missing
par_rep <- list(log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), eta = rep(-2, 2))
nll_rep <- function(par) {
  getAll(par, dat_rep)
  Gamma <- tpm(eta)
  mu <- exp(log_mu); sigma <- exp(log_sigma)
  allprobs <- matrix(1, length(step), 2)
  ind <- which(!is.na(step))
  for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], mu[j], sigma[j])
  -forward(stationary(Gamma), Gamma, allprobs)
}
obj_rep <- RTMB::MakeADFun(nll_rep, par_rep, silent = TRUE)
opt_rep <- nlminb(obj_rep$par, obj_rep$fn, obj_rep$gr)

test_that("report: stores obj itself, get_obj() finds it", {
  mod <- report(obj_rep)
  expect_s3_class(mod, "LaMaModel")
  expect_identical(mod$obj$env, obj_rep$env) # same environment, not a copy
  expect_identical(get_obj(mod), get_obj(obj_rep))
  expect_identical(get_obj(structure(list(obj = obj_rep), class = "qremlModel"))$env, obj_rep$env)
  expect_error(get_obj(list(a = 1)), "Expected an RTMB object")
})

test_that("report: log-likelihood, missing observations not counted, gradient", {
  mod <- report(obj_rep)
  expect_equal(mod$llk, -opt_rep$objective)
  expect_equal(mod$nobs, sum(!is.na(dat_rep$step)))
  expect_equal(attr(logLik(mod), "nobs"), mod$nobs)
  expect_lt(mod$mgc, 1e-3)
})

test_that("MCreport: model object gives the same as obj", {
  mod <- report(obj_rep)
  set.seed(1); a <- suppressMessages(MCreport(obj_rep, nSamples = 20))
  set.seed(1); b <- suppressMessages(MCreport(mod, nSamples = 20))
  expect_equal(a, b)
})
