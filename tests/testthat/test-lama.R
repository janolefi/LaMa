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

test_that("coef, vcov: estimates and inverse Hessian as from sdreport()", {
  mod <- lama(obj_rep)
  expect_equal(coef(mod), opt_rep$par, ignore_attr = TRUE)
  expect_named(coef(mod), c("log_mu.1", "log_mu.2", "log_sigma.1", "log_sigma.2", "eta.1", "eta.2"))
  V <- vcov(mod)
  expect_equal(V, RTMB::sdreport(obj_rep)$cov.fixed, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(rownames(V), names(coef(mod)))
  ci <- confint(mod)
  expect_equal(unname(rowMeans(ci)), unname(coef(mod))) # each row its own parameter
})

test_that("coef, vcov: random effects excluded, marginal Hessian", {
  nll <- function(par) {
    getAll(par, dat_rep)
    Gamma <- tpm(eta)
    mu <- exp(log_mu + u); sigma <- exp(log_sigma)
    allprobs <- matrix(1, length(step), 2)
    ind <- which(!is.na(step))
    for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], mu[j], sigma[j])
    -forward(stationary(Gamma), Gamma, allprobs) - sum(dnorm(u, 0, exp(log_tau), log = TRUE))
  }
  obj <- RTMB::MakeADFun(nll, c(par_rep, list(u = c(0, 0), log_tau = log(0.1))), random = "u", silent = TRUE)
  opt <- nlminb(obj$par, obj$fn, obj$gr)
  mod <- lama(obj)
  expect_equal(coef(mod), opt$par, tolerance = 1e-6, ignore_attr = TRUE)
  expect_true("log_tau" %in% names(coef(mod))) # single entries keep their name
  expect_equal(vcov(mod), RTMB::sdreport(obj)$cov.fixed, ignore_attr = TRUE, tolerance = 1e-4)
})

test_that("vcov: qreml models use the conditional Hessian, warning if not positive definite", {
  H <- obj_rep$he(opt_rep$par)
  mod <- structure(list(obj = obj_rep, parname = "par", par_vec = opt_rep$par, Hessian_conditional = H, mgc = 0),
                   class = c("qremlModel", "LaMaModel"))
  expect_equal(vcov(mod), solve(H), ignore_attr = TRUE)
  mod$Hessian_conditional <- -H
  expect_warning(vcov(mod), "not positive definite")
})

test_that("vcov HAC: single track as sandwich::meatHAC, both kernels, all lags", {
  skip_if_not_installed("sandwich")
  mod <- lama(obj_rep)
  U <- scores(mod)
  registerS3method("estfun", "lama_scores", function(x, ...) unclass(x), envir = asNamespace("sandwich"))
  x <- structure(U, class = "lama_scores")
  Hi <- solve(obj_rep$he(opt_rep$par))
  for(kernel in c("QS", "Bartlett")) {
    sw_kernel <- if(kernel == "QS") "Quadratic Spectral" else "Bartlett"
    for(bw in c(0.5, 1, 3.7, 25)) {
      # tol = 0: sandwich keeps all lags (by default it drops those with weight below 1e-7)
      w <- sandwich::weightsAndrews(x, bw = bw, kernel = sw_kernel, prewhite = 0, tol = 0)
      M <- nrow(U) * sandwich::meatHAC(x, weights = w, prewhite = 0, adjust = FALSE)
      expect_equal(hac_meat(U, kernel = kernel, bw = bw), M, ignore_attr = TRUE, tolerance = 1e-10)
      expect_equal(vcov(mod, type = "HAC", kernel = kernel, bw = bw), Hi %*% M %*% Hi, ignore_attr = TRUE, tolerance = 1e-6)
    }
    # default: Andrews bandwidth, as sandwich::bwAndrews without prewhitening and with equal weights for all
    # score columns (by default, sandwich guesses which columns belong to an intercept and gives them weight 0)
    b <- sandwich::bwAndrews(x, kernel = sw_kernel, prewhite = 0, weights = rep(1, ncol(U)))
    expect_equal(bw_andrews(U, list(seq_len(nrow(U))), kernel), b, tolerance = 1e-10)
    expect_equal(vcov(mod, type = "HAC", kernel = kernel), vcov(mod, type = "HAC", kernel = kernel, bw = b))
  }
  expect_equal(vcov(mod, type = "HAC"), vcov(mod, type = "HAC", kernel = "QS"))
  expect_equal(rownames(vcov(mod, type = "HAC")), names(coef(mod)))
  # Bartlett with bw <= 1, or any kernel with bw = 0: lag 0 only, the sandwich
  expect_equal(vcov(mod, type = "HAC", kernel = "Bartlett", bw = 1), vcov(mod, type = "sandwich"))
  expect_equal(vcov(mod, type = "HAC", bw = 0), vcov(mod, type = "sandwich"))
  expect_equal(vcov(mod, type = "sandwich"), Hi %*% crossprod(U) %*% Hi, ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("vcov HAC: FFT equals the direct sum over lags", {
  set.seed(1)
  for(n in c(1, 2, 5, 40)) {
    U <- matrix(rnorm(3 * n), n, 3)
    for(kernel in c("QS", "Bartlett")) {
      w <- hac_kernel((0:(n - 1)) / 2.3, kernel)
      M <- w[1] * crossprod(U)
      for(l in seq_len(n - 1)) {
        G <- crossprod(U[(l + 1):n, , drop = FALSE], U[1:(n - l), , drop = FALSE])
        M <- M + w[l + 1] * (G + t(G))
      }
      expect_equal(fft_meat(U, w), M, tolerance = 1e-10)
    }
  }
})

test_that("vcov HAC: autocovariances only within tracks", {
  dat2 <- list(step = trex$step[1:300], ID = rep(c("a", "b"), c(120, 180)))
  nll2 <- function(par) {
    getAll(par, dat2)
    Gamma <- tpm(eta)
    allprobs <- matrix(1, length(step), 2)
    ind <- which(!is.na(step))
    for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], exp(log_mu[j]), exp(log_sigma[j]))
    -forward(stationary(Gamma), Gamma, allprobs, trackID = ID)
  }
  obj2 <- RTMB::MakeADFun(nll2, par_rep, silent = TRUE)
  opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
  mod <- lama(obj2)
  U <- scores(mod)
  Hi <- solve(obj2$he(opt2$par))
  ind <- split(seq_len(nrow(U)), dat2$ID)
  for(kernel in c("QS", "Bartlett")) {
    M <- hac_meat(U[ind$a, ], kernel = kernel, bw = 5) + hac_meat(U[ind$b, ], kernel = kernel, bw = 5)
    expect_equal(hac_meat(U, dat2$ID, kernel, bw = 5), M)
    expect_equal(vcov(mod, type = "HAC", kernel = kernel, bw = 5), Hi %*% M %*% Hi, ignore_attr = TRUE, tolerance = 1e-6)
    # a single series would include pairs across the track boundary
    expect_false(isTRUE(all.equal(hac_meat(U, kernel = kernel, bw = 5), M)))
    # default bandwidth: AR(1) fits on pairs within tracks, mean track length 150
    expect_equal(vcov(mod, type = "HAC", kernel = kernel),
                 vcov(mod, type = "HAC", kernel = kernel, bw = bw_andrews(U, ind, kernel)))
  }
  # the Andrews bandwidth uses the mean track length
  U2 <- rbind(U, U)
  expect_equal(bw_andrews(U2, list(1:300, 301:600), "QS"), bw_andrews(U, list(1:300), "QS"), tolerance = 1e-10)
})

test_that("vcov HAC: errors for qreml models, random effects, bad bw, kernel and trackID", {
  mod <- lama(obj_rep)
  expect_error(vcov(mod, type = "HAC", bw = -1), "non-negative number")
  expect_error(vcov(mod, type = "HAC", bw = c(1, 2)), "non-negative number")
  expect_error(vcov(mod, type = "HAC", bw = Inf), "non-negative number")
  expect_error(vcov(mod, type = "HAC", bw = "a"), "non-negative number")
  expect_error(vcov(mod, type = "HAC", kernel = "Parzen"), "should be one of")
  expect_error(hac_meat(matrix(0, 10, 2), trackID = 1:5), "does not match")
  qmod <- structure(list(obj = obj_rep), class = c("qremlModel", "LaMaModel"))
  expect_error(vcov(qmod, type = "HAC"), "qreml")
  obj_re <- RTMB::MakeADFun(function(par) nll_rep(par[names(par_rep)]) - sum(dnorm(par$u, log = TRUE)),
                            c(par_rep, list(u = 0)), random = "u", silent = TRUE)
  expect_error(vcov(structure(list(obj = obj_re), class = "LaMaModel"), type = "HAC"), "random effects")
  expect_error(vcov(qmod, type = "cluster"), "qreml")
  expect_error(vcov(qmod, type = "sandwich"), "qreml")
  expect_error(vcov(structure(list(obj = obj_re), class = "LaMaModel"), type = "cluster"), "random effects")
})

# 8 tracks of 50 steps, character IDs
dat_cl <- list(step = trex$step[1:400], ID = rep(paste0("track", 1:8), each = 50))
nll_cl <- function(par) {
  getAll(par, dat_cl)
  Gamma <- tpm(eta)
  allprobs <- matrix(1, length(step), 2)
  ind <- which(!is.na(step))
  for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], exp(log_mu[j]), exp(log_sigma[j]))
  -forward(stationary(Gamma), Gamma, allprobs, trackID = ID)
}
obj_cl <- RTMB::MakeADFun(nll_cl, par_rep, silent = TRUE)
opt_cl <- nlminb(obj_cl$par, obj_cl$fn, obj_cl$gr)

test_that("vcov cluster: CR0 and CR1 from the score sums of the tracks", {
  mod <- lama(obj_cl)
  U <- scores(mod)
  Hi <- solve(obj_cl$he(opt_cl$par))
  S <- rowsum(U, dat_cl$ID)
  V0 <- Hi %*% crossprod(S) %*% Hi
  expect_equal(vcov(mod, type = "cluster", cr = "CR0"), V0, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(vcov(mod, type = "cluster", cr = "CR1"), 8 / 7 * V0, ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(vcov(mod, type = "cluster"), vcov(mod, type = "cluster", cr = "CR1")) # CR1 is the default
  expect_equal(rownames(vcov(mod, type = "cluster")), names(coef(mod)))
  # explicit cluster = tracks gives the same, whatever the type of the IDs
  expect_equal(vcov(mod, type = "cluster", cluster = factor(dat_cl$ID)), vcov(mod, type = "cluster"))
})

test_that("vcov cluster: clusters combining tracks, but not splitting them", {
  mod <- lama(obj_cl)
  U <- scores(mod)
  Hi <- solve(obj_cl$he(opt_cl$par))
  # 4 "individuals" with 2 tracks each (G = 4 < 7: singular, hence the warning)
  indiv <- rep(c("A", "B", "C", "D"), each = 100)
  expect_warning(V <- vcov(mod, type = "cluster", cluster = indiv, cr = "CR0"), "rank at most 3")
  expect_equal(V, Hi %*% crossprod(rowsum(U, indiv)) %*% Hi, ignore_attr = TRUE, tolerance = 1e-6)
  # splitting a track between clusters is an error
  expect_error(vcov(mod, type = "cluster", cluster = rep(1:16, each = 25)), "split between clusters")
})

test_that("vcov cluster: errors without tracks, for one cluster and bad cluster vectors", {
  mod <- lama(obj_rep) # single track
  expect_error(vcov(mod, type = "cluster"), "needs several tracks")
  expect_error(vcov(mod, type = "cluster", cluster = rep(1, 200)), "at least two clusters")
  expect_error(vcov(mod, type = "cluster", cluster = 1:10), "one value per observation")
  expect_error(vcov(mod, type = "cluster", cluster = c(NA, rep(1:2, each = 100)[-1])), "missing values")
  expect_error(vcov(mod, type = "cluster", cr = "CR2"), "should be one of")
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

test_that("MCreport: qreml() models, conditional and joint", {
  skip_if_not("par" %in% names(formals(RTMBdist::mcreport)))
  modmat <- modmat(~ s(tod, bs = "cp"), data = data.frame(tod = 1:24), knots = list(tod = c(0, 24)))
  dat <- list(step = trex$step[1:500], tod = trex$tod[1:500],
              Z = modmat$Z, S = modmat$S, lambda = rep(100, 2))
  pnll <- function(par) {
    getAll(par, dat)
    Gamma <- tpm(cbind(beta0, betaSpline), Z)
    Delta <- stationary_p(Gamma)
    allprobs <- matrix(1, length(step), 2)
    ind <- which(!is.na(step))
    for(j in 1:2) allprobs[ind, j] <- dgamma2(step[ind], exp(logmu[j]), exp(logsigma[j]))
    -forward(Delta[tod[1], ], Gamma[, , tod], allprobs) + penalty(betaSpline, S, lambda)
  }
  par <- list(logmu = log(c(0.3, 2.5)), logsigma = log(c(0.3, 2)),
              beta0 = c(-2, -2), betaSpline = matrix(0, 2, ncol(modmat$Z) - 1))
  mod <- suppressMessages(qreml(pnll, par, dat, random = "betaSpline", silent = 2))
  sample_mat <- function(s) t(sapply(seq_along(s$beta0), function(i) unlist(lapply(s, `[[`, i))))

  # conditional: centred at the final estimate with covariance J^{-1}
  set.seed(1); s <- sample_mat(suppressMessages(MCreport(mod, nSamples = 4000)))
  expect_equal(ncol(s), length(mod$par_vec))
  expect_equal(unname(colMeans(s)), unname(mod$par_vec), tolerance = 0.05)

  # joint: the same as passing Q from joint_uncertainty(), wider than conditional
  Q <- mod$joint_uncertainty()
  set.seed(1); a <- suppressMessages(MCreport(mod, nSamples = 4000, joint = TRUE))
  set.seed(1); b <- suppressMessages(MCreport(mod, nSamples = 4000, Q = Q))
  expect_equal(a, b)
  sj <- sample_mat(a)
  expect_equal(ncol(sj), length(mod$par_vec))
  V <- mod$joint_uncertainty(precision = FALSE)
  spl <- grep("betaSpline", names(mod$par_vec))
  expect_equal(unname(apply(sj[, spl], 2, var)), unname(diag(V)[spl]), tolerance = 0.15)
  expect_true(all(diag(V)[spl] >= diag(solve(mod$Hessian_conditional))[spl] - 1e-10))

  # the Hessian returned is that of the final fit, and outer_hessian() leaves no copy behind
  expect_equal(mod$Hessian_conditional, mod$obj$he(mod$par_vec), ignore_attr = TRUE)
  expect_null(environment(mod$outer_hessian)$fit_state$J)

  expect_error(MCreport(obj_rep, joint = TRUE), "qreml")
})

test_that("vcov CR3: one-step jackknife from leave-one-cluster-out Hessians", {
  mod <- lama(obj_cl)
  tp <- weighted_tape(mod)
  U <- scores(mod)
  S <- rowsum(U, dat_cl$ID, reorder = FALSE)
  D <- t(vapply(unique(dat_cl$ID), function(k) solve(obj_cl$he(opt_cl$par) - weighted_hessian(tp, dat_cl$ID == k), S[k, ]),
                numeric(6)))
  expect_equal(vcov(mod, type = "cluster", cr = "CR3"), 7 / 8 * crossprod(D), ignore_attr = TRUE, tolerance = 1e-6)
  expect_equal(rownames(vcov(mod, type = "cluster", cr = "CR3")), names(coef(mod)))
})

# linear regression as a 1-state model: Gaussian log-likelihood with known sd is quadratic,
# so the one-step leave-one-cluster-out estimates are exact
set.seed(1)
dat_lm <- list(y = rnorm(60), x = runif(60), ID = rep(c("a", "b", "c", "d"), c(10, 10, 15, 25)))
nll_lm <- function(par) {
  getAll(par, dat_lm)
  -forward(1, matrix(1, 1, 1), matrix(dnorm(y, b0 + b1 * x, 1), ncol = 1), trackID = ID)
}
obj_lm <- RTMB::MakeADFun(nll_lm, list(b0 = 0, b1 = 0), silent = TRUE)
opt_lm <- nlminb(obj_lm$par, obj_lm$fn, obj_lm$gr)

test_that("vcov CR3: exact jackknife of MacKinnon, Nielsen and Webb for linear regression", {
  mod <- lama(obj_lm)
  b <- coef(lm(y ~ x, dat_lm))
  D <- t(vapply(unique(dat_lm$ID), function(k) coef(lm(y ~ x, dat_lm, subset = ID != k)) - b, numeric(2)))
  expect_equal(vcov(mod, type = "cluster", cr = "CR3"), 3 / 4 * crossprod(D), ignore_attr = TRUE, tolerance = 1e-6)
})

test_that("vcov CR3: equals CR1 with equal leverage", {
  # mean only, equal cluster sizes: H_k = H / G, so (G - 1) / G * (G / (G - 1))^2 = G / (G - 1)
  y <- dat_lm$y; ID <- rep(1:4, each = 15)
  nll <- function(par) -forward(1, matrix(1, 1, 1), matrix(dnorm(y, par$b0, 1), ncol = 1), trackID = ID)
  obj <- RTMB::MakeADFun(nll, list(b0 = 0), silent = TRUE)
  nlminb(obj$par, obj$fn, obj$gr)
  mod <- lama(obj)
  expect_equal(vcov(mod, type = "cluster", cr = "CR3"), vcov(mod, type = "cluster", cr = "CR1"), tolerance = 1e-8)
})

test_that("vcov CR3: error for parameters informed by a single cluster", {
  # intercept per track: leaving out a track leaves its intercept unidentified
  nll <- function(par) {
    getAll(par, dat_lm)
    -forward(1, matrix(1, 1, 1), matrix(dnorm(y, b0[match(ID, unique(ID))] + b1 * x, 1), ncol = 1), trackID = ID)
  }
  obj <- RTMB::MakeADFun(nll, list(b0 = rep(0, 4), b1 = 0), silent = TRUE)
  nlminb(obj$par, obj$fn, obj$gr)
  expect_error(suppressWarnings(vcov(lama(obj), type = "cluster", cr = "CR3")), "Hessian without cluster 1 is singular or not positive definite")
})

test_that("vcov CR3: error if the Hessian without a cluster is indefinite", {
  # mock tape: H - H_k = diag(c(1, -1)) is indefinite but well conditioned
  tp <- list(p = 2, nObs = 4, names = c("a", "b"), theta = c(0, 0),
             G = function(x) c(x[1] * sum(x[3:6]), -x[2] * sum(x[3:6]), rep(0, 4)))
  U <- matrix(1, 4, 2)
  expect_error(cr3_vcov(tp, U, c(1, 1, 2, 2)), "not positive definite")
})
