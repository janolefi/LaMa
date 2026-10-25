lp <- logP_test

test_that("forward: single track, homogeneous", {
  ref <- ref_forward(d_test, function(t, i) G1, lp)
  expect_forward(function(ap, ...) forward(d_test, G1, ap, report = FALSE, ...), lp, ref)
})

test_that("forward: multiple tracks, common and track-specific delta and Gamma", {
  ID <- rep(1:3, c(5, 8, 7))
  expect_forward(function(ap, ...) forward(d_test, G1, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(d_test, function(t, i) G1, lp, ID))
  expect_forward(function(ap, ...) forward(D_test, G1, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(D_test, function(t, i) G1, lp, ID))
  expect_forward(function(ap, ...) forward(D_test, Gk_test, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(D_test, function(t, i) Gk_test[, , i], lp, ID))
})

test_that("forward: track-specific inputs follow the order of appearance of IDs", {
  ref <- ref_forward(D_test, function(t, i) Gk_test[, , i], lp, rep(1:3, c(5, 8, 7)))
  for(ID in list(rep(c("c", "a", "b"), c(5, 8, 7)),
                 rep(c(3, 1, 2), c(5, 8, 7)),
                 rep(c("ID2", "ID10", "ID9"), c(5, 8, 7)))) {
    expect_forward(function(ap, ...) forward(D_test, Gk_test, ap, trackID = ID, report = FALSE, ...), lp, ref)
  }
})

test_that("forward: tracks of length one and nObs = 1", {
  ID <- rep(1:5, c(1, 7, 1, 10, 1))
  expect_forward(function(ap, ...) forward(d_test, G1, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(d_test, function(t, i) G1, lp, ID))
  expect_forward(function(ap, ...) forward(d_test, Gt_test, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(d_test, function(t, i) Gt_test[, , t], lp, ID))
  lp1 <- lp[1, , drop = FALSE]
  expect_forward(function(ap, ...) forward(d_test, G1, ap, report = FALSE, ...), lp1,
                 ref_forward(d_test, function(t, i) G1, lp1))
})

test_that("forward: track-specific Gamma with nTracks = nObs - 1", {
  lp10 <- lp[1:10, ]
  ID <- c(1, 1, 2:9)
  Gk <- array(G1, dim = c(2, 2, 9)); Gk[, , 1] <- G2
  expect_forward(function(ap, ...) forward(d_test, Gk, ap, trackID = ID, report = FALSE, ...), lp10,
                 ref_forward(d_test, function(t, i) Gk[, , i], lp10, ID))
})

test_that("forward: time-varying Gamma with nObs and nObs - 1 slices", {
  ref <- ref_forward(d_test, function(t, i) Gt_test[, , t], lp)
  expect_forward(function(ap, ...) forward(d_test, Gt_test, ap, report = FALSE, ...), lp, ref)
  expect_forward(function(ap, ...) forward(d_test, Gt_test[, , -1], ap, report = FALSE, ...), lp, ref)
  expect_forward(function(ap, ...) forward(d_test, Gt_test, ap, report = FALSE, ...), lp, ref)
  # nObs = 2
  lp2 <- lp[1:2, ]
  ref2 <- ref_forward(d_test, function(t, i) Gt_test[, , t], lp2)
  expect_forward(function(ap, ...) forward(d_test, Gt_test[, , 1:2], ap, report = FALSE, ...), lp2, ref2)
  expect_forward(function(ap, ...) forward(d_test, Gt_test[, , 2, drop = FALSE], ap, report = FALSE, ...), lp2, ref2)
  # multiple tracks
  ID <- rep(c("b", "a", "c"), c(5, 8, 7))
  expect_forward(function(ap, ...) forward(D_test, Gt_test, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(D_test, function(t, i) Gt_test[, , t], lp, ID))
})

test_that("forward: periodic Gamma with tod, also if tod does not cover the full cycle", {
  Gp <- tpm_p(1:24, L = 24, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2))
  for(tod in list(rep(1:4, 5), c(3:24, 1:2)[1:20], rep(c(2, 5, 24), length.out = 20))) {
    expect_forward(function(ap, ...) forward(d_test, Gp, ap, tod = tod, report = FALSE, ...), lp,
                   ref_forward(d_test, function(t, i) Gp[, , tod[t]], lp))
  }
  ID <- rep(1:3, c(5, 8, 7)); tod <- rep(1:4, 5)
  expect_forward(function(ap, ...) forward(D_test, Gp, ap, tod = tod, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(D_test, function(t, i) Gp[, , tod[t]], lp, ID))
})

test_that('forward: delta = "stationary" for all types of Gamma', {
  ID <- rep(1:3, c(5, 8, 7)); s <- c(1, 6, 14) # track starts
  expect_forward(function(ap, ...) forward("stationary", G1, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(d_test, function(t, i) G1, lp, ID))
  expect_forward(function(ap, ...) forward("stationary", Gk_test, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(stationary(Gk_test), function(t, i) Gk_test[, , i], lp, ID))
  # time-varying: stationary distribution of the first slice of each track
  expect_forward(function(ap, ...) forward("stationary", Gt_test, ap, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(stationary(Gt_test[, , s]), function(t, i) Gt_test[, , t], lp, ID))
  expect_forward(function(ap, ...) forward("stationary", Gt_test[, , -1], ap, report = FALSE, ...), lp,
                 ref_forward(stationary(Gt_test[, , 2]), function(t, i) Gt_test[, , t], lp))
  # periodic: Gp[,,s] leads into s, so delta at s is stationary for the cycle Gp[,,s+1] ... Gp[,,s]
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2))
  cyc <- function(s) Reduce(`%*%`, lapply((s + 0:3) %% 4 + 1, function(k) Gp[, , k]))
  tod <- rep(c(2:4, 1), 5)
  expect_forward(function(ap, ...) forward("stationary", Gp, ap, tod = tod, trackID = ID, report = FALSE, ...), lp,
                 ref_forward(t(sapply(tod[s], function(k) stationary(cyc(k)))), function(t, i) Gp[, , tod[t]], lp, ID))
  expect_error(forward("stat", G1, exp(lp), ad = FALSE), "stationary")
  expect_error(forward("stationary", 0.9 * G1, exp(lp)), "sub-stochastic")
  expect_error(forward("stationary", Gt_test * 0.99, exp(lp)), "sub-stochastic")
})

test_that("forward: one state, with and without logspace", {
  lp1 <- lp[, 1, drop = FALSE]
  expect_forward(function(ap, ...) forward(1, matrix(1), ap, report = FALSE, ...), lp1, sum(lp1))
  expect_forward(function(ap, ...) forward(1, matrix(1), ap, logspace = TRUE, report = FALSE, ...), lp1,
                 sum(lp1), logspace = TRUE)
  expect_forward(function(ap, ...) forward(1, array(1, c(1, 1, 20)), ap, logspace = TRUE, report = FALSE, ...), lp1,
                 sum(lp1), logspace = TRUE)
})

test_that("forward: logspace is exact and stable", {
  # delta must not be modified
  ref <- ref_forward(c(1, 0), function(t, i) G1, lp)
  expect_forward(function(ap, ...) forward(c(1, 0), G1, ap, logspace = TRUE, report = FALSE, ...), lp,
                 ref, logspace = TRUE, tolerance = 1e-12)
  # every density underflows exp()
  lpu <- lp - 800
  expect_forward(function(ap, ...) forward(d_test, G1, ap, logspace = TRUE, report = FALSE, ...), lpu,
                 ref_forward(d_test, function(t, i) G1, lpu), logspace = TRUE)
  expect_forward(function(ap, ...) forward(d_test, Gt_test, ap, logspace = TRUE, report = FALSE, ...), lpu,
                 ref_forward(d_test, function(t, i) Gt_test[, , t], lpu), logspace = TRUE)
  # very different densities in one row
  lpr <- lp; lpr[, 2] <- lpr[, 2] - 3000
  expect_forward(function(ap, ...) forward(d_test, G1, ap, logspace = TRUE, report = FALSE, ...), lpr,
                 ref_forward(d_test, function(t, i) G1, lpr), logspace = TRUE)
})

test_that("forward: banded approximation is close to exact, non-AD is exact", {
  set.seed(2)
  x <- rnorm(200, rep(c(0, 3), each = 10), 1)
  lpb <- cbind(dnorm(x, 0, 1, log = TRUE), dnorm(x, 3, 1, log = TRUE))
  ID <- rep(1:4, each = 50)
  ref <- ref_forward(d_test, function(t, i) G1, lpb, ID)
  F <- RTMB::MakeTape(function(x) forward(d_test, G1, exp(matrix(x, 200)), trackID = ID, bw = 15, report = FALSE),
                      as.vector(lpb))
  expect_equal(F(as.vector(lpb)), ref, tolerance = 1e-8)
  expect_false(F(as.vector(lpb)) == ref)
  expect_equal(forward(d_test, G1, exp(lpb), trackID = ID, bw = 15, ad = FALSE), ref, tolerance = 1e-10)
  # bandwidth larger than tracks: exact
  expect_forward(function(ap, ...) forward(d_test, G1, ap, trackID = ID, bw = 60, report = FALSE, ...), lpb, ref)
})

test_that("forward: gradients match finite differences", {
  ID <- rep(1:3, c(5, 8, 7)); tod <- rep(1:4, 5)
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2))
  expect_forward_grad(function(ap, ...) forward(D_test, Gk_test, ap, trackID = ID, report = FALSE, ...), lp,
                      function(x) ref_forward(D_test, function(t, i) Gk_test[, , i], x, ID))
  expect_forward_grad(function(ap, ...) forward(d_test, Gt_test, ap, trackID = ID, report = FALSE, ...), lp,
                      function(x) ref_forward(d_test, function(t, i) Gt_test[, , t], x, ID))
  expect_forward_grad(function(ap, ...) forward(d_test, Gp, ap, tod = tod, report = FALSE, ...), lp,
                      function(x) ref_forward(d_test, function(t, i) Gp[, , tod[t]], x))
  expect_forward_grad(function(ap, ...) forward(d_test, G1, ap, logspace = TRUE, report = FALSE, ...), lp,
                      function(x) ref_forward(d_test, function(t, i) G1, x), logspace = TRUE)
})

test_that("forward: informative errors for bad input", {
  ap <- exp(lp)
  ID <- rep(rep(1:2, each = 5), 2) # interleaved tracks
  expect_error(forward(d_test, G1, ap, trackID = ID, ad = FALSE), "contiguous")
  expect_error(RTMB::MakeTape(function(x) forward(d_test, G1, ap * x, trackID = ID, report = FALSE), 1), "contiguous")
  expect_error(forward(d_test, array(G1, c(2, 2, 18)), ap, ad = FALSE), "third dimension")
  expect_error(forward(d_test, array(G1, c(2, 2, 19)), ap, trackID = rep(1:2, each = 10), ad = FALSE), "third dimension")
  expect_error(forward(d_test, G1, ap, trackID = 1:3, ad = FALSE), "trackID")
  expect_error(forward(1, matrix(1), ap[, 1, drop = FALSE], trackID = 1:3, ad = FALSE), "trackID")
  expect_error(forward(c(0.5, 0.3, 0.2), G1, ap, ad = FALSE), "delta")
  expect_error(forward(d_test, array(G1, c(2, 2, 3)), ap, tod = rep(1:4, 5), ad = FALSE), "tod")
})

test_that("forward: value checks in the non-AD version", {
  ap <- exp(lp); ap[3, 1] <- NA
  expect_error(forward(d_test, G1, ap, ad = FALSE), "NA")
  expect_error(forward(d_test, G1, lp, ad = FALSE), "logspace") # log-densities without logspace = TRUE
  expect_no_warning(forward(d_test, 0.9 * G1, exp(lp), ad = FALSE)) # sub-stochastic Gamma is allowed (MMPPs)
})

test_that("forward: reported quantities", {
  dat <- list(x = rnorm(20), ID = rep(c("b", "a"), each = 10))
  fit <- function(f) {
    nll <- function(par) {
      getAll(par, dat)
      ap <- cbind(dnorm(x, mu[1], 1), dnorm(x, mu[2], 1))
      -f(ap)
    }
    RTMB::MakeADFun(nll, list(mu = c(-1, 1)), silent = TRUE)$report()
  }
  r <- fit(function(ap) forward(d_test, G1, ap, trackID = dat$ID))
  expect_equal(r$type, "homogeneous")
  expect_equal(r$trackID, rep(1:2, each = 10))
  expect_equal(r$Gamma, G1)
  expect_equal(r$delta, d_test)
  r <- fit(function(ap) forward("stationary", G1, ap, trackID = dat$ID))
  expect_equal(r$delta, rbind(d_test, d_test), ignore_attr = TRUE)
  r <- fit(function(ap) forward(d_test, array(G1, c(2, 2, 2)), ap, trackID = dat$ID))
  expect_equal(r$type, "track-specific")
  r <- fit(function(ap) forward(d_test, Gt_test, ap))
  expect_equal(r$type, "inhomogeneous")
  Gp <- tpm_p(1:4, L = 4, beta = matrix(c(-1, -1, 0.5, 0.5, 0, 0), 2))
  r <- fit(function(ap) forward(d_test, Gp, ap, tod = rep(1:4, 5)))
  expect_equal(r$type, "periodic")
  expect_equal(r$Gamma, Gp)
  expect_equal(r$tod, rep(1:4, 5))
  r <- fit(function(ap) forward(1, matrix(1), ap[, 1, drop = FALSE], logspace = FALSE))
  expect_equal(r$type, "homogeneous")
})

test_that("forward: TapeConfig tip is shown only once per session", {
  old <- RTMB::TapeConfig()$matmul
  on.exit(RTMB::TapeConfig(matmul = old))
  RTMB::TapeConfig(matmul = "atomic")
  assign("tape_msg_shown", FALSE, envir = getFromNamespace(".lama_env", "LaMa"))
  tape <- function() RTMB::MakeTape(function(x) forward(d_test, G1, exp(lp) * x, report = FALSE), 1)
  expect_message(tape(), "TapeConfig")
  expect_no_message(tape())
})
