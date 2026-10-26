source(file.path(Sys.getenv("SP"), "helpers.R"))
G1 <- matrix(c(.9, .1, .2, .8), 2, byrow = TRUE)
d0 <- c(2/3, 1/3)
err <- function(e) paste("ERR:", conditionMessage(e))

cat("\n=== E. logspace = TRUE, no underflow (bug 5: delta + 1e-7) ===\n")
logP <- sim_logP(200)
cmp("moderate logP, stationary-ish delta",
    ref = ref_forward(d0, function(t, i) G1, logP),
    ad_log = ad_eval(forward, logP, logspace = TRUE, delta = d0, Gamma = G1),
    ad_lin = ad_eval(forward, logP, delta = d0, Gamma = G1),
    cpp_log = nonad_eval(forward, logP, logspace = TRUE, delta = d0, Gamma = G1))
cat("diff ad_log - ref:", ad_eval(forward, logP, logspace = TRUE, delta = d0, Gamma = G1) -
      ref_forward(d0, function(t, i) G1, logP), "\n")
## delta with a zero entry
dz <- c(1, 0)
cmp("delta = c(1, 0)",
    ref = ref_forward(dz, function(t, i) G1, logP),
    ad_log = ad_eval(forward, logP, logspace = TRUE, delta = dz, Gamma = G1),
    ad_lin = ad_eval(forward, logP, delta = dz, Gamma = G1))
## delta = c(1,0) and first obs impossible-ish under state 1
lp <- logP; lp[1, ] <- c(-1000, 0)
cmp("delta = c(1,0), logP[1,] = c(-1000, 0)",
    ref = ref_forward(dz, function(t, i) G1, lp),
    ad_log = ad_eval(forward, lp, logspace = TRUE, delta = dz, Gamma = G1),
    ad_lin = ad_eval(forward, lp, delta = dz, Gamma = G1))

cat("\n=== F. logspace = TRUE with underflow (bug 4) ===\n")
lp <- logP - 800 # every density underflows exp()
cmp("logP - 800",
    ref = ref_forward(d0, function(t, i) G1, lp),
    ad_log = ad_eval(forward, lp, logspace = TRUE, delta = d0, Gamma = G1),
    cpp_log = nonad_eval(forward, lp, logspace = TRUE, delta = d0, Gamma = G1))
## a realistic case: within-row range large (well separated tight states)
set.seed(3); s <- rep(1:2, each = 50); x <- rnorm(100, c(0, 50)[s], 0.5)
lp2 <- cbind(dnorm(x, 0, .5, log = TRUE), dnorm(x, 50, .5, log = TRUE))
cat("range of logP rows:", range(apply(lp2, 1, function(r) diff(range(r)))), "\n")
cmp("tight separated states (row range ~5000)",
    ref = ref_forward(d0, function(t, i) G1, lp2),
    ad_log = ad_eval(forward, lp2, logspace = TRUE, delta = d0, Gamma = G1),
    ad_lin = ad_eval(forward, lp2, delta = d0, Gamma = G1),
    cpp_log = nonad_eval(forward, lp2, logspace = TRUE, delta = d0, Gamma = G1))
## forward_g logspace too
Gt <- array(G1, c(2, 2, 200))
cmp("forward_g, logP - 800",
    ref = ref_forward(d0, function(t, i) G1, lp),
    ad_log = ad_eval(forward_g, lp, logspace = TRUE, delta = d0, Gamma = Gt),
    cpp_log = nonad_eval(forward_g, lp, logspace = TRUE, delta = d0, Gamma = Gt))
## user-level view: nll(par) vs obj$fn(par)
dat <- list(x = x)
nll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta); delta <- stationary(Gamma)
  lap <- cbind(dnorm(x, mu[1], exp(ls[1]), log = TRUE), dnorm(x, mu[2], exp(ls[2]), log = TRUE))
  -forward(delta, Gamma, lap, logspace = TRUE)
}
par <- list(mu = c(0, 50), ls = log(c(.5, .5)), eta = c(-2, -2))
obj <- MakeADFun(nll, par, silent = TRUE)
cmp("user nll(par) vs obj$fn(par), logspace", plain_nll = nll(par), obj_fn = obj$fn(obj$par))

cat("\n=== G. one state (bug 3) ===\n")
lp1 <- logP[, 1, drop = FALSE]
cmp("nStates = 1, logspace = FALSE", ref = sum(lp1),
    ad = ad_eval(forward, lp1, delta = 1, Gamma = matrix(1)),
    cpp = nonad_eval(forward, lp1, delta = 1, Gamma = matrix(1)))
cmp("nStates = 1, logspace = TRUE", ref = sum(lp1),
    ad = suppressWarnings(ad_eval(forward, lp1, logspace = TRUE, delta = 1, Gamma = matrix(1))),
    cpp = suppressWarnings(nonad_eval(forward, lp1, logspace = TRUE, delta = 1, Gamma = matrix(1))))
cmp("forward_g, nStates = 1, logspace = TRUE", ref = sum(lp1),
    ad = suppressWarnings(ad_eval(forward_g, lp1, logspace = TRUE, delta = 1, Gamma = array(1, c(1, 1, 200)))))
nll1 <- function(par) { getAll(par); -forward(1, matrix(1), matrix(dnorm(dat$x, mu, 1), ncol = 1)) }
obj1 <- MakeADFun(nll1, list(mu = 0), silent = TRUE)
cat("REPORT names for one-state model:", paste(names(obj1$report()), collapse = ","), "(empty = bypassed)\n")
## bad inputs silently accepted by shortcut
cmp("nStates = 1, wrong-length trackID",
    out = tryCatch(nonad_eval(forward, lp1, delta = 1, Gamma = matrix(1), trackID = 1:3), error = err))
