source(file.path(Sys.getenv("SP"), "helpers.R"))
source(file.path(Sys.getenv("SP"), "proto.R"))
G1 <- matrix(c(.9, .1, .2, .8), 2, byrow = TRUE)
d0 <- c(2/3, 1/3)
n <- 500
logP <- sim_logP(n)
ntape <- function(F) nrow(F$data.frame())
tape <- function(fun) MakeTape(function(x) fun(matrix(x, n)), as.vector(logP))

cat("\n=== P1. row-max shift: exactness of value and gradient ===\n")
D1 <- matrix(d0, 1); G3 <- array(G1, c(2, 2, 1))
ex <- ref_forward(d0, function(t, i) G1, logP)
Fs <- tape(function(lp) proto_shift(D1, G3, lp, rep(1, n), 1, n))
gnum <- numDeriv::grad(function(x) ref_forward(d0, function(t, i) G1, matrix(x, n)), as.vector(logP))
cat("value diff:", Fs(as.vector(logP)) - ex, " max grad diff:", max(abs(Fs$jacobian(as.vector(logP)) - gnum)), "\n")
## gradient wrt a parameter entering through allprobs AND max (ties too)
lpt <- logP; lpt[5, ] <- c(-1, -1) # exact tie in row 5
Ft <- tape(function(lp) proto_shift(D1, G3, lp, rep(1, n), 1, n))
gnum_t <- numDeriv::grad(function(x) ref_forward(d0, function(t, i) G1, matrix(x, n)), as.vector(lpt))
cat("with tie: max grad diff:", max(abs(Ft$jacobian(as.vector(lpt)) - gnum_t)), "\n")

cat("\n=== P2. row-max shift: stability cases ===\n")
stab <- function(label, lp, delta = d0, G = G1) {
  r <- ref_forward(delta, function(t, i) G, lp)
  Fs <- MakeTape(function(x) proto_shift(matrix(delta, 1), array(G, c(2, 2, 1)), matrix(x, nrow(lp)), rep(1, nrow(lp)), 1, nrow(lp)), as.vector(lp))
  old <- ad_eval(forward, lp, logspace = TRUE, delta = delta, Gamma = G)
  cmp(label, ref = r, shift = Fs(as.vector(lp)), old_logspace = old)
}
stab("logP - 800", logP - 800)
stab("logP - 1e5", logP - 1e5)
set.seed(3); s <- rep(1:2, each = 50); x <- rnorm(100, c(0, 50)[s], 0.5)
stab("tight separated states", cbind(dnorm(x, 0, .5, log = TRUE), dnorm(x, 50, .5, log = TRUE)))
lp <- logP; lp[1, ] <- c(-1000, 0)
stab("delta=c(1,0), logP[1,]=c(-1000,0)", lp, delta = c(1, 0))
## structural zero in Gamma: left-to-right chain, observation conflicts strongly with forced state
Gz <- matrix(c(.9, .1, 0, 1), 2, byrow = TRUE)
lp <- logP[1:50, ]; lp[30:50, ] <- cbind(rep(0, 21), rep(-1000, 21)); lp[29, ] <- c(-1000, 0)
stab("left-to-right Gamma, data says state 1 after absorption", lp, delta = c(1, 0), G = Gz)
## near-zero (not structural) transition probability
Gs <- tpm(c(-800, -2)) # off-diagonal 1->2 ~ exp(-800) underflows to 0 in double
cat("tpm(c(-800,-2))[1,2] =", Gs[1, 2], "\n")

cat("\n=== P3. tape sizes (n =", n, ", N = 2) ===\n")
trackS <- 1; trackE <- n
F_cur_h <- tape(function(lp) forward(d0, G1, exp(lp), report = FALSE))
F_new_h <- tape(function(lp) proto_kernel(D1, G3, exp(lp), rep(1, n), 1, n))
Gt <- array(G1, c(2, 2, n))
F_cur_g <- MakeTape(function(x) forward_g(d0, AD(Gt) * x[1] / x[1], exp(logP), report = FALSE), 1)
F_new_g <- MakeTape(function(x) proto_kernel(D1, AD(Gt) * x[1] / x[1], exp(logP), 1:n, 1, n), 1)
Gp <- tpm_p(1:24, L = 24, beta = matrix(c(-1, -1, .5, .5, 0, 0), 2)); tod <- rep(1:24, length.out = n)
F_cur_p <- MakeTape(function(x) forward_p(d0, AD(Gp) * x[1] / x[1], exp(logP), tod, report = FALSE), 1)
F_new_p <- MakeTape(function(x) proto_kernel(D1, AD(Gp) * x[1] / x[1], exp(logP), tod, 1, n), 1)
F_cur_log <- tape(function(lp) forward(d0, G1, lp, logspace = TRUE, report = FALSE))
F_new_log <- tape(function(lp) proto_shift(D1, G3, lp, rep(1, n), 1, n))
## AD Gamma tapes: x*Gamma/x adds n*N^2*2 ops on input side; subtract baseline
base_g <- ntape(MakeTape(function(x) sum(AD(Gt) * x[1] / x[1]), 1)) ; base_p <- ntape(MakeTape(function(x) sum(AD(Gp) * x[1] / x[1]), 1))
print(data.frame(
  case = c("homogeneous", "inhomog (nObs slices)", "periodic (L=24)", "logspace"),
  current = c(ntape(F_cur_h), ntape(F_cur_g), ntape(F_cur_p), ntape(F_cur_log)),
  proto   = c(ntape(F_new_h), ntape(F_new_g), ntape(F_new_p), ntape(F_new_log)),
  check_equal = c(F_cur_h(as.vector(logP)) - F_new_h(as.vector(logP)), F_cur_g(1) - F_new_g(1),
                  F_cur_p(1) - F_new_p(1), F_cur_log(as.vector(logP)) - F_new_log(as.vector(logP)))))
cat("(periodic/inhomog include", base_g, "/", base_p, "ops for building Gamma inputs)\n")

cat("\n=== P4. taping time, n = 20000, 4 tracks, track-specific delta ===\n")
n2 <- 20000; lp2 <- sim_logP(n2); ID <- rep(1:4, each = n2 / 4); D4 <- matrix(d0, 4, 2, byrow = TRUE)
t_cur <- system.time(MakeTape(function(x) forward(D4, G1, exp(matrix(x, n2)), trackID = ID, report = FALSE), as.vector(lp2)))[3]
t_new <- system.time(MakeTape(function(x) proto_kernel(D4, G3, exp(matrix(x, n2)), rep(1, n2), c(1, 5001, 10001, 15001), c(5000, 10000, 15000, 20000)), as.vector(lp2)))[3]
cat("current:", t_cur, "s   proto:", t_new, "s\n")
