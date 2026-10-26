source(file.path(Sys.getenv("SP"), "helpers.R"))
G1 <- matrix(c(.9, .1, .2, .8), 2, byrow = TRUE)
d0 <- c(2/3, 1/3)
err <- function(e) paste("ERR:", conditionMessage(e))
logP <- sim_logP(200)

cat("\n=== H. banded (bug 6: bw ignored in non-AD) ===\n")
ex <- ref_forward(d0, function(t, i) G1, logP)
for (bw in c(5, 10, 20, 50)) {
  a <- suppressWarnings(ad_eval(forward, logP, delta = d0, Gamma = G1, bw = bw))
  cc <- suppressWarnings(nonad_eval(forward, logP, delta = d0, Gamma = G1, bw = bw))
  cmp(paste("bw =", bw), exact = ex, ad_banded = a, cpp = cc, ad_minus_exact = a - ex)
}
a <- ad_eval(forward, logP, delta = d0, Gamma = G1, bw = 20, trackID = rep(1:4, each = 50))
cmp("bw = 20 with 4 tracks", exact = ref_forward(d0, function(t, i) G1, logP, rep(1:4, each = 50)), ad = a)
cmp("bw larger than a track (bw=60, tracks of 50)",
    exact = ref_forward(d0, function(t, i) G1, logP, rep(1:4, each = 50)),
    ad = ad_eval(forward, logP, delta = d0, Gamma = G1, bw = 60, trackID = rep(1:4, each = 50)))
cat("warning count for bw = 5 in non-AD loop of 3 calls: ")
w <- 0; withCallingHandlers(for (i in 1:3) nonad_eval(forward, logP, delta = d0, Gamma = G1, bw = 5),
                            warning = function(x) { w <<- w + 1; invokeRestart("muffleWarning") }); cat(w, "\n")

cat("\n=== I. REPORT type in forward_p (bug 7) ===\n")
set.seed(4); n <- 96; tod <- rep(1:24, 4); x <- rnorm(n)
dat <- list(x = x, tod = tod)
nllp <- function(par) {
  getAll(par, dat)
  Gamma <- tpm_p(1:24, L = 24, beta = matrix(beta, 2))
  delta <- stationary_p(Gamma, tod[1])
  ap <- cbind(dnorm(x, mu[1], 1), dnorm(x, mu[2], 1))
  -forward_p(delta, Gamma, ap, tod)
}
objp <- MakeADFun(nllp, list(mu = c(-1, 1), beta = rep(0, 6)), silent = TRUE)
rep <- objp$report()
cat("reported type:", rep$type, " | dim(Gamma) reported:", dim(rep$Gamma), " | tod reported:", !is.null(rep$tod), "\n")
cat("pseudo_res would dispatch on type =", rep$type, "\n")

cat("\n=== J. forward_p: tod not covering all L values ===\n")
tod2 <- rep(1:12, 8)  # data only covers first half of the cycle
Gp <- tpm_p(1:24, L = 24, beta = matrix(c(-1, -1, .5, .5, 0, 0), 2))
cmp("tod in 1:12, Gamma has 24 slices",
    ref = ref_forward(d0, function(t, i) Gp[, , tod2[t]], logP[1:96, ]),
    out = tryCatch(forward_p(d0, Gp, exp(logP[1:96, ]), tod2, ad = FALSE, report = FALSE), error = err))
tod3 <- rep(c(2, 5, 24), 32)
cmp("tod = c(2,5,24) repeated (3 unique values)",
    ref = ref_forward(d0, function(t, i) Gp[, , tod3[t]], logP[1:96, ]),
    out = tryCatch(forward_p(d0, Gp, exp(logP[1:96, ]), tod3, ad = FALSE, report = FALSE), error = err))

cat("\n=== K. Gamma shape classification ===\n")
lp <- logP[1:10, ]
## track-specific Gamma when nTracks == nObs - 1 -> misclassified as time-varying?
ID <- c(1, 1, 2:9)  # 9 tracks, 10 obs
Gk <- array(0, c(2, 2, 9)); set.seed(5); for (k in 1:9) Gk[, , k] <- tpm(rnorm(2, -1))
cmp("track-specific Gamma, nTracks = nObs - 1",
    ref = ref_forward(d0, function(t, i) Gk[, , i], lp, ID),
    cpp = tryCatch(nonad_eval(forward, lp, delta = d0, Gamma = Gk, trackID = ID), error = err))
## wrong-length time-varying Gamma
cmp("time-varying Gamma with nObs - 2 slices",
    out = tryCatch(nonad_eval(forward, lp, delta = d0, Gamma = array(G1, c(2, 2, 8))), error = err))
cmp("time-varying Gamma nObs - 1 slices + trackID",
    out = tryCatch(nonad_eval(forward, lp, delta = d0, Gamma = array(G1, c(2, 2, 9)), trackID = rep(1:2, each = 5)), error = err))
## transposed / non-stochastic Gamma silently accepted
cmp("Gamma rows not summing to 1 (transposed)",
    out = tryCatch(nonad_eval(forward, lp, delta = d0, Gamma = matrix(c(.9, .2, .1, .3), 2, byrow = TRUE)), error = err))
cmp("NA in allprobs",
    out = tryCatch({ap <- exp(lp); ap[3, 1] <- NA; forward(d0, G1, ap, ad = FALSE)}, error = err))
cmp("allprobs data.frame",
    out = tryCatch(forward(d0, G1, as.data.frame(exp(lp)), ad = FALSE), error = err))
cmp("delta wrong length",
    out = tryCatch(forward(c(.5, .3, .2), G1, exp(lp), ad = FALSE), error = err))

cat("\n=== L. TapeConfig messages per taping ===\n")
msgs <- 0
withCallingHandlers({
  for (i in 1:3) ad_eval(forward, lp, delta = d0, Gamma = G1)
}, message = function(m) { msgs <<- msgs + 1; invokeRestart("muffleMessage") })
cat("messages from 3 tapings:", msgs, " (TapeConfig()$matmul =", TapeConfig()$matmul, ")\n")

cat("\n=== M. gradients: AD tape vs numDeriv on exact reference ===\n")
## d/dlogP of log-likelihood, multitrack track-specific Gamma, contiguous sorted IDs
ID <- rep(1:2, each = 5); D <- rbind(c(.7, .3), c(.2, .8)); Gk2 <- Gk[, , 1:2]
F <- MakeTape(function(x) forward(D, Gk2, exp(matrix(x, 10)), trackID = ID, report = FALSE), as.vector(lp))
gnum <- numDeriv::grad(function(x) ref_forward(D, function(t, i) Gk2[, , i], matrix(x, 10), ID), as.vector(lp))
cat("max |grad AD - grad numDeriv|:", max(abs(F$jacobian(as.vector(lp)) - gnum)), "\n")
F <- MakeTape(function(x) forward(D, Gk2, matrix(x, 10), trackID = ID, logspace = TRUE, report = FALSE), as.vector(lp))
cat("logspace: max |grad AD - grad numDeriv|:", max(abs(F$jacobian(as.vector(lp)) - gnum)), "\n")
