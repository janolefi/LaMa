source(file.path(Sys.getenv("SP"), "helpers.R"))
N <- 2
G1 <- matrix(c(.9, .1, .2, .8), 2, byrow = TRUE)
G2 <- matrix(c(.6, .4, .3, .7), 2, byrow = TRUE)
G3 <- matrix(c(.5, .5, .1, .9), 2, byrow = TRUE)
D <- rbind(c(.7, .3), c(.2, .8), c(.5, .5))

cat("\n=== A. contiguous, sorted IDs (control) ===\n")
logP <- sim_logP(30); ID <- rep(1:3, each = 10)
Gk <- array(c(G1, G2, G3), c(2, 2, 3))
cmp("homog Gamma, track-specific delta", ref = ref_forward(D, function(t, i) G1, logP, ID),
    ad = ad_eval(forward, logP, delta = D, Gamma = G1, trackID = ID),
    cpp = nonad_eval(forward, logP, delta = D, Gamma = G1, trackID = ID))
cmp("track-specific Gamma + delta", ref = ref_forward(D, function(t, i) Gk[, , i], logP, ID),
    ad = ad_eval(forward, logP, delta = D, Gamma = Gk, trackID = ID),
    cpp = nonad_eval(forward, logP, delta = D, Gamma = Gk, trackID = ID))

cat("\n=== B. contiguous but NOT sorted IDs (new suspect: uID vs recoded trackID) ===\n")
ID <- rep(c("c", "a", "b"), each = 10) # appearance order c, a, b
cmp("char IDs c,a,b: track-specific Gamma + delta",
    ref = ref_forward(D, function(t, i) Gk[, , i], logP, ID),
    ad = ad_eval(forward, logP, delta = D, Gamma = Gk, trackID = ID),
    cpp = nonad_eval(forward, logP, delta = D, Gamma = Gk, trackID = ID))
ID <- rep(c(3, 1, 2), each = 10)
cmp("numeric IDs 3,1,2: track-specific delta only",
    ref = ref_forward(D, function(t, i) G1, logP, ID),
    ad = ad_eval(forward, logP, delta = D, Gamma = G1, trackID = ID),
    cpp = nonad_eval(forward, logP, delta = D, Gamma = G1, trackID = ID))
## forward_g, time-varying Gamma, unsorted IDs, track-specific delta
set.seed(2); Gt <- array(0, c(2, 2, 30)); for (t in 1:30) Gt[, , t] <- tpm(rnorm(2, -1))
cmp("forward_g, IDs 3,1,2, track-specific delta",
    ref = ref_forward(D, function(t, i) Gt[, , t], logP, ID),
    ad = ad_eval(forward_g, logP, delta = D, Gamma = Gt, trackID = ID),
    cpp = nonad_eval(forward_g, logP, delta = D, Gamma = Gt, trackID = ID))
cmp("forward_g, IDs 3,1,2, common delta",
    ref = ref_forward(D[1, ], function(t, i) Gt[, , t], logP, ID),
    ad = ad_eval(forward_g, logP, delta = D[1, ], Gamma = Gt, trackID = ID),
    cpp = nonad_eval(forward_g, logP, delta = D[1, ], Gamma = Gt, trackID = ID))

cat("\n=== C. interleaved (non-contiguous) IDs (bug 1) ===\n")
ID <- rep(rep(1:3, each = 5), 2) # 1..1 2..2 3..3 1..1 2..2 3..3 -> 6 runs, 3 tracks
cat("non-AD (C++, in subprocess):\n")
res <- sapply(list(
  h_common = quote(nonad_eval(forward, logP, delta = D[1, ], Gamma = G1, trackID = ID)),
  h_trackD = quote(nonad_eval(forward, logP, delta = D, Gamma = G1, trackID = ID)),
  h_trackG = quote(nonad_eval(forward, logP, delta = D, Gamma = Gk, trackID = ID)),
  h_trackD6 = quote(nonad_eval(forward, logP, delta = rbind(D, D), Gamma = G1, trackID = ID)),
  g_trackD = quote(nonad_eval(forward_g, logP, delta = D, Gamma = Gt, trackID = ID)),
  g_common = quote(nonad_eval(forward_g, logP, delta = D[1,], Gamma = Gt, trackID = ID))),
  function(ex) tryCatch(format(callr::r(function(SP, ex, ...) { source(file.path(SP, "helpers.R")); list2env(list(...), environment()); eval(ex) },
    args = list(SP = Sys.getenv("SP"), ex = ex, logP = logP, ID = ID, D = D, G1 = G1, Gk = Gk, Gt = Gt)), digits = 10),
    error = function(e) paste("ERR:", conditionMessage(e$parent %||% e))))
print(res)
cat("refs (treating interleaved runs as the same track, AD semantics):\n")
print(c(h_common = ref_forward(D[1, ], function(t, i) G1, logP, ID),
        h_trackD = ref_forward(D, function(t, i) G1, logP, ID),
        h_trackG = ref_forward(D, function(t, i) Gk[, , i], logP, ID),
        g_trackD = ref_forward(D, function(t, i) Gt[, , t], logP, ID)))
cat("AD:\n")
print(c(h_common = ad_eval(forward, logP, delta = D[1, ], Gamma = G1, trackID = ID),
        h_trackD = ad_eval(forward, logP, delta = D, Gamma = G1, trackID = ID),
        h_trackG = ad_eval(forward, logP, delta = D, Gamma = Gk, trackID = ID),
        g_trackD = ad_eval(forward_g, logP, delta = D, Gamma = Gt, trackID = ID)))
cat("refs treating each run as its own track (C++ run semantics, common delta):\n")
runID <- rep(seq_along(rle(ID)$lengths), rle(ID)$lengths)
print(c(h_common = ref_forward(D[1, ], function(t, i) G1, logP, runID)))

cat("\n=== D. length-1 tracks (bug 2) ===\n")
ID <- c(1, 2, 2, 2, 2, 3, 4, 4, 4, 4) ; lp <- logP[1:10, ]
for (f in c("forward", "forward_g")) {
  Gam <- if (f == "forward") G1 else Gt[, , 1:10]
  r <- ref_forward(D[1, ], if (f == "forward") function(t, i) G1 else function(t, i) Gt[, , t], lp, ID)
  a <- tryCatch(ad_eval(get(f), lp, delta = D[1, ], Gamma = Gam, trackID = ID), error = function(e) paste("ERR:", conditionMessage(e)))
  cc <- tryCatch(nonad_eval(get(f), lp, delta = D[1, ], Gamma = Gam, trackID = ID), error = function(e) paste("ERR:", conditionMessage(e)))
  cmp(paste(f, "with length-1 tracks"), ref = r, ad = a, cpp = cc)
}
## single observation, single track
a <- tryCatch(ad_eval(forward, logP[1, , drop = FALSE], delta = D[1, ], Gamma = G1), error = function(e) paste("ERR:", conditionMessage(e)))
cc <- tryCatch(nonad_eval(forward, logP[1, , drop = FALSE], delta = D[1, ], Gamma = G1), error = function(e) paste("ERR:", conditionMessage(e)))
cmp("forward, nObs = 1", ref = ref_forward(D[1, ], function(t, i) G1, logP[1, , drop = FALSE]), ad = a, cpp = cc)
## nObs = 2 single track, forward_g with nObs slices (Gamma[,,-1] drops to matrix?)
lp2 <- logP[1:2, ]
a <- tryCatch(ad_eval(forward_g, lp2, delta = D[1, ], Gamma = Gt[, , 1:2]), error = function(e) paste("ERR:", conditionMessage(e)))
cc <- tryCatch(nonad_eval(forward_g, lp2, delta = D[1, ], Gamma = Gt[, , 1:2]), error = function(e) paste("ERR:", conditionMessage(e)))
cmp("forward_g, nObs = 2, 2 slices", ref = ref_forward(D[1, ], function(t, i) Gt[, , t], lp2), ad = a, cpp = cc)
