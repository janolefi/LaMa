# exact reference forward algorithm on the log scale (plain R, no scaling tricks)
# Gam(t, i): transition matrix into row t of logP, for track i
ref_forward <- function(delta, Gam, logP, trackID = rep(1, nrow(logP))) {
  lse <- function(x) { m <- max(x); if(!is.finite(m)) return(m); m + log(sum(exp(x - m))) }
  uID <- unique(trackID)
  if(is.null(dim(delta))) delta <- matrix(delta, length(uID), length(delta), byrow = TRUE)
  l <- 0
  for(i in seq_along(uID)) {
    ind <- which(trackID == uID[i])
    la <- log(delta[i, ]) + logP[ind[1], ]
    for(t in ind[-1]) {
      G <- Gam(t, i)
      la <- sapply(seq_along(la), function(j) lse(la + log(G[, j]))) + logP[t, ]
    }
    l <- l + lse(la)
  }
  l
}

# compare AD path (taped R code) and non-AD path (C++) of f with ref
# f(allprobs, ...) calls a forward function and passes ... on (ad = TRUE/FALSE)
expect_forward <- function(f, logP, ref, logspace = FALSE, tolerance = 1e-10) {
  F <- RTMB::MakeTape(function(x) {
    x <- matrix(x, nrow(logP))
    f(if(logspace) x else exp(x), ad = TRUE)
  }, as.vector(logP))
  expect_equal(F(as.vector(logP)), ref, tolerance = tolerance)
  expect_equal(f(if(logspace) logP else exp(logP), ad = FALSE), ref, tolerance = tolerance)
}

# compare AD gradient w.r.t. logP with central differences of ref_fun(logP)
expect_forward_grad <- function(f, logP, ref_fun, logspace = FALSE) {
  x <- as.vector(logP)
  F <- RTMB::MakeTape(function(x) {
    x <- matrix(x, nrow(logP))
    f(if(logspace) x else exp(x), ad = TRUE)
  }, x)
  h <- 1e-5
  num <- sapply(seq_along(x), function(k) {
    e <- replace(0 * x, k, h)
    (ref_fun(matrix(x + e, nrow(logP))) - ref_fun(matrix(x - e, nrow(logP)))) / (2 * h)
  })
  expect_equal(as.vector(F$jacobian(x)), num, tolerance = 1e-6)
}

# exact references for the decoding functions, same arguments as ref_forward()
# most probable state sequence by enumerating all paths of each track (small tracks only)
ref_viterbi <- function(delta, Gam, logP, trackID = rep(1, nrow(logP))) {
  uID <- unique(trackID)
  if(is.null(dim(delta))) delta <- matrix(delta, length(uID), length(delta), byrow = TRUE)
  states <- integer(0)
  for(i in seq_along(uID)) {
    ind <- which(trackID == uID[i])
    paths <- as.matrix(expand.grid(rep(list(seq_len(ncol(logP))), length(ind))))
    lp <- apply(paths, 1, function(s) {
      l <- log(delta[i, s[1]]) + logP[ind[1], s[1]]
      for(k in seq_along(ind)[-1]) l <- l + log(Gam(ind[k], i)[s[k - 1], s[k]]) + logP[ind[k], s[k]]
      l
    })
    states <- c(states, paths[which.max(lp), ])
  }
  unname(states)
}

# Pr(S_t = j | all observations): likelihood with only state j allowed at t
ref_stateprobs <- function(delta, Gam, logP, trackID = rep(1, nrow(logP))) {
  l <- ref_forward(delta, Gam, logP, trackID)
  outer(seq_len(nrow(logP)), seq_len(ncol(logP)), Vectorize(function(t, j) {
    logP[t, -j] <- -Inf
    exp(ref_forward(delta, Gam, logP, trackID) - l)
  }))
}

# Pr(S_t = j | observations before t in the track): drop the observations from t on
ref_forecast <- function(delta, Gam, logP, trackID = rep(1, nrow(logP))) {
  outer(seq_len(nrow(logP)), seq_len(ncol(logP)), Vectorize(function(t, j) {
    logP[trackID == trackID[t] & seq_len(nrow(logP)) >= t, ] <- 0
    l <- ref_forward(delta, Gam, logP, trackID)
    logP[t, -j] <- -Inf
    exp(ref_forward(delta, Gam, logP, trackID) - l)
  }))
}

# log-likelihood contributions log p(x_t | observations before t in the track):
# differences of ref_forward() with the observations from t (or t + 1) on dropped
ref_contrib <- function(delta, Gam, logP, trackID = rep(1, nrow(logP))) {
  sapply(seq_len(nrow(logP)), function(t) {
    drop_from <- function(s) replace(logP, trackID == trackID[t] & seq_len(nrow(logP)) >= s, 0)
    ref_forward(delta, Gam, drop_from(t + 1), trackID) - ref_forward(delta, Gam, drop_from(t), trackID)
  })
}

# forward() taped with weights w for its log-likelihood contributions, as set by scores()
# f(allprobs) calls a forward function, evaluated at logP
forward_weighted <- function(f, logP, w, logspace = FALSE) {
  env <- getFromNamespace(".lama_env", "LaMa")
  env$weights <- w
  env$pos <- 0
  on.exit(env$weights <- NULL)
  F <- RTMB::MakeTape(function(x) {
    x <- matrix(x, nrow(logP))
    f(if(logspace) x else exp(x))
  }, as.vector(logP))
  F(as.vector(logP))
}

# small test data: 2 states, 20 observations
set.seed(1)
logP_test <- local({
  x <- rnorm(20, rep(c(0, 3), each = 5), 1)
  cbind(dnorm(x, 0, 1, log = TRUE), dnorm(x, 3, 1, log = TRUE))
})
G1 <- matrix(c(0.9, 0.1, 0.2, 0.8), 2, byrow = TRUE)
G2 <- matrix(c(0.6, 0.4, 0.3, 0.7), 2, byrow = TRUE)
G3 <- matrix(c(0.5, 0.5, 0.1, 0.9), 2, byrow = TRUE)
Gk_test <- array(c(G1, G2, G3), dim = c(2, 2, 3)) # track-specific
Gt_test <- array(apply(matrix(seq(-2, 0, length.out = 40), 2), 2, function(e) tpm(e)), dim = c(2, 2, 20)) # time-varying
D_test <- rbind(c(0.7, 0.3), c(0.2, 0.8), c(0.5, 0.5)) # track-specific delta
d_test <- c(2/3, 1/3)
