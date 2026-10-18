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
