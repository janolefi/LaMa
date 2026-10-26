suppressPackageStartupMessages({library(RTMB); library(LaMa)})

## exact reference: log-space forward algorithm, plain R, no tricks
## Gam_t(t, i): transition matrix into time t (row index within data), track i
ref_forward <- function(delta, Gam_t, logP, trackID = rep(1, nrow(logP))) {
  lse <- function(x) { m <- max(x); if (!is.finite(m)) return(m); m + log(sum(exp(x - m))) }
  uID <- unique(trackID)
  if (is.null(dim(delta))) delta <- matrix(delta, length(uID), length(delta), byrow = TRUE)
  l <- 0
  for (i in seq_along(uID)) {
    ind <- which(trackID == uID[i])
    la <- log(delta[i, ]) + logP[ind[1], ]
    for (t in ind[-1]) {
      G <- Gam_t(t, i)
      la <- sapply(seq_along(la), function(j) lse(la + log(G[, j]))) + logP[t, ]
    }
    l <- l + lse(la)
  }
  l
}

## evaluate forward-type function through the AD (R) path by taping
## allprobs is passed as exp(x) of a tape input so it is a genuine advector
ad_eval <- function(f, logP, logspace = FALSE, ...) {
  F <- MakeTape(function(x) {
    lp <- matrix(x, nrow(logP))
    f(allprobs = if (logspace) lp else exp(lp), logspace = logspace, report = FALSE, ...)
  }, as.vector(logP))
  F(as.vector(logP))
}
nonad_eval <- function(f, logP, logspace = FALSE, ...) {
  f(allprobs = if (logspace) logP else exp(logP), logspace = logspace, report = FALSE, ad = FALSE, ...)
}

## run expr in a fresh R process (protects against crashes / UB)
safe <- function(expr_fun) {
  r <- tryCatch(callr::r(expr_fun, error = "error", timeout = 60),
                error = function(e) paste("ERROR/CRASH:", conditionMessage(e)))
  r
}

cmp <- function(label, ...) {
  v <- list(...)
  cat(sprintf("%-55s %s\n", label,
              paste(sprintf("%s=%s", names(v), sapply(v, function(z)
                if (is.numeric(z)) format(z, digits = 10) else as.character(z))), collapse = "  ")))
}

sim_logP <- function(n, N = 2, seed = 1) {
  set.seed(seed)
  x <- rnorm(n, rep(c(0, 3), length.out = n), 1)
  mu <- seq(0, 3, length.out = N)
  sapply(1:N, function(j) dnorm(x, mu[j], 1, log = TRUE))
}
