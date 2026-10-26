source(file.path(Sys.getenv("SP"), "helpers.R"))
ids <- paste0("ID", 1:10)
cat("unique(as.integer(as.factor(ids))) =", unique(as.integer(as.factor(ids))), "\n")
## realistic fit: track-specific Gamma via sex covariate, IDs ID1..ID10 in natural order
set.seed(7); K <- 10; nT <- 200; sex <- rep(0:1, 5)
x <- unlist(lapply(1:K, function(k) { G <- tpm(c(-2.5, -2.5) + c(1.5, -1.5) * sex[k]); s <- numeric(nT); s[1] <- 1
  for (t in 2:nT) s[t] <- sample(1:2, 1, prob = G[s[t-1], ]); rnorm(nT, c(0, 3)[s], 1) }))
dat <- list(x = x, ID = rep(ids, each = nT), sex = sex)
nll <- function(par) { getAll(par, dat)
  Gk <- array(NA, c(2, 2, 10)); Gk <- AD(Gk)
  for (k in 1:10) Gk[, , k] <- tpm(beta0 + beta1 * sex[k])
  delta <- t(sapply(1:10, function(k) c(1, 0))) * 0.98 + 0.01
  ap <- cbind(dnorm(x, mu[1], 1), dnorm(x, mu[2], 1))
  -forward(delta, Gk, ap, trackID = ID) }
par <- list(mu = c(0, 3), beta0 = c(-2, -2), beta1 = c(0, 0))
obj <- MakeADFun(nll, par, silent = TRUE); opt <- nlminb(obj$par, obj$fn, obj$gr)
cat("obj$fn(opt$par) =", opt$objective, "  plain nll(opt) =", nll(relist(opt$par, par)), "\n")
cat("beta1 estimate (truth 1.5, -1.5):", round(opt$par[5:6], 2), "\n")
