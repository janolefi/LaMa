## Prototype: pseudo_res(mod) without details, by rerunning the likelihood function
## with the densities swapped for CDFs. Nothing in the package is changed: the
## likelihood function is taken from obj and its densities and forward() are
## shadowed from outside.
##
## Idea: in the run for data stream s
##   - densities of s return their CDF F_j(x_t), all other densities return 1
##   - forward() records the "allprobs" matrix it receives
## so the user's own code puts F_j(x_t) into exactly the cells of allprobs.
## Then u_t = sum_j Pr(S_t = j | x_1, ..., x_t-1) F_j(x_t) and residual = qnorm(u_t).
##
## Check: the swap is valid iff each cell of allprobs is a combination of the
## stream's densities with weights summing to one (single density, products with
## other streams, mixtures within a state). Hence, with all densities set to 1 every
## cell must be 1, and with the densities of s set to 2 (others 1) every cell must
## be 1 (s not observed) or 2 (s observed).

pkgload::load_all(".") # LaMa from this repo (report() stores obj)

# discrete families, as in pseudo_res()
discrete_fams <- c("pois","binom","geom","nbinom","nbinom2", "betabinom","genpois","zibinom",
                   "zibinom2", "zipois","ztbinom","ztnbinom","ztnbinom2", "ztpois")

# name of the data stream of a density call, e.g. step from step[ind] or dat$step
stream_name <- function(e) {
  if(is.call(e) && identical(e[[1]], as.name("$"))) return(as.character(e[[3]]))
  if(is.call(e) && identical(e[[1]], as.name("["))) return(stream_name(e[[2]]))
  all.vars(e)[1]
}

# evaluates the likelihood function of obj at the estimates, with
# - each density replaced according to mode(stream): "one", "two", "cdf", "cdf_lower"
# - forward() recording the allprobs matrix it receives
run_mode <- function(obj, dens, mode) {
  f <- environment(obj$report)$func
  fenv <- environment(f)
  rec <- new.env()
  rec$P <- NULL
  rec$fams <- list() # families of each stream
  shadow <- new.env(parent = fenv)
  for(d in dens) local({
    dfun <- get(d, envir = fenv)
    pfun <- get(sub("^d", "p", d), envir = fenv)
    assign(d, function(x, ...) {
      s <- stream_name(substitute(x))
      rec$fams[[s]] <- union(rec$fams[[s]], sub("^d", "", d))
      args <- list(...)
      lg <- isTRUE(args$log)
      args$log <- NULL
      v <- dfun(x, ...) # only for the shape of the result
      fill <- function(val) { r <- rep(val, length(v)); dim(r) <- dim(v); r }
      out <- switch(mode(s),
        one = fill(1),
        two = fill(2),
        cdf = do.call(pfun, c(list(x), args)),
        cdf_lower = do.call(pfun, c(list(x - 1), args)))
      if(lg) log(out) else out
    }, envir = shadow)
  })
  rec_forward <- function(delta, Gamma, allprobs, ..., logspace = FALSE) {
    if(!is.null(rec$P)) stop("Only one call of forward() is supported for now.")
    rec$P <- if(logspace) exp(allprobs) else allprobs
    0
  }
  for(fw in c("forward", "forward_g", "forward_p")) assign(fw, rec_forward, envir = shadow)
  environment(f) <- shadow
  pl <- obj$env$parList(obj$env$last.par.best)
  attr(pl, "check.passed") <- NULL
  f(pl)
  list(P = as.matrix(rec$P), fams = rec$fams)
}

pseudo_res_auto <- function(mod, randomise = TRUE) {
  obj <- get_obj(mod)
  f <- environment(obj$report)$func
  # densities: functions d* called in the likelihood function with a CDF p*
  fns <- codetools::findGlobals(f, merge = FALSE)$functions
  dens <- Filter(function(d) exists(sub("^d", "p", d), envir = environment(f), mode = "function"),
                 grep("^d", fns, value = TRUE))
  if(!length(dens)) stop("No density with a CDF found in the likelihood function.")

  # check 1: all densities 1 -> every cell 1
  one <- run_mode(obj, dens, function(s) "one")
  bad <- abs(one$P - 1) > 1e-10
  if(any(bad)) stop("With all densities set to 1, allprobs is not 1 in ", sum(bad), " cells, ",
                    "e.g. row ", which(bad, arr.ind = TRUE)[1, 1], ". The densities do not enter ",
                    "allprobs as a plain product/mixture, so CDFs cannot be swapped in. ",
                    "Use pseudo_res(obs, dist, par, mod = mod) instead.")

  fc <- stateprobs(mod, forecast = TRUE)
  res <- list()
  for(s in names(one$fams)) {
    # check 2: densities of s set to 2 -> cells 1 (s missing) or 2 (s observed)
    two <- run_mode(obj, dens, function(x) if(x == s) "two" else "one")$P
    if(all(abs(two - 1) < 1e-10)) next # stream does not enter allprobs, e.g. a prior
    ok <- abs(two - 1) < 1e-10 | abs(two - 2) < 1e-10
    if(any(!ok)) stop("Stream '", s, "' enters allprobs non-linearly, e.g. row ",
                      which(!ok, arr.ind = TRUE)[1, 1], ". Use pseudo_res(obs, dist, par, mod = mod) instead.")
    observed <- rowSums(abs(two - 2) < 1e-10) == ncol(two)
    if(any(rowSums(abs(two - 2) < 1e-10) %% ncol(two) != 0))
      stop("Stream '", s, "' is observed in some states only.")

    # predictive CDF
    u <- rep(NA_real_, nrow(fc))
    upper <- rowSums(fc * run_mode(obj, dens, function(x) if(x == s) "cdf" else "one")$P)
    if(any(one$fams[[s]] %in% discrete_fams)) {
      lower <- rowSums(fc * run_mode(obj, dens, function(x) if(x == s) "cdf_lower" else "one")$P)
      u[observed] <- if(randomise) runif(sum(observed), lower[observed], upper[observed]) else
        (lower[observed] + upper[observed]) / 2
    } else {
      u[observed] <- upper[observed]
    }
    res[[s]] <- qnorm(u)
  }
  if(length(res) == 1) res[[1]] else res
}


## ---------------------------------------------------------------------------
## Good cases: compare with residuals built by hand
## ---------------------------------------------------------------------------

## 1) loop over states, covariate-dependent mean, missing values, prior on eta
dat <- list(step = trex$step[1:500], z = rnorm(500))
dat$step[c(10, 50)] <- NA
par <- list(log_mu = log(c(0.3, 2.5)), b = c(0, 0), log_sigma = log(c(0.2, 1.5)), eta = rep(-2, 2))
nll <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta)
  delta <- stationary(Gamma)
  sigma <- exp(log_sigma)
  allprobs <- matrix(1, length(step), 2)
  ind <- which(!is.na(step))
  for(j in 1:2) {
    mu <- exp(log_mu[j] + b[j] * z[ind])
    allprobs[ind, j] <- dgamma2(step[ind], mu, sigma[j])
  }
  -forward(delta, Gamma, allprobs) - sum(dnorm(eta, -2, 2, log = TRUE)) # prior: not in allprobs
}
obj <- MakeADFun(nll, par, silent = TRUE)
opt <- nlminb(obj$par, obj$fn, obj$gr)
mod <- report(obj)
r_auto <- pseudo_res_auto(mod)
mu_t <- sapply(1:2, function(j) exp(mod$par$log_mu[j] + mod$par$b[j] * dat$z))
r_hand <- qnorm(rowSums(stateprobs(mod, forecast = TRUE) *
                        cbind(pgamma2(dat$step, mu_t[, 1], exp(mod$par$log_sigma[1])),
                              pgamma2(dat$step, mu_t[, 2], exp(mod$par$log_sigma[2])))))
cat("1) loop over states + covariate + NA + prior: max diff",
    max(abs(r_auto - r_hand), na.rm = TRUE), "| NA rows", which(is.na(r_auto)), "\n")

## 2) two streams via cbind and product, angle missing at some times
dat2 <- list(step = trex$step[1:500], angle = trex$angle[1:500])
par2 <- list(log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)), log_kappa = log(c(0.5, 2)),
             eta = rep(-2, 2))
nll2 <- function(par) {
  getAll(par, dat2)
  Gamma <- tpm(eta)
  mu <- exp(log_mu); sigma <- exp(log_sigma); kappa <- exp(log_kappa)
  allprobs <- cbind(dgamma2(step, mu[1], sigma[1]), dgamma2(step, mu[2], sigma[2]))
  ind <- which(!is.na(angle))
  for(j in 1:2) allprobs[ind, j] <- allprobs[ind, j] * dvm(angle[ind], 0, kappa[j])
  -forward(stationary(Gamma), Gamma, allprobs)
}
obj2 <- MakeADFun(nll2, par2, silent = TRUE)
opt2 <- nlminb(obj2$par, obj2$fn, obj2$gr)
mod2 <- report(obj2)
r2 <- pseudo_res_auto(mod2)
fc2 <- stateprobs(mod2, forecast = TRUE)
p2 <- mod2$par
h_step <- qnorm(rowSums(fc2 * cbind(pgamma2(dat2$step, exp(p2$log_mu[1]), exp(p2$log_sigma[1])),
                                    pgamma2(dat2$step, exp(p2$log_mu[2]), exp(p2$log_sigma[2])))))
h_angle <- qnorm(rowSums(fc2 * cbind(pvm(dat2$angle, 0, exp(p2$log_kappa[1])),
                                     pvm(dat2$angle, 0, exp(p2$log_kappa[2])))))
cat("2) two streams: names", names(r2), "| max diff step", max(abs(r2$step - h_step), na.rm = TRUE),
    "| angle", max(abs(r2$angle - h_angle), na.rm = TRUE),
    "| NA angle rows match", identical(which(is.na(r2$angle)), which(is.na(dat2$angle))), "\n")

## 3) loop over time, densities on the log scale, logspace = TRUE
nll3 <- function(par) {
  getAll(par, dat)
  Gamma <- tpm(eta)
  sigma <- exp(log_sigma)
  lallprobs <- matrix(0, length(step), 2)
  for(t in which(!is.na(step))) lallprobs[t, ] <- dnorm(step[t], exp(log_mu), sigma, log = TRUE)
  -forward(stationary(Gamma), Gamma, lallprobs, logspace = TRUE)
}
obj3 <- MakeADFun(nll3, par[-2], silent = TRUE)
opt3 <- nlminb(obj3$par, obj3$fn, obj3$gr)
mod3 <- report(obj3)
r3 <- pseudo_res_auto(mod3)
p3 <- mod3$par
h3 <- qnorm(rowSums(stateprobs(mod3, forecast = TRUE) *
                    sapply(1:2, function(j) pnorm(dat$step, exp(p3$log_mu[j]), exp(p3$log_sigma[j])))))
cat("3) loop over time, log scale: max diff", max(abs(r3 - h3), na.rm = TRUE), "\n")

## 4) discrete: Poisson counts, randomised
set.seed(1)
cnt <- rpois(300, rep(c(2, 8), each = 30))
par4 <- list(log_lambda = log(c(1, 5)), eta = rep(-2, 2))
nll4 <- function(par) {
  getAll(par)
  Gamma <- tpm(eta)
  lambda <- exp(log_lambda)
  allprobs <- cbind(dpois(cnt, lambda[1]), dpois(cnt, lambda[2]))
  -forward(stationary(Gamma), Gamma, allprobs)
}
obj4 <- MakeADFun(nll4, par4, silent = TRUE)
opt4 <- nlminb(obj4$par, obj4$fn, obj4$gr)
mod4 <- report(obj4)
set.seed(2); r4 <- pseudo_res_auto(mod4)
set.seed(2); h4 <- suppressMessages(pseudo_res(cnt, "pois", list(lambda = exp(mod4$par$log_lambda)), mod = mod4))
cat("4) discrete Poisson (randomised, same seed): max diff", max(abs(r4 - unclass(h4))), "\n")


## 5) two streams summed on the log scale: dnorm(..., log = TRUE) + dgamma(..., log = TRUE)
set.seed(3)
s5 <- rep(1:2, each = 25, length.out = 400)
d5 <- list(x = rnorm(400, c(0, 3)[s5], 1), y = rgamma(400, shape = c(2, 6)[s5], rate = 1))
par5 <- list(mu = c(0, 3), log_sd = c(0, 0), log_shape = log(c(2, 6)), eta = rep(-2, 2))
nll5 <- function(par) {
  getAll(par, d5)
  Gamma <- tpm(eta)
  logallprobs <- matrix(0, length(x), 2)
  for(j in 1:2) logallprobs[, j] <- dnorm(x, mu[j], exp(log_sd[j]), log = TRUE) +
    dgamma(y, shape = exp(log_shape[j]), rate = 1, log = TRUE)
  -forward(stationary(Gamma), Gamma, logallprobs, logspace = TRUE)
}
obj5 <- MakeADFun(nll5, par5, silent = TRUE)
opt5 <- nlminb(obj5$par, obj5$fn, obj5$gr)
mod5 <- report(obj5)
r5 <- pseudo_res_auto(mod5)
fc5 <- stateprobs(mod5, forecast = TRUE); p5 <- mod5$par
hx <- qnorm(rowSums(fc5 * sapply(1:2, function(j) pnorm(d5$x, p5$mu[j], exp(p5$log_sd[j])))))
hy <- qnorm(rowSums(fc5 * sapply(1:2, function(j) pgamma(d5$y, shape = exp(p5$log_shape[j]), rate = 1))))
cat("5) log-scale sum of two streams: names", names(r5), "| max diff x", max(abs(r5$x - hx)),
    "| y", max(abs(r5$y - hy)), "\n")


## ---------------------------------------------------------------------------
## Bad cases: the check must stop
## ---------------------------------------------------------------------------
try_stop <- function(label, mod) {
  msg <- tryCatch({ pseudo_res_auto(mod); "NO ERROR (bad!)" }, error = function(e) conditionMessage(e))
  cat(label, "->", substr(msg, 1, 90), "\n")
}
fit <- function(nll, par) { o <- MakeADFun(nll, par, silent = TRUE); nlminb(o$par, o$fn, o$gr); report(o) }
pos <- trex$step[1:300]; pos <- pos[!is.na(pos)]

## a) truncation by hand: dnorm / (1 - pnorm(0))
try_stop("a) truncation by hand ", fit(function(par) {
  getAll(par); Gamma <- tpm(eta)
  allprobs <- sapply(1:2, function(j) dnorm(pos, mu[j], 1) / (1 - pnorm(0, mu[j], 1)))
  -forward(stationary(Gamma), Gamma, allprobs)
}, list(mu = c(0.3, 2), eta = rep(-2, 2))))

## b) Jacobian by hand: lognormal via dnorm(log(x)) / x
try_stop("b) Jacobian by hand   ", fit(function(par) {
  getAll(par); Gamma <- tpm(eta)
  allprobs <- cbind(dnorm(log(pos), m[1], 1) / pos, dnorm(log(pos), m[2], 1) / pos)
  -forward(stationary(Gamma), Gamma, allprobs)
}, list(m = c(-1, 0.5), eta = rep(-2, 2))))

## c) zero inflation via an indicator
zcnt <- replace(cnt, sample(300, 40), 0)
try_stop("c) zero inflation     ", fit(function(par) {
  getAll(par); Gamma <- tpm(eta); p0 <- plogis(lp0)
  allprobs <- sapply(1:2, function(j) p0 * (zcnt == 0) + (1 - p0) * dpois(zcnt, exp(log_lambda[j])))
  -forward(stationary(Gamma), Gamma, allprobs)
}, list(log_lambda = log(c(1, 5)), lp0 = -1, eta = rep(-2, 2))))

## d) measurement error integrated out by quadrature inside allprobs
g <- seq(-3, 6, length.out = 40); h <- diff(g)[1]
y <- rnorm(200, rep(c(0, 3), each = 20), 1.2)
try_stop("d) quadrature         ", fit(function(par) {
  getAll(par); Gamma <- tpm(eta)
  allprobs <- sapply(1:2, function(j)
    sapply(y, function(yt) sum(dnorm(yt, g, 0.5) * dnorm(g, mu[j], 1) * h)))
  -forward(stationary(Gamma), Gamma, allprobs)
}, list(mu = c(0, 3), eta = rep(-2, 2))))
