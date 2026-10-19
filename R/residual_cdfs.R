# Residuals from the model object alone, used by pseudo_res(mod):
# the likelihood function is evaluated again with its densities swapped for CDFs,
# so the user's own code puts F_j(x_t) into the cells of allprobs.
# This is valid iff each cell is a combination of the stream's densities with weights
# summing to one (single density, products with other streams, mixtures within a state),
# which is checked by setting all densities to 1 and those of each stream to 2.

# discrete families, for which pseudo-residuals use F(x) and F(x - 1)
.discrete_fams <- c("pois","binom","geom","nbinom","nbinom2", "betabinom","genpois","zibinom",
                    "zibinom2", "zipois","ztbinom","ztnbinom","ztnbinom2", "ztpois")

# name of the data stream of a density call, e.g. step from step[ind] or dat$step
stream_name <- function(e) {
  if(is.call(e) && identical(e[[1]], as.name("$"))) return(as.character(e[[3]]))
  if(is.call(e) && identical(e[[1]], as.name("["))) return(stream_name(e[[2]]))
  all.vars(e)[1]
}

# densities d* called in f with a CDF p*
density_names <- function(f) {
  d <- grep("^d", unique(all.names(body(f))), value = TRUE)
  Filter(function(x) exists(x, envir = environment(f), mode = "function") &&
           exists(sub("^d", "p", x), envir = environment(f), mode = "function"), d)
}

# allprobs from the likelihood function of obj at the estimates, with each density
# replaced according to mode(stream): "one", "two", "cdf" or "cdf_lower" (F(x - 1))
swapped_allprobs <- function(obj, dens, mode) {
  f <- obj_frame(obj)$func
  fenv <- environment(f)
  fams <- list() # families of each stream
  shadow <- new.env(parent = fenv)
  for(d in dens) local({
    dfun <- get(d, envir = fenv, mode = "function")
    pfun <- get(sub("^d", "p", d), envir = fenv, mode = "function")
    assign(d, function(x, ...) {
      s <- stream_name(substitute(x))
      fams[[s]] <<- union(fams[[s]], sub("^d", "", d))
      args <- list(...)
      lg <- isTRUE(args$log)
      args$log <- NULL
      v <- dfun(x, ...) # for the shape of the result
      fill <- function(val) { r <- rep(val, length(v)); dim(r) <- dim(v); r }
      out <- switch(mode(s),
                    one = fill(1),
                    two = fill(2),
                    cdf = do.call(pfun, c(list(x), args)),
                    cdf_lower = do.call(pfun, c(list(x - 1), args)))
      if(lg) log(out) else out
    }, envir = shadow)
  })
  environment(f) <- shadow

  .lama_env$record <- list()
  on.exit(.lama_env$record <- NULL)
  par <- obj$env$parList(obj$env$last.par.best)
  attr(par, "check.passed") <- NULL
  f(par)
  rec <- .lama_env$record
  if(length(rec) == 0) stop("The likelihood function contains no call of forward().")
  if(length(rec) > 1) stop("Only one call of forward() is supported, use trackID for several tracks.")
  list(P = as.matrix(rec[[1]]), fams = fams)
}

# CDF matrices of each data stream that enters allprobs, after checking that swapping is valid:
# list of streams, each with upper (F(x)), lower (F(x - 1) or NULL) and observed (rows)
stream_cdfs <- function(obj) {
  dens <- density_names(obj_frame(obj)$func)
  if(!length(dens)) stop("No density with a CDF found in the likelihood function.")
  hint <- " Use pseudo_res(obs, dist, par, mod = mod) instead."
  near <- function(P, val) abs(P - val) < 1e-10

  # check 1: all densities 1 -> every cell 1
  one <- swapped_allprobs(obj, dens, function(s) "one")
  if(any(!near(one$P, 1))) {
    stop("With all densities set to 1, allprobs is not 1 in ", sum(!near(one$P, 1)), " cells. ",
         "The densities do not enter allprobs as a product or mixture, so their CDFs cannot be used.", hint)
  }

  out <- list()
  for(s in names(one$fams)) {
    # check 2: densities of s set to 2 -> cells 1 (s missing) or 2 (s observed)
    two <- swapped_allprobs(obj, dens, function(x) if(x == s) "two" else "one")$P
    if(all(near(two, 1))) next # stream does not enter allprobs, e.g. a prior
    if(any(!near(two, 1) & !near(two, 2))) stop("Stream '", s, "' enters allprobs non-linearly.", hint)
    n2 <- rowSums(near(two, 2))
    if(any(n2 > 0 & n2 < ncol(two))) stop("Stream '", s, "' is observed in some states only.", hint)

    cdf <- function(m) swapped_allprobs(obj, dens, function(x) if(x == s) m else "one")$P
    discrete <- any(one$fams[[s]] %in% .discrete_fams)
    out[[s]] <- list(upper = cdf("cdf"), lower = if(discrete) cdf("cdf_lower"),
                     observed = n2 > 0)
  }
  out
}

# pseudo-residuals of each data stream from the model object alone (pseudo_res(mod))
pseudo_res_mod <- function(mod, normal = TRUE, randomise = TRUE) {
  obj <- get_obj(mod)
  # forecast probabilities from the reported quantities
  if(!inherits(mod, c("LaMaModel", "qremlModel"))) mod <- obj$report(obj$env$last.par.best)
  fc <- stateprobs(mod, forecast = TRUE)
  res <- lapply(stream_cdfs(obj), function(s) {
    lower <- if(!is.null(s$lower)) rowSums(fc * s$lower)
    finish_res(rowSums(fc * s$upper), lower, s$observed, normal, randomise)
  })
  if(length(res) == 1) res[[1]] else res
}
