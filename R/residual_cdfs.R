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

# is x a density d* with a CDF p*, seen from env
has_cdf <- function(x, env) {
  grepl("^d", x) && exists(x, envir = env, mode = "function") &&
    exists(sub("^d", "p", x), envir = env, mode = "function")
}

# user-defined functions called in f, recursively, e.g. pnll inside qreml() or helper
# functions (no package functions and no densities), so that their densities can be swapped as well
user_functions <- function(f, found = list()) {
  for(n in unique(all.names(body(f)))) {
    if(n %in% names(found) || has_cdf(n, environment(f)) ||
       !exists(n, envir = environment(f), mode = "function")) next
    g <- get(n, envir = environment(f), mode = "function")
    if(is.primitive(g) || isNamespace(environment(g))) next
    found[[n]] <- g
    found <- user_functions(g, found)
  }
  found
}

# densities d* called in f or the user functions it calls, with a CDF p*
density_names <- function(f) {
  fs <- c(list(f), user_functions(f))
  d <- unique(unlist(lapply(fs, function(g) all.names(body(g)))))
  Filter(function(x) has_cdf(x, environment(f)), d)
}

# allprobs from the likelihood function of obj at par (like obj$par), with each density
# replaced according to mode(stream): "one", "two", "cdf" or "cdf_lower" (F(x - 1))
swapped_allprobs <- function(obj, par, dens, mode) {
  f <- obj_frame(obj)$func
  fenv <- environment(f)
  fams <- list() # families of each stream
  shadow <- new.env(parent = fenv)
  for(d in dens) local({
    dfun <- get(d, envir = fenv, mode = "function")
    pfun <- get(sub("^d", "p", d), envir = fenv, mode = "function")
    assign(d, function(x, ...) {
      # set by statedens(), which calls the density with its own variable
      s <- if(is.null(.lama_env$stream)) stream_name(substitute(x)) else .lama_env$stream
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
  # copies of the user functions called by f, each enclosed by a child of its own environment
  # that holds the swapped densities and the other copies
  ufs <- user_functions(f)
  envs <- lapply(ufs, function(g) new.env(parent = environment(g)))
  copies <- Map(function(g, e) { environment(g) <- e; g }, ufs, envs)
  for(e in envs) for(d in dens) assign(d, get(d, envir = shadow), envir = e)
  for(e in c(list(shadow), envs)) for(n in names(copies)) assign(n, copies[[n]], envir = e)
  environment(f) <- shadow

  .lama_env$record <- list()
  on.exit(.lama_env$record <- NULL)
  parlist <- obj$env$parList(par = par) # par is the full vector, incl. random effects
  attr(parlist, "check.passed") <- NULL
  f(parlist)
  rec <- .lama_env$record
  if(length(rec) == 0) stop("The likelihood function contains no call of forward().")
  if(length(rec) > 1) stop("Only one call of forward() is supported, use trackID for several tracks.")
  list(P = as.matrix(rec[[1]]), fams = fams)
}

# CDF matrices of each data stream that enters allprobs, after checking that swapping is valid:
# list of streams, each with upper (F(x)), lower (F(x - 1) or NULL) and observed (rows)
stream_cdfs <- function(obj, par = obj$env$last.par.best) {
  dens <- density_names(obj_frame(obj)$func)
  if(!length(dens)) stop("No density with a CDF found in the likelihood function.")
  hint <- " Use pseudo_res(obs, dist, par, mod = mod) instead."
  near <- function(P, val) abs(P - val) < 1e-10

  # check 1: all densities 1 -> every cell 1
  one <- swapped_allprobs(obj, par, dens, function(s) "one")
  if(any(!near(one$P, 1))) {
    stop("With all densities set to 1, allprobs is not 1 in ", sum(!near(one$P, 1)), " cells. ",
         "The densities do not enter allprobs as a product or mixture, so their CDFs cannot be used.", hint)
  }

  out <- list()
  for(s in names(one$fams)) {
    # check 2: densities of s set to 2 -> cells 1 (s missing) or 2 (s observed)
    two <- swapped_allprobs(obj, par, dens, function(x) if(x == s) "two" else "one")$P
    if(all(near(two, 1))) next # stream does not enter allprobs, e.g. a prior
    if(any(!near(two, 1) & !near(two, 2))) stop("Stream '", s, "' enters allprobs non-linearly.", hint)
    n2 <- rowSums(near(two, 2))
    if(any(n2 > 0 & n2 < ncol(two))) stop("Stream '", s, "' is observed in some states only.", hint)

    cdf <- function(m) swapped_allprobs(obj, par, dens, function(x) if(x == s) m else "one")$P
    discrete <- any(one$fams[[s]] %in% .discrete_fams)
    out[[s]] <- list(upper = cdf("cdf"), lower = if(discrete) cdf("cdf_lower"),
                     observed = n2 > 0)
  }
  out
}

# pseudo-residuals of each data stream from the model object alone (pseudo_res(mod))
pseudo_res_mod <- function(mod, normal = TRUE, randomise = TRUE) {
  obj <- get_obj(mod)
  par <- est_par(mod)
  # forecast probabilities from the reported quantities
  if(!inherits(mod, c("LaMaModel", "qremlModel"))) mod <- obj$report(par)
  fc <- stateprobs(mod, forecast = TRUE)
  res <- lapply(stream_cdfs(obj, par), function(s) {
    lower <- if(!is.null(s$lower)) rowSums(fc * s$lower)
    finish_res(rowSums(fc * s$upper), lower, s$observed, normal, randomise)
  })
  if(length(res) == 1) return(res[[1]])
  class(res) <- "LaMaResiduals" # list of streams, plotted together
  res
}
