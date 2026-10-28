# Fields of the model object set by lama() or qreml(), not reported by the likelihood function
model_fields <- function(x){
  fields <- c("par", "llk", "df", "nobs", "mgc", "obs", "obj")
  if(inherits(x, "qremlModel")){
    sp <- x$spname
    p <- x$parname
    fields <- c(fields, sp, paste0(c("all_", "map_"), sp), p, paste0("relist_", p), paste0(p, "_vec"),
                "spname", "parname", "n_fixpar", "edf", "edf_valid", "Hessian_conditional",
                "obj_joint", "llk_restricted", "allmods", "converged", "iter", "best_iter",
                "outer_grad", "outer_hessian", "tol", "hessian_repaired", "lambda_at_bound", "max_drop")
  }
  fields
}

# Is Q a generator matrix (or array of them) for N states?
is_generator <- function(Q, N){
  d <- dim(Q)
  if(!is.numeric(Q) || !length(d) %in% 2:3 || any(d[1:2] != N)) return(FALSE)
  rs <- if(length(d) == 2) rowSums(Q) else apply(Q, c(1, 3), sum)
  all(abs(rs) < 1e-8)
}

# State process as reported by forward() and tpm_ct(): number of states, number of tracks,
# continuous time (if a generator Q was reported) and type
state_process <- function(x){
  N <- ncol(x$allprobs)
  if(is.null(N)) return(NULL)
  nTracks <- if(is.null(x$trackID)) 1 else length(unique(x$trackID))
  ct <- is_generator(x$Q, N)
  type <- x$type
  if(ct){ # Gamma holds exp(Q dt) for all time differences, so the type comes from Q
    dQ <- dim(x$Q)
    type <- if(length(dQ) == 2 || dQ[3] == 1) "homogeneous" else
      if(dQ[3] == nTracks) "track-specific" else "inhomogeneous"
  }
  list(N = N, nTracks = nTracks, ct = ct, type = type)
}

# First lines of print() and summary(): state process, observations and parameters
model_header <- function(x){
  sp <- state_process(x)
  line1 <- "Latent Markov model"
  if(!is.null(sp)){
    line1 <- paste0(line1, " with ", sp$N, if(sp$N == 1) " state" else " states",
                    " (", if(sp$ct) "continuous" else "discrete", " time",
                    if(!is.null(sp$type)) paste0(", ", sp$type), ")")
  }

  line2 <- NULL
  if(!is.null(x$nobs)){
    line2 <- paste(x$nobs, "observations")
    nmiss <- nrow(x$allprobs) - x$nobs
    if(nmiss > 0) line2 <- paste0(line2, " (", nmiss, " missing)")
    if(sp$nTracks > 1) line2 <- paste0(line2, " in ", sp$nTracks, " tracks")
  }
  # number of parameters, for qreml() fixed ones and effective degrees of freedom of the smooths
  npar <- if(is.null(x$edf)) paste(x$df, "parameters") else
    paste0(x$n_fixpar, " fixed parameters + smooths with ", round(x$df - x$n_fixpar, 2), " edf")
  nrandom <- length(x$obj$env$random)
  if(nrandom > 0) npar <- paste0(npar, " + ", nrandom, " random effects")

  c(line1, paste(c(line2, npar), collapse = ", "))
}

# Log-likelihood, AIC and BIC (if the number of observations is known) in one line
fit_line <- function(x){
  ll <- logLik(x)
  vals <- c("Log-likelihood" = ll, AIC = AIC(ll), BIC = if(!is.null(x$nobs)) BIC(ll))
  paste0(names(vals), ": ", format(round(vals, 2), nsmall = 2, trim = TRUE), collapse = "   ")
}

#' Print a \code{LaMaModel} object
#'
#' Prints a short overview of a model from \code{\link{lama}} or \code{\link{qreml}}:
#' the state process, observations, parameters, log-likelihood, AIC and BIC, and the names of the reported quantities.
#' Use \code{summary()} for the estimates.
#'
#' @param x model object of class \code{"LaMaModel"}
#' @param ... additional arguments (not used)
#'
#' @returns \code{x}, invisibly
#' @export
print.LaMaModel <- function(x, ...) {
  cat(model_header(x), sep = "\n")

  cat(fit_line(x), sep = "\n")

  # only worth a line if the optimiser has clearly not converged
  if(isTRUE(x$mgc > 0.1)){
    cat("! Maximum gradient component is", signif(x$mgc, 3), "- the optimiser may not have converged\n")
  }

  reported <- setdiff(names(x), c(model_fields(x), "type", "streams", "tpm_formula")) # type is in the first line
  if(length(reported) > 0) cat(strwrap(paste("Reported:", paste(reported, collapse = ", ")), exdent = 2), sep = "\n")

  invisible(x)
}

#' Plot a \code{LaMaModel} object
#'
#' Plots the pseudo-residuals of a model from \code{\link{lama}} or \code{\link{qreml}},
#' computed by \code{\link{pseudo_res}} from the model object alone, with one row of plots per data stream,
#' or the estimated state-dependent densities of each data stream built with \code{\link{statedens}}.
#'
#' @details
#' For \code{type = "statedens"}, each data stream is shown as a histogram with the state-dependent densities
#' (probabilities for discrete distributions) weighted by the state probabilities, and their sum, the marginal density, as a dashed line.
#' The weighted density of state \eqn{j} is the mean over the observations of \eqn{p_{tj} f_j(x; \theta_{tj})},
#' which allows for parameters that depend on covariates, several tracks and inhomogeneous state processes.
#' This needs \code{allprobs} to be the product of matrices from \code{\link{statedens}}, which is checked.
#'
#' Further arguments for \code{type = "statedens"}:
#' \describe{
#'   \item{\code{weights}}{state probabilities \eqn{p_{tj}}: \code{"stateprobs"} (default) for the local state probabilities from \code{\link{stateprobs}},
#'   \code{"viterbi"} for the states decoded by \code{\link{viterbi}}, \code{"stationary"} for the stationary distribution
#'   (averaged over tracks, or over the cycle for periodic models), or a vector of length \code{nStates}}
#'   \item{\code{col}}{colours of the states, defaults to \code{\link{LaMaColors}}}
#'   \item{\code{lwd}}{line width}
#'   \item{\code{n}}{number of grid points for continuous distributions}
#'   \item{\code{legend}}{position of the legend, e.g. \code{"topright"} (default), or \code{FALSE} for none}
#'   \item{\dots}{arguments passed to \code{plot()} of the histogram, e.g. \code{breaks}, \code{main}, \code{xlab}, \code{xlim}, \code{ylim} or \code{border}}
#' }
#' Several streams are shown side by side, unless \code{par(mfrow)} has been set.
#'
#' @param x model object of class \code{"LaMaModel"}
#' @param type \code{"residuals"} (default) for the pseudo-residuals, \code{"statedens"} for the state-dependent densities
#' @param stream optional names of the data streams to plot, by default all
#' @param ... further arguments passed to \code{\link{plot.LaMaResiduals}}, or for \code{type = "statedens"} see \sQuote{Details}
#'
#' @returns for \code{type = "residuals"} the pseudo-residuals, for \code{type = "statedens"} a list with the grid \code{x},
#' the weighted state-dependent densities \code{density} (one column per state), the state \code{weights} and the observations \code{data}
#' (or a list of these for several streams), invisibly
#' @importFrom graphics legend
#' @export
#'
#' @examples
#' dat <- list(step = trex$step[1:500], angle = trex$angle[1:500])
#' nll <- function(par) {
#'   getAll(par, dat)
#'   Gamma <- tpm(eta)
#'   allprobs <- statedens(step, dgamma2, exp(log_mu), exp(log_sigma)) *
#'     statedens(angle, dvm, 0, exp(log_kappa))
#'   -forward(stationary(Gamma), Gamma, allprobs)
#' }
#' par <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)),
#'             log_kappa = log(c(0.2, 1.5)))
#' obj <- MakeADFun(nll, par, silent = TRUE)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' mod <- lama(obj)
#'
#' plot(mod) # pseudo-residuals
#' plot(mod, type = "statedens") # state-dependent densities of both streams
#' plot(mod, type = "statedens", stream = "step", weights = "viterbi", 
#'      breaks = 40, xlab = "step length")
plot.LaMaModel <- function(x, type = c("residuals", "statedens"), stream = NULL, ...) {
  type <- match.arg(type)
  if(type == "statedens") return(invisible(plot_statedens(x, stream, ...)))

  pres <- pseudo_res(x)
  # only the given streams, if there are several
  if(!is.null(stream) && is.list(pres) && all(sapply(pres, inherits, "LaMaResiduals"))) {
    pres <- if(length(stream) == 1) pres[[stream]] else structure(pres[stream], class = "LaMaResiduals")
  }
  plot(pres, ...)
  invisible(pres) # no need to compute (or randomise) them again
}

#' Summary of a \code{LaMaModel} object
#'
#' @description
#' Summarises a model from \code{\link{lama}} or \code{\link{qreml}}.
#' Besides the overview of \code{print()}, this contains
#' \itemize{
#'   \item for a homogeneous state process (or one per track) with at most 10 states: the transition probability matrix
#'   (or generator matrix \code{Q} in continuous time, if reported, e.g. by \code{\link{tpm_ct}}), its stationary distribution and mean dwell times,
#'   averaged over tracks if there is one per track,
#'   \item for a periodic state process: the stationary distribution averaged over the cycle, i.e. the long-run proportion of time spent in each state,
#'   \item the quantities reported by the likelihood function via \code{REPORT()}: numeric vectors of length equal to the number of states
#'   are shown as one table with a column per state, other small objects separately, and large ones only by name,
#'   \item log-likelihood, AIC, BIC and the maximum gradient component at the estimate.
#' }
#'
#' @param object model object of class \code{"LaMaModel"}
#' @param x object of class \code{"summary.LaMaModel"}
#' @param digits number of significant digits to print
#' @param ... additional arguments (not used)
#'
#' @importFrom stats AIC BIC logLik setNames
#' @returns \code{summary()} returns an object of class \code{"summary.LaMaModel"}, a list containing the quantities shown,
#' e.g. \code{Gamma}, \code{stationary}, \code{dwell} and \code{bystate}. \code{print()} returns \code{x} invisibly.
#' @export
summary.LaMaModel <- function(object, ...) {
  sp <- state_process(object)
  s <- list(header = model_header(object), fit = fit_line(object), mgc = object$mgc,
            tpm_formula = object$tpm_formula)
  
  ### state process: matrix, stationary distribution and dwell times, only for few states
  if(!is.null(sp) && sp$N > 10){
    s$note <- paste0("State process not shown for more than 10 states")
  } else if(!is.null(sp) && sp$N > 1){
    N <- sp$N
    states <- paste0("S", 1:N)
    if(isTRUE(sp$type %in% c("homogeneous", "track-specific"))){
      M <- if(sp$ct) object$Q else object$Gamma
      M <- array(M, c(N, N, length(M) / N^2)) # one slice per track
      stat <- if(sp$ct) stationary_ct(M) else stationary(M)
      d <- apply(M, 3, diag)
      dwell <- if(sp$ct) -1 / d else 1 / (1 - d)
      
      # averages over tracks
      M <- rowMeans(M, dims = 2)
      dimnames(M) <- list(states, states)
      if(sp$ct) s$Q <- M else s$Gamma <- M
      s$stationary <- setNames(colMeans(stat), states)
      s$dwell <- setNames(rowMeans(dwell), states)
      s$nTracks <- dim(d)[2]
    } else if(identical(sp$type, "periodic")){
      s$occupancy <- setNames(colMeans(stationary_p(object$Gamma)), states)
    }
  }
  
  ### quantities reported by the likelihood function
  internal <- c("delta", "Gamma", "allprobs", "trackID", "type", "tod", "streams",
                "delta_sparse", "Gamma_sparse", "Delta_sparse", "tpm_formula", if(isTRUE(sp$ct)) "Q")
  reported <- object[setdiff(names(object), c(model_fields(object), internal))]
  
  # vectors with one entry per state form one table
  N <- if(is.null(sp)) 0 else sp$N
  bystate <- vapply(reported, function(v) is.numeric(v) && is.null(dim(v)) && length(v) == N, logical(1))
  if(N > 1 && N <= 10 && any(bystate)){
    s$bystate <- do.call(rbind, reported[bystate])
    colnames(s$bystate) <- paste0("S", 1:N)
    reported <- reported[!bystate]
  }
  
  # other small objects are shown, large ones only by name
  small <- vapply(reported, function(v) {
    d <- dim(v)
    is.atomic(v) && (if(is.null(d)) length(v) <= 10 else length(d) == 2 && all(d <= 10))
  }, logical(1))
  s$other <- reported[small]
  s$hidden <- vapply(names(reported)[!small], function(nm) {
    v <- reported[[nm]]
    paste0(nm, " [", if(is.list(v)) "list" else paste(if(is.null(dim(v))) length(v) else dim(v), collapse = " x "), "]")
  }, character(1), USE.NAMES = FALSE)
  
  class(s) <- "summary.LaMaModel"
  s
}

#' @rdname summary.LaMaModel
#' @export
print.summary.LaMaModel <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat(x$header, sep = "\n")
  if(!is.null(x$tpm_formula)) cat("State process: ", x$tpm_formula, "\n", sep = "")
  
  ### state process
  avg <- if(isTRUE(x$nTracks > 1)) paste0(" (average over ", x$nTracks, " tracks)") else ""
  if(!is.null(x$Gamma)) cat("\nTransition probability matrix", avg, ":\n", sep = "")
  if(!is.null(x$Q)) cat("\nGenerator matrix", avg, ":\n", sep = "")
  if(!is.null(x$Gamma) || !is.null(x$Q)) print(if(is.null(x$Q)) x$Gamma else x$Q, digits = digits)
  if(!is.null(x$stationary)){
    cat("\nStationary distribution", avg, ":\n", sep = "")
    print(x$stationary, digits = digits)
    cat("\nMean dwell times", avg, ":\n", sep = "")
    print(x$dwell, digits = digits)
  }
  if(!is.null(x$occupancy)){
    cat("\nLong-run proportion of time in each state (stationary distribution averaged over the cycle):\n")
    print(x$occupancy, digits = digits)
  }
  if(!is.null(x$note)) cat("\n", x$note, "\n", sep = "")
  
  ### reported quantities
  if(!is.null(x$bystate) || length(x$other) > 0 || length(x$hidden) > 0) cat("\nReported quantities:\n")
  if(!is.null(x$bystate)) print(x$bystate, digits = digits)
  for(nm in names(x$other)){
    v <- x$other[[nm]]
    if(is.null(dim(v)) && is.null(names(v))){
      cat(nm, ": ", paste(format(v, digits = digits), collapse = " "), "\n", sep = "")
    } else{
      cat(nm, ":\n", sep = "")
      print(v, digits = digits)
    }
  }
  if(length(x$hidden) > 0) cat(strwrap(paste("Not shown:", paste(x$hidden, collapse = ", ")), exdent = 2), sep = "\n")
  
  ### fit
  cat("\n", x$fit, "\n", sep = "")
  if(!is.null(x$mgc)) cat("Maximum gradient component: ", signif(x$mgc, 3), "\n", sep = "")
  
  invisible(x)
}

#' Summary of a \code{qremlModel} object
#'
#' @description
#' Summarises a model from \code{\link{qreml}}: the summary of \code{\link{summary.LaMaModel}},
#' plus the effective degrees of freedom and smoothing parameters of the smooth terms and diagnostics of the smoothness selection.
#'
#' @param object model object of class \code{"qremlModel"}
#' @param x object of class \code{"summary.qremlModel"}
#' @param digits number of significant digits to print
#' @param ... additional arguments (not used)
#'
#' @importFrom utils tail
#' @returns \code{summary()} returns an object of class \code{"summary.qremlModel"}, 
#' extending \code{"summary.LaMaModel"} by \code{edf}, \code{lambda} and \code{selection}. \code{print()} returns \code{x} invisibly.
#' @export
summary.qremlModel <- function(object, ...) {
  s <- summary.LaMaModel(object)
  
  # one edf per smoothed coefficient vector, named like the smoothing parameters
  smooths <- names(object$edf)
  if(is.null(smooths)) smooths <- paste0("s", seq_along(object$edf))
  s$edf <- unlist(object$edf, use.names = FALSE)
  names(s$edf) <- unlist(Map(function(e, nm) if(length(e) == 1) nm else paste0(nm, ".", seq_along(e)),
                             object$edf, smooths), use.names = FALSE)
  s$lambda <- object[[object$spname]]
  if(is.null(names(s$lambda))) names(s$lambda) <- paste0(object$spname, ".", seq_along(s$lambda))
  s$spname <- object$spname
  
  ### diagnostics of the smoothness selection, if the fit reports them
  if(!is.null(object$converged)){
    sel <- paste0("Converged: ", object$converged, " (", object$iter, " outer iterations)")
    if(!is.null(object$outer_grad)){
      sel <- c(sel, paste("Maximum absolute outer gradient:", signif(max(abs(object$outer_grad)), 3)))
    }
    # the rest are flags, and only worth a line when they are actually set
    if(isFALSE(object$edf_valid)){
      sel <- c(sel, "! Effective degrees of freedom outside [0, block dimension]; the inverse Hessian is unreliable")
    }
    if(isTRUE(object$hessian_repaired)){
      sel <- c(sel, "! Hessian was not positive definite; a ridge was added")
    }
    if(length(object$lambda_at_bound) > 0){
      sel <- c(sel, paste("!", length(object$lambda_at_bound), "of", length(s$lambda),
                          "smoothing parameters ended at the upper bound"))
    }
    # both of the following are judged against the fit's own 'tol', below which a
    # change in the criterion counts as negligible: after convergence it wobbles,
    # and that is not worth a line. isTRUE() so that a fit without 'tol' is quiet
    small <- object[["tol"]]
    if(isTRUE(object$max_drop > small)){
      sel <- c(sel, paste("! Restricted likelihood decreased by up to", signif(object$max_drop, 4),
                          "during the iteration"))
    }
    llk_r <- object$llk_restricted
    shortfall <- if(length(llk_r) > 0) max(llk_r) - tail(llk_r, 1)
    if(isTRUE(shortfall > small)){
      sel <- c(sel, paste0("! Returned the fit at iteration ", object$best_iter, "; the run ended ",
                           signif(shortfall, 4), " worse than that"))
    }
    s$selection <- sel
  }
  
  class(s) <- c("summary.qremlModel", class(s))
  s
}

#' @rdname summary.qremlModel
#' @export
print.summary.qremlModel <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  NextMethod()
  
  cat("\nSmooth terms:\n")
  if(identical(names(x$edf), names(x$lambda))){ # one smoothing parameter per edf, i.e. no tensor products
    print(cbind(edf = x$edf, "smoothing parameters" = x$lambda), digits = digits)
  } else{
    cat("Effective degrees of freedom:\n")
    print(x$edf, digits = digits)
    cat("Smoothing parameters:\n")
    print(x$lambda, digits = digits)
  }
  
  if(!is.null(x$selection)) cat("\nSmoothness selection:\n", paste0(x$selection, "\n"), sep = "")
  
  invisible(x)
}
