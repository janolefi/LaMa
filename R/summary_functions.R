# Fields of the model object set by report() or qreml(), not reported by the likelihood function
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
#' Prints a short overview of a model from \code{\link{report}} or \code{\link{qreml}}:
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

  reported <- setdiff(names(x), c(model_fields(x), "type")) # type is in the first line
  if(length(reported) > 0) cat(strwrap(paste("Reported:", paste(reported, collapse = ", ")), exdent = 2), sep = "\n")

  invisible(x)
}
