#' Sample from a multivariate Gaussian with a sparse precision matrix
#'
#' Draw samples from a multivariate Gaussian distribution specified by a sparse precision matrix.
#' This is numerically efficient for high-dimensional but sparse systems.
#'
#' @param n Number of samples to draw.
#' @param mean Mean vector (or scalar, which will be recycled to match the dimension of \code{Q}).
#' @param Q Sparse precision matrix (\eqn{\Sigma^{-1}}).
#'
#' @returns A matrix of samples with rows corresponding to samples and columns to dimensions.
#' @export
#'
#' @importFrom Matrix Cholesky solve
#' @importFrom stats rnorm
#'
#' @examples
#' rgmrf(3, mean = c(1, 2, 3), Q = Matrix::Diagonal(3))
rgmrf <- function(n, mean = 0, Q) {
  d <- nrow(Q)
  if (length(mean) == 1) mean <- rep(mean, d)
  stopifnot(length(mean) == d)
  
  # efficient sampling (taken from RTMBdist:::rgmrf0)
  L <- Matrix::Cholesky(Q, super = TRUE, LDL = FALSE)
  u <- matrix(rnorm(ncol(L) * n), ncol(L), n)
  u <- Matrix::solve(L, u, system = "Lt")
  u <- Matrix::solve(L, u, system = "Pt")
  
  # adding mean
  samples <- as.matrix(u) + mean
  
  # returning transposed (rows = samples, columsn = dimensions)
  t(samples)
}

# helper functions to MCreport()
transpose_samples <- function(x) {
  nm <- unique(unlist(lapply(x, names)))
  out <- stats::setNames(vector("list", length(nm)), nm)
  for(k in nm) {
    out[[k]] <- lapply(x, function(el) el[[k]])
  }
  out
}
append_report_suffix <- function(report_names, par_names, suffix = ".report") {
  ind <- report_names %in% par_names
  report_names[ind] <- paste0(report_names[ind], suffix)
  report_names
}

#' Sample parameters from approximate Gaussian posterior distribution
#' 
#' Efficient Monte Carlo sampling of parameters (and \code{REPORT}ed quantities) from the approximate posterior of an \code{RTMB} model.
#' See \link[TMB]{sdreport} for details on posterior variance-covariance in random effects models.
#'
#' @param obj Optimised \code{RTMB} object
#' @param nSamples Number of samples to draw
#' @param include_random_pars Logical; Should random parameters be included in the output?
#' @param report Logical; Should reported quantities be samples as well? 
#' Defaults to \code{FALSE} because this may be slow depending on your model.
#' @param ... For internal use only
#'
#' @returns A list of parameter samples, each structured like the initial parameter list from \link[RTMB]{MakeADFun}
#' @export
#' 
#' @importFrom stats optimHess
#' @importFrom Matrix Matrix
#' @importFrom RTMB sdreport
#'
#' @examples
#' step <- trex$step[1:1000] # subsetting trex data
#' N <- 2                    # 2 states
#' 
#' # custom likelihood
#' nll <- function(par) {
#'   getAll(par)
#'   Gamma <- tpm(eta)
#'   delta <- stationary(Gamma)
#'   mu <- exp(log_mu); REPORT(mu)
#'   sigma <- exp(log_sigma); REPORT(sigma)
#'   allprobs <- matrix(1, length(step), N)
#'   for(j in 1:N) allprobs[,j] <- dgamma2(step, mu[j], sigma[j])
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' # initial parameters in named list
#' par0 <- list(eta = rep(-2,2), 
#'              log_mu = log(c(0.3, 1)), 
#'              log_sigma = log(c(0.2, 0.7)))
#'     
#' # constructing AD object        
#' obj <- MakeADFun(nll, par0, silent = TRUE)
#' 
#' # optimising
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' 
#' # sampling from distribution of the MLE
#' samples <- MCreport(obj, nSamples = 10, report = TRUE)
MCreport <- function(obj, 
                     nSamples = 1000,
                     include_random_pars = TRUE,
                     report = FALSE,
                     ...) {
  if(nSamples < 1) stop("At least one sample needs to be drawn")
  
  # Function to map parameter vector to list
  relist_par <- obj$env$parList
  
  # MLE or c(MLE, predicted random effects)
  p_hat <- obj$env$last.par.best
  
  # Are random effects present?
  random_ind <- obj$env$random # index of random effects
  random <- !is.null(random_ind)
  random_names <- unique(names(p_hat[random_ind]))
  
  # Quantities to exclude when REPORT sampling
  excl <- c("allprobs", "type", "trackID")
  # Should allprobs be included?
  dots <- list(...)
  include_allprobs <- isTRUE(dots$include_allprobs)
  if(include_allprobs) excl <- setdiff(excl, "allprobs")
  
  if(random) {
    # model with random effects -> use joint precision approx from TMB
    message("Computing joint precision...")
    Q <- sdreport(obj, getJointPrecision = TRUE,
                  skip.delta.method = TRUE, getReportCovariance = FALSE)$jointPrecision
  } else {
    # model without random effects -> use Hessian
    message("Evaluating Hessian...")
    Hessian <- obj$he(p_hat)
    Q <- Matrix(Hessian, sparse = TRUE)
  }
  
  # Sampling parameter vectors
  samples <- rgmrf(nSamples, p_hat, Q)
  
  # Relisting samples
  par_samples <- lapply(seq_len(nSamples), function(i) {
    p <- samples[i, ]
    relist_par(par = p)
  })
  # nicer shape: named list outside
  par_samples <- transpose_samples(par_samples)
  
  if(!include_random_pars && random) {
    par_samples <- par_samples[!names(par_samples) %in% random_names]
  }
  
  if(report) {
    message("Computing reported quantities...")
    report_samples <- lapply(seq_len(nSamples), function(i) {
      p <- samples[i, ]
      r <- obj$report(par = p)
      r[!names(r) %in% excl] # exclude large, unnecessary quantities
    })
    # nicer shape: named list outside
    report_samples <- transpose_samples(report_samples)
    # adding .report to names when par and report share names to avoid collision when merging
    names(report_samples) <- append_report_suffix(names(report_samples),names(par_samples))
    # merging
    par_samples <- c(par_samples, report_samples)
  }
  
  return(par_samples)
}


