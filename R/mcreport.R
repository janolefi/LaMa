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

# Compute joint precision matrix from RTMB model
getJointPrecision <- function(obj) {
  q_hat <- obj$env$last.par.best  # joint mode
  par_names <- names(q_hat)       # parmeter names
  n <- length(q_hat)
  r <- obj$env$random # index of random effects
  nonr <- setdiff(seq_along(q_hat), r) # index of fixed effects
  theta_hat <- q_hat[nonr] # mode of fixed effects
  
  # Hessian block for fixed effects using finite differences
  message("Computing marginal Hessian...")
  H_Bhat <- stats::optimHess(theta_hat, obj$fn, obj$gr) # Hessian of marginal posterior
  
  # Hessian of random effects at joint mode using AD.
  H_AA <- obj$env$spHess(q_hat, random = TRUE)
  
  # Second derivatives of the joint posterior at the joint
  # mode for the fixed:random effect elements only. Uses AD.
  message("\nEvaluating cross-derivatives...")
  H_AB <- obj$env$f(q_hat, order = 1, type = "ADGrad", keepx=nonr, keepy=r) ## TMBad only !!!
  H_BA <- t(H_AB)
  
  # Numerically efficient way to compute H_BA %*% solve(t(H_AA)) %*% H_AB + H_Bhat
  # way faster than solve(t(H_AA))
  message("\nSolving system...")
  X <- Matrix::solve(H_AA, H_AB)
  H_BB <- H_BA %*% X + H_Bhat
  
  # Building joint precision
  Q <- rbind(
    cbind(H_AA, H_AB),
    cbind(H_BA, H_BB)
  )
  
  rownames(Q) <- colnames(Q) <- par_names
  
  gc()
  
  return(Q)
}


#' Sample parameters from approximate Gaussian posterior distribution
#' 
#' Efficient Monte Carlo sampling of parameters from the approximate posterior of an \code{RTMB} model.
#' See \link[TMB]{sdreport} for details on posterior variance-covariance in random effects models.
#'
#' @param obj Optimised \code{RTMB} object
#' @param nSamples Number of samples to draw
#' @param sample_random_effects Logical; should random effects be sampled? Ignored if the model has no random effects.
#'
#' @returns A list of parameter samples, each structured like the initial parameter list from \link[RTMB]{MakeADFun}
#' @export
#' 
#' @importFrom stats optimHess
#' @importFrom Matrix Matrix
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
#'   mu <- exp(log_mu)
#'   sigma <- exp(log_sigma)
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
#' par_samples <- MCreport(obj, nSamples = 10)
#' # each entry has same structure as par0
#' 
#' # e.g. extracting mean samples
#' mus <- lapply(par_samples, function(p) exp(p$log_mu))
MCreport <- function(obj, 
                     nSamples = 1000, 
                     sample_random_effects = TRUE) {
  # Function to map parameter vector to list
  relist_par <- obj$env$parList
  
  # Are random effects present?
  random_ind <- obj$env$random # index of random effects
  random <- !is.null(random_ind)
  
  # MLE or c(MLE, predicted random effects)
  p_hat <- obj$env$last.par.best
  
  if(random) {
    if(sample_random_effects) {
      Q <- getJointPrecision(obj) # build joint precision matrix
    } else{
      # only sample using marginal Hessian
      random_names <- unique(names(p_hat[random_ind]))
      fixed_ind <- setdiff(seq_along(p_hat), random_ind) # index of fixed effects
      p_hat <- p_hat[fixed_ind]
      message("Computing marginal Hessian...")
      Hessian <- optimHess(p_hat, obj$fn, obj$gr) # Hessian of marginal posterior
      Q <- Matrix(Hessian, sparse = TRUE)
      par_samples <- lapply(seq_len(nSamples),
                            function(i) {
                              out <- relist_par(x = samples[i, ])
                              out[!names(out) %in% random_names]
                            })
      return(par_samples)
    }
    
  } else {
    Hessian <- obj$he(obj$env$last.par.best)
    Q <- Matrix(Hessian, sparse = TRUE)
  }
  
  # Sampling parameter vectors
  samples <- rgmrf(nSamples, p_hat, Q)
  
  # Relisting samples
  par_samples <- lapply(seq_len(nSamples),
                        function(i) {
                          relist_par(par = samples[i, ])
                        })
  
  gc() # cleaning up
  
  return(par_samples)
}
