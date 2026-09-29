#' Computes penalty based on quadratic form
#'
#' @description
#' This function computes quadratic penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
#' with smoothing parameters \eqn{\lambda_i}, coefficient vectors \eqn{b_i}, and fixed penalty matrices \eqn{S_i}.
#' 
#' It is intended to be used inside the \strong{penalised negative log-likelihood function} when fitting models with penalised splines or simple random effects via \strong{quasi restricted maximum likelihood} (qREML) with the \code{\link{qreml}} function.
#' For \code{\link{qreml}} to work, the likelihood function needs to be compatible with the \code{RTMB} R package to enable automatic differentiation.
#' 
#' @seealso \code{\link{qreml}} for the \strong{qREML} algorithm
#' 
#' @details
#' \strong{Caution:} The formatting of \code{re_coef} needs to match the structure of the parameter list in your penalised negative log-likelihood function, 
#' i.e. you cannot have two random effect vectors of different names (different list elements in the parameter list), combine them into a matrix inside your likelihood and pass the matrix to \code{penalty}.
#' If these are seperate random effects, each with its own name, they need to be passed as a list to \code{penalty}. Moreover, the ordering of \code{re_coef} needs to match the character vector \code{random} specified in \code{\link{qreml}}.
#' 
#' @references Koslik, J. O. (2024). Efficient smoothness selection for nonparametric Markov-switching models via quasi restricted maximum likelihood. arXiv preprint arXiv:2411.11498.
#'
#' @param re_coef coefficient vector/ matrix or list of coefficient vectors/ matrices
#'
#' Each list entry corresponds to a different smooth/ random effect with its own associated penalty matrix in \code{S}.
#' When several smooths/ random effects of the same kind are present, it is convenient to pass them as a matrix, where each row corresponds to one smooth/ random effect. 
#' This way all rows can use the same penalty matrix.
#' @param S fixed penalty matrix or list of penalty matrices matching the structure of \code{re_coef} and also the dimension of the individuals smooths/ random effects
#' @param lambda penalty strength parameter vector that has a length corresponding to the \strong{total number} of random effects/ spline coefficients in \code{re_coef}
#'
#' E.g. if \code{re_coef} contains one vector and one matrix with 4 rows, then \code{lambda} needs to be of length 5.
#'
#' @return returns the penalty value and reports to \code{\link{qreml}}.
#' @export
#' 
#' @import RTMB
#'
#' @examples
#' # Example with a single random effect
#' re = rep(0, 5)
#' S = diag(5)
#' lambda = 1
#' penalty(re, S, lambda)
#'
#' # Example with two random effects, 
#' # where one element contains two random effects of similar structure
#' re = list(matrix(0, 2, 5), rep(0, 4))
#' S = list(diag(5), diag(4))
#' lambda = c(1,1,2) # length = total number of random effects
#' penalty(re, S, lambda)
#' 
#' # Full model-fitting example
#' \donttest{
#' data = trex[1:1000,] # subset
#'
#' # initial parameter list
#' par = list(logmu = log(c(0.3, 2.5)), # step mean
#'            logsigma = log(c(0.3, 1.5)), # step sd
#'            beta0 = c(-2,-2), # state process intercept
#'            betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
#'           
#' # data object with initial penalty strength lambda
#' dat = list(step = data$step, # step length
#'            tod = data$tod, # time of day covariate
#'            N = 2, # number of states
#'            lambda = rep(10,2)) # initial penalty strength
#'
#' # building model matrices
#' modmat = make_matrices(~ s(tod, bs = "cp"), 
#'                        data = data.frame(tod = 1:24), 
#'                        knots = list(tod = c(0,24))) # wrapping points
#' dat$Z = modmat$Z # spline design matrix
#' dat$S = modmat$S # penalty matrix
#'
#' # penalised negative log-likelihood function
#' pnll = function(par) {
#'   getAll(par, dat) # makes everything contained available without $
#'   Gamma = tpm_g(Z, cbind(beta0, betaspline)) # transition probabilities
#'   delta = stationary_p(Gamma, t = 1) # initial distribution
#'   mu = exp(logmu) # step mean
#'   sigma = exp(logsigma) # step sd
#'   # calculating all state-dependent densities
#'   allprobs = matrix(1, nrow = length(step), ncol = N)
#'   ind = which(!is.na(step)) # only for non-NA obs.
#'   for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
#'   -forward_g(delta, Gamma[,,tod], allprobs) +
#'       penalty(betaspline, S, lambda) # this does all the penalization work
#' }
#'
#' # model fitting
#' mod = qreml(pnll, par, dat, random = "betaspline")
#' }
penalty = function(re_coef, S, lambda) {
  # Convert re_coef to a list of matrices (even if originally a vector)
  if (!is.list(re_coef)) {
    re_coef = list(re_coef)
  }
  
  re_coef = lapply(re_coef, function(x) {
    if (is.null(dim(x))) {
      matrix(x, nrow = 1)  # Convert vectors to 1-row matrices
    } else {
      x  # Leave matrices unchanged
    }
  })
  
  # Get number of distinct random effects (of the same structure)
  n_re = length(re_coef)
  
  # Ensure S is a list of length n_re, replicating it if necessary
  if (!is.list(S)) {
    S = list(S)
  }
  if (length(S) == 1) {
    S = rep(S, n_re)
  }
  
  # transpose if necessary to match S
  re_coef <- lapply(seq_len(n_re), function(i) {
    if (ncol(re_coef[[i]]) != nrow(S[[i]])) {
      t(re_coef[[i]])
    } else if (nrow(re_coef[[i]]) != nrow(S[[i]])) {
      re_coef[[i]]
    } else{
      stop("The coefficient structure does not match the provided penalty matrices.")
    }
  })
  
  # Get the number of similar random effects for each distinct random effect
  re_lengths = sapply(re_coef, nrow)  # All elements are matrices now
  
  # check if lambda has the correct length
  n_lambdas <- sum(re_lengths)
  if(length(lambda) != n_lambdas){
    stop(paste0("The length of lambda (", length(lambda), ") does not match the total number of random effects (", n_lambdas, ")."))
  }
  
  # Precompute start and end indices for lambda
  end = cumsum(re_lengths)
  start = c(1, end[-length(end)] + 1)
  
  RTMB::REPORT(S) # Report penalty matrix list
  
  # Initialize penalty variables
  Pen = vector("list", n_re)
  pen = 0
  
  # Loop over distinct random effects - each now a matrix
  for (i in 1:n_re) {
    current_re = re_coef[[i]]  # current_re is always a matrix now
    
    # Vectorized calculation of penalty for each random effect
    quadform = rowSums(current_re %*% S[[i]] * current_re)
    Pen[[i]] = quadform
    
    # Apply lambda directly using precomputed indices
    pen = pen + sum(lambda[start[i]:end[i]] * quadform)
  }
  
  RTMB::REPORT(Pen) # Report the penalty list for qreml update
  
  pen = 0.5 * pen
  RTMB::REPORT(pen)
  pen
}

#' Computes generalised quadratic-form penalties
#'
#' @description
#' This function computes a quadratic penalty of the form
#' \deqn{0.5 \sum_{i} \lambda_i b^T S_i b,}
#' with smoothing parameters \eqn{\lambda_i}, coefficient vector \eqn{b}, and fixed penalty matrices \eqn{S_i}.
#' This generalises the \code{\link{penalty}} by allowing subsets of the coefficient vector  \eqn{b} to be penalised multiple times with different smoothing parameters, which is necessary for \strong{tensor products}, \strong{functional random effects} or \strong{adaptive smoothing}.
#' 
#' It is intended to be used inside the \strong{penalised negative log-likelihood function} when fitting models with penalised splines or simple random effects via \strong{quasi restricted maximum likelihood} (qREML) with the \code{\link{qreml}} function.
#' For \code{\link{qreml}} to work, the likelihood function needs to be compatible with the \code{RTMB} R package to enable automatic differentiation.
#' 
#' @seealso \code{\link{qreml}} for the \strong{qREML} algorithm
#' 
#' @details
#' \strong{Caution:} The formatting of \code{re_coef} needs to match the structure of the parameter list in your penalised negative log-likelihood function, 
#' i.e. you cannot have two random effect vectors of different names (different list elements in the parameter list), combine them into a matrix inside your likelihood and pass the matrix to \code{penalty}.
#' If these are seperate random effects, each with its own name, they need to be passed as a list to \code{penalty}. Moreover, the ordering of \code{re_coef} needs to match the character vector \code{random} specified in \code{\link{qreml}}.
#' 
#'
#' @param re_coef list of coefficient vectors/ matrices
#'
#' Each list entry corresponds to a different smooth/ random effect with its own associated penalty matrix or penalty-matrix list in \code{S}.
#' When several smooths/ random effects of the same kind are present, it is convenient to pass them as a matrix, where each row corresponds to one smooth/ random effect. 
#' This way all rows can use the same penalty matrix.
#' @param S list of fixed penalty matrices matching the structure of \code{re_coef}. 
#' 
#' This means if \code{re_coef} is of length 3, then \code{S} needs to be a list of length 3. Each entry needs to be either a penalty matrix, matching the dimension of the corresponding entry in \code{re_coef}, or a list with multiple penalty matrices for tensor products.
#' @param lambda penalty strength parameter vector that has a length corresponding to the provided \code{re_coef} and \code{S}. 
#' 
#' Specifically, for entries with one penalty matrix, \code{nrow(re_coef[[i]])} parameters are needed. For entries with \code{k} penalty matrices, \code{k * nrow(re_coef[[i]])} parameters are needed.
#'
#' E.g. if \code{re_coef[[1]]} is a vector and \code{re_coef[[2]]} a matrix with 4 rows, 
#' \code{S[[1]]} is a list of length 2 and \code{S[[2]]} is a matrix, then \code{lambda} needs to be of length 1 * 2 + 4 = 6.
#'
#' @return returns the penalty value and reports to \code{\link{qreml}}.
#' @export
#' 
#' @import RTMB
#'
#' @examples
#' # Example with a single random effect
#' re = rep(0, 5)
#' S = diag(5)
#' lambda = 1
#' penalty(re, S, lambda)
#'
#' # Example with two random effects, 
#' # where one element contains two random effects of similar structure
#' re = list(matrix(0, 2, 5), rep(0, 4))
#' S = list(diag(5), diag(4))
#' lambda = c(1,1,2) # length = total number of random effects
#' penalty(re, S, lambda)
#' 
#' # Full model-fitting example
#' \donttest{
#' data = trex[1:1000,] # subset
#'
#' # initial parameter list
#' par = list(logmu = log(c(0.3, 2.5)), # step mean
#'            logsigma = log(c(0.3, 1.5)), # step sd
#'            beta0 = c(-2,-2), # state process intercept
#'            betaspline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
#'           
#' # data object with initial penalty strength lambda
#' dat = list(step = data$step, # step length
#'            tod = data$tod, # time of day covariate
#'            N = 2, # number of states
#'            lambda = rep(10,2)) # initial penalty strength
#'
#' # building model matrices
#' modmat = make_matrices(~ s(tod, bs = "cp"), 
#'                        data = data.frame(tod = 1:24), 
#'                        knots = list(tod = c(0,24))) # wrapping points
#' dat$Z = modmat$Z # spline design matrix
#' dat$S = modmat$S # penalty matrix
#'
#' # penalised negative log-likelihood function
#' pnll = function(par) {
#'   getAll(par, dat) # makes everything contained available without $
#'   Gamma = tpm_g(Z, cbind(beta0, betaspline)) # transition probabilities
#'   delta = stationary_p(Gamma, t = 1) # initial distribution
#'   mu = exp(logmu) # step mean
#'   sigma = exp(logsigma) # step sd
#'   # calculating all state-dependent densities
#'   allprobs = matrix(1, nrow = length(step), ncol = N)
#'   ind = which(!is.na(step)) # only for non-NA obs.
#'   for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
#'   -forward_g(delta, Gamma[,,tod], allprobs) +
#'       penalty(betaspline, S, lambda) # this does all the penalisation work
#' }
#'
#' # model fitting
#' mod = qreml(pnll, par, dat, random = "betaspline")
#' }
penalty2 = function(re_coef, # coefficient vector/ matrix or list of coefficient vectors/ matrices
                    S, # always needs to be a list: matrix entries for smooths with one penalty matrix, list (length 2) entries for 2D tensorproducts
                    lambda)
{
  # RTMB stuff to avoid annyoing problems
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  
  ## If re_coef is not a list -> list it
  if (!is.list(re_coef)) {
    re_coef = list(re_coef)
  }
  
  ## Convert re_coef to a list of matrices (even if originally a vector)
  re_coef = lapply(re_coef, function(x) {
    if (is.null(dim(x))) {
      matrix(x, nrow = 1)  # Convert vectors to 1-row matrices
    } else {
      x  # Leave matrices unchanged
    }
  })
  
  ## Get number of distinct random effects (of the same structure)
  # i.e. number of random effects with own penalty matrix (list)
  n_re = length(re_coef)
  
  ## Get the number of similar random effects for each distinct random effect
  re_lengths = sapply(re_coef, nrow)  # All elements are matrices now
  
  ## find how many penalty strength pars are needed for each random effect
  # 1: univariate smooth
  # >1: tensorproduct
  n_penalties = sapply(S, function(x){
    if(is.matrix(x)){
      return(1)
    } else{
      return(length(x))
    }
  })
  
  ## Compte indices of simple univariate smooths and of tensorproduct smooths
  simple_ind = which(n_penalties == 1)
  tp_ind = which(n_penalties > 1)
  
  ## total number of lambdas for each random effect with one penalty matrix/list
  lambda_lengths = n_penalties * re_lengths
  
  RTMB::REPORT(S) # Report penalty matrix list (potentially nested)
  
  ## reshape lambdas to list of vectors
  Lambda = reshape_lambda(lambda_lengths, lambda)
  
  ## Initialise penalty variables
  Pen = vector("list", length(n_re))
  # this will get filled either by vector, or by matrix of evaluated penalties b^t S b
  pen = 0
  
  ## Loop over distinct random effects - each now a matrix
  # first simple random effects
  for(ind in seq_along(simple_ind)){
    i = simple_ind[ind] # get original index
    thislambda = Lambda[[i]] # extract lambdas for this smooth
    # extract coefficients for this smooth
    current_re = re_coef[[i]]  # current_re is always a matrix now
    
    # Vectorised calculation of penalty for each random effect
    quadform = rowSums(current_re %*% S[[i]] * current_re)
    Pen[[i]] = quadform
    
    # Apply lambda directly using precomputed indices
    pen = pen + sum(thislambda * quadform)
  }
  # then tensorproducts
  for(ind in seq_along(tp_ind)){
    i = tp_ind[ind] # extract original index
    # extract coefficients for this smooth: these will be penalised by multiple matrices
    this_tp = re_coef[[i]]  # current_re is always a matrix now
    
    # extract penalty matrix
    thisS = S[[i]] # currently penalty matrix list
    n_pen = length(thisS) # number of lambdas for this smooth
    
    # initialise penalty vector for this re (i)
    thispen = vector("list", re_lengths[i])
    counter = 0
    for(j in seq_len(re_lengths[i])) {
      # extract sub-vector of n_pen lambdas
      thislambda = Lambda[[i]][counter + 1:n_pen]
      counter = counter + n_pen
      
      # calculate all penalties separately
      subthispen = numeric(n_pen)
      for(k in seq_len(n_pen)) {
        # compute quadratic form b^t S_k b (WITHOUT lambda_k, same as for simple smooths above)
        subthispen[k] = t(this_tp[j,]) %*% thisS[[k]] %*% this_tp[j,]
      }
      # add to overall penalty (this is where lambda enters)
      pen = pen + sum(thislambda * subthispen)
      
      thispen[[j]] <- subthispen # in this case we have multiple penalties we need for the update
    }
    
    Pen[[i]] <- thispen
  }
  
  RTMB::REPORT(Pen) # Report the penalty list for qreml update
  pen = 0.5 * pen
  RTMB::REPORT(pen)
  pen
}


#' Extract log-likelihood from qremlModel object
#' @param object A fitted model of class "qremlModel"
#' @param ... Additional arguments (not used)
#' @return An object of class "logLik"
#' @export
logLik.qremlModel <- function(object, ...) {
  ll <- object$llk  # your stored log-likelihood
  df <- object$df # number of free parameters
  nobs <- object$nobs  # number of observations
  
  val <- as.numeric(ll)
  attr(val, "df") <- df
  attr(val, "nobs") <- nobs
  class(val) <- "logLik"
  val
}


#' Summary method for \code{qremlModel} objects
#'
#' @description
#' Prints a summary of a model object created by \code{\link{qreml}}.
#'
#' @param object \code{qremlModel} object created by \code{\link{qreml}}
#' @param ... additional arguments
#'
#' @returns prints a summary of the model object
#' 
#' @importFrom stats AIC
#' @importFrom stats BIC
#' 
#' @export
#'
#' @examples
#' # no examples
summary.qremlModel <- function(object, ...) {

  ### Printing state process parameters
  if(!is.null(object$Gamma) | !is.null(object$delta)){
    
    if(ncol(object$Gamma) <= 15){ # don't print for very large state-spaces
      cat("State process parameters:\n")
      if (!is.null(object$Gamma)) {
        if(is.matrix(object$Gamma)){
          cat("\nTransition probability matrix:\n")
          print(object$Gamma)
        } else if(length(dim(object$Gamma)) == 3){
          cat("\nFirst transition probability matrix (t = 1):\n")
          print(object$Gamma[,,1])
        }
      }
      if (!is.null(object$delta)) {
        if(is.vector(object$delta)){
          cat("\nInitial state distribution:\n")
          print(object$delta)
        } else if(is.matrix(object$delta)){
          cat("\nFirst initial state distribution:\n")
          print(object$delta[1,])
        }
      }
      
      cat("\n---")
    }
  }
  
  
  ### Printing effective degrees of freedom if present
  if (!is.null(object$edf)) {
    smoothnames <- names(object$edf)
    if(is.null(smoothnames)){
      smoothnames <- paste0("s.", 1:length(object$edf))
    }
    
    cat("\nEffective degrees of freedom:\n")
    
    cat("\nFixed effects:", object$n_fixpar)
    for(i in seq_along(object$edf)){
      cat("\n", smoothnames[i], ": ", sep = "")
      cat(object$edf[[i]])
    }
    cat("\nTotal: ", object$df, "\n")
    
    cat("\n---")
  }
  
  ### Printing log-likelihood
  cat("\nLog-Likelihood:", object$llk, "\n")
  
  ### Printing Print AIC and BIC
  # AIC
  suppressMessages(aic <- AIC(object))
  # BIC
  bic <- tryCatch(
    suppressMessages(BIC(object)),
    error = function(e) "could not be determined"
  )
  
  cat("AIC:", aic, "  ")
  cat("BIC:", bic, "\n")
  cat("\n---")
  
  ### Printing smoothing parameter estimates
  cat("\nSmoothing parameters:")
  lambdas <- object[[object$spname]]
  lambda_names <- names(lambdas)
  if(is.null(lambda_names)){
    lambda_names <- paste0("lambda.", 1:length(lambdas))
  }
  for(i in seq_along(lambdas)){
    cat("\n", lambda_names[i], ": ", sep = "")
    cat(lambdas[[i]])
  }
  cat("\n")
  
  ### Printing diagnostics of the smoothness selection, if the fit reports them
  if(!is.null(object$converged)){
    cat("\n---")
    cat("\nSmoothness selection:\n")
    cat("Converged:", object$converged, paste0("(", object$iter, " outer iterations)\n"))
    if(!is.null(object$outer_grad)){
      cat("Max. absolute outer gradient:", signif(max(abs(object$outer_grad)), 4), "\n")
    }
    
    # the rest are flags, and only worth a line when they are actually set
    if(isFALSE(object$edf_valid)){
      cat("! Effective degrees of freedom outside [0, block dimension];",
          "the inverse Hessian is unreliable\n")
    }
    if(isTRUE(object$hessian_repaired)){
      cat("! Hessian was not positive definite; a ridge was added\n")
    }
    if(length(object$lambda_at_bound) > 0){
      cat("!", length(object$lambda_at_bound), "of", length(lambdas),
          "smoothing parameters ended at the upper bound\n")
    }
    # both of the following are judged against the fit's own 'tol', below which a
    # change in the criterion counts as negligible: after convergence it wobbles,
    # and that is not worth a line. isTRUE() so that a fit without 'tol' is quiet
    small <- object[["tol"]]
    if(isTRUE(object$max_drop > small)){
      cat("! Restricted likelihood decreased by up to", signif(object$max_drop, 4),
          "during the iteration\n")
    }
    shortfall <- max(object$llk_restricted) - tail(object$llk_restricted, 1)
    if(isTRUE(shortfall > small)){
      cat(paste0("! Returned the fit at iteration ", object$best_iter,
                 "; the run ended ", signif(shortfall, 4), " worse than that\n"))
    }
  }
  
  # Print additional user-specified objects, excluding unwanted ones
  excluded <- c("allprobs", "trackID", "type", "obj", "outer_gr", 
                paste0("all_", object$spname), "parname", object$parname, paste0("relist_", object$parname), 
                paste0("map_", object$spname), "spname", paste0(object$parname, "_vec"), 
                "edf", "Hessian_conditional", "obj_joint",
                "beta", "delta", "Gamma", "lambda", "llk", "n_fixpar", "df", "nobs",
                "llk_restricted", "allmods",
                # reported in the smoothness selection block above
                "converged", "iter", "best_iter", "outer_grad", "outer_hessian", "tol",
                "edf_valid",
                "hessian_repaired", "lambda_at_bound", "max_drop")
  
  remaining_names <- setdiff(names(object), excluded)
  
  if(length(remaining_names) > 0){
    cat("\n---")
    cat("\nOther reported quantities:\n")
  }
  count = 1
  for (name in remaining_names) {
    if (!is.null(object[[name]]) & count <= 10) {
      this <- object[[name]]
      
      # check if object is a matrix, if so, not print if too large
      if(is.matrix(this)){
        if(nrow(this) <= 10 & ncol(this) <= 10){
          cat(name, ":\n", sep = "")
          print(round(this, 4))
        } else{
          cat(name, ": [large matrix, not displayed]\n")
        }
        
      # check if object is a vector, if so, not print if too long
      } else if(is.vector(this)){
        if(length(this) <= 20){
          cat(name, ":\n", sep = "")
          print(round(this, 4))
        } else{
          cat(name, ": [large vector, not displayed]\n")
        }
        
      }
    }
    count = count + 1
  }
  
  # invisible(object)
}


#' Report uncertainty of the estimated smoothing parameters or variances
#' 
#' Computes standard deviations for the smoothing parameters of a model object returned by \code{qreml} using the delta method.
#' 
#' @details
#' The computations are based on the approximate gradient of the restricted log likelihood. The outer Hessian is computed by finite differencing of this gradient. If the inverse smoothing parameters are requested, the standard deviations are transformed to the variances using the delta method.
#' 
#'
#' @param mod model objects as returned by \code{\link{qreml}}
#' @param invert optional logical; if \code{TRUE}, the inverse smoothing paramaters (variances) are returned along with the transformed standard deviations obtained via the delta method.
#'
#' @return list containing \code{report} matrix summarising parameters and standard deviations as well as the outer \code{Hessian} matrix.
#' @export
#' 
#' @importFrom numDeriv jacobian
#' @importFrom MASS ginv
#'
#' @examples
#' ## no examples
sdreport_outer <- function(mod, invert = FALSE){
  if(!inherits(mod, "qremlModel")){
    stop("Model object is not of class 'qremlModel'")
  }
  
  spname <- mod$spname
  map_lambda <- mod[[paste0("map_", spname)]]

  # map lambda
  lambda_mapped <- map_lambda(mod[[spname]])
  
  # map names of lambda
  lambda_names <- names(mod[[spname]])
  names(lambda_names) <- lambda_names
  mapped_names <- map_lambda(lambda_names)

  # qreml() reports the outer Hessian itself; qreml_old() only the gradient.
  # Indexed with [[ ]], since $ would partially match 'outer_grad'
  H <- if(!is.null(mod[["outer_hessian"]])){
    mod[["outer_hessian"]]()
  } else if(!is.null(mod[["outer_gr"]])){
    -jacobian(mod[["outer_gr"]], lambda_mapped, method = "simple")
  } else stop("Model object carries neither 'outer_hessian' nor 'outer_gr'")

  I <- ginv(H)
  vars <- diag(I)
  
  self <- list()
  
  if(invert){
    vars <- vars * lambda_mapped^(-4)
    sds <- sqrt(vars)
    self$report <- rbind(par = 1/lambda_mapped, sd = sds)
  } else{
    sds <- sqrt(vars)
    self$report <- rbind(par = lambda_mapped, sd = sds)
  }
  colnames(self$report) <- mapped_names
  
  self$Hessian <- H
  class(self) = "sdreport_outer"
  self
}


# detect_cycling <- function(param_matrix, threshold = 1, window_size = 10) {
#   n_params <- ncol(param_matrix)
#   cycling_flags <- rep(FALSE, n_params)
#   
#   for (j in 1:n_params) {
#     rolling_var <- sapply(window_size:nrow(param_matrix), function(i) var(param_matrix[(i-window_size+5):i, j]))
#     if (mean(rolling_var) > threshold) {
#       cycling_flags[j] <- TRUE
#     }
#   }
#   return(any(cycling_flags))
# }


#' Penalty approximation of unimodality constraints for univariates smooths
#'
#' @param coef coefficient vector of matrix on which to apply the unimodality penalty
#' @param m vector of indices for the position of the coefficient mode. 
#' If \code{coef} is a vector, must be of length 1. Otherwise, must be of length equal to nrow(coef)
#' @param kappa global scaling factor for the penalty
#' @param concave logical; if \code{TRUE} (default), the penalty enforces increasing until the mode then decreasing. If the coefficients should decrease until the mode, then increase, set \code{concave = FALSE}.
#' @param rho control parameter for smooth approximation to \code{min(x, 0)} used internally. 
#' For large values, gets closer to true minimum function but less stable. 
#'
#' @returns a numeric value of the penalty for the given coefficients
#' @export
#'
#' @examples
#' ## coefficient vector
#' coef <- c(1, 2, 3, 2, 1)
#' # mode at position 3
#' penalty_uni(coef, m = 3) # basically zero
#' #' # mode at position 2
#' penalty_uni(coef, m = 2) # large positive penalty
#' 
#' ## coefficient matrix
#' coef <- rbind(coef, coef)
#' m <- c(1, 4)
#' penalty_uni(coef, m)
penalty_uni <- function(coef, 
                        m, 
                        kappa = 1e3, 
                        concave = TRUE,
                        rho = 20) {
  
  N <- length(m) # number of states
  if(is.null(dim(coef))){
    coef <- matrix(coef, nrow = 1, ncol = length(coef))
  }
  
  if(nrow(coef) != N) {
    stop("Coefficient matrix must have as many rows as there are states.")
  }
  k <- ncol(coef) + 1 # number of coefficients
  
  if(!concave) coef <- -coef # if concave == FALSE, flip coefficients to get convexity penalty
  
  # set up constraint matrices
  C <- construct_C(m, k, exclude_last = TRUE)
  
  # compute penalty by summing over states
  pen <- 0
  for(i in 1:N) {
    pen <- pen - sum(min0_smooth(C[[i]] %*% coef[i,], rho = rho))
  }
  
  # return result scaled by kappa
  kappa * pen
}
