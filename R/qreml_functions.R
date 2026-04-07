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
#' par = list(logmu = log(c(0.3, 1)), # step mean
#'            logsigma = log(c(0.2, 0.7)), # step sd
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
  # Capture the argument name used in the call to `penalty`
  # current_name <- as.character(substitute(lambda))
  # 
  # # Try to recover the original name by searching the parent frame
  # recover_original_name <- function(value, env) {
  #   # Search for all objects in the parent frame
  #   all_objects <- ls(env)
  #   for (obj in all_objects) {
  #     # Check if the object matches the value of `lambda`
  #     if (identical(get(obj, envir = env), value)) {
  #       return(obj)  # Return the name of the matching object
  #     }
  #   }
  #   NULL  # Return NULL if no match is found
  # }
  # 
  # # Attempt to recover the original name
  # original_name <- recover_original_name(lambda, parent.frame())
  # if (!is.null(original_name)) {
  #   argname_lambda <- original_name
  # } else {
  #   argname_lambda <- current_name
  # }
  # 
  # # Store the recovered name in the penalty metadata environment
  # assign("argname_lambda", argname_lambda, envir = penalty_metadata)
  
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
  
  # check if re_coef and S match
  # if(any(sapply(re_coef, ncol) != sapply(S, nrow))){
  #   stop("The coefficient structure does not match the provided penalty matrices.\n Are the coefficients arranged by row?")
  # }
  
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

#' Quasi restricted maximum likelihood (qREML) algorithm for models with penalised splines or simple i.i.d. random effects
#'
#' @description
#' This algorithm can be used very flexibly to fit statistical models that involve \strong{penalised splines} or simple \strong{i.i.d. random effects}, i.e. that have penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
#' with smoothing parameters \eqn{\lambda_i}, coefficient vectors \eqn{b_i}, and fixed penalty matrices \eqn{S_i}.
#'
#' The \strong{qREML} algorithm is typically much faster than REML or marginal ML using the full Laplace approximation method, but may be slightly less accurate regarding the estimation of the penalty strength parameters.
#'
#' Under the hood, \code{qreml} uses the R package \code{RTMB} for automatic differentiation in the inner optimisation.
#' The user has to specify the \strong{penalised negative log-likelihood function} \code{pnll} structured as dictated by \code{RTMB} and use the \code{\link{penalty}} function to compute the quadratic-form penalty inside the likelihood.
#' 
#' @seealso \code{\link{penalty}} to compute the penalty inside the likelihood function
#' 
#' @references Koslik, J. O. (2024). Efficient smoothness selection for nonparametric Markov-switching models via quasi restricted maximum likelihood. arXiv preprint arXiv:2411.11498.
#'
#' @param pnll penalised negative log-likelihood function that is structured as dictated by \code{RTMB} and uses the \code{\link{penalty}} function from \code{LaMa} to compute the penalty
#'
#' Needs to be a function of the named list of initial parameters \code{par} only.
#' @param par named list of initial parameters
#'
#' The random effects/ spline coefficients can be vectors or matrices, the latter summarising several random effects of the same structure, each one being a row in the matrix.
#' @param dat initial data list that contains the data used in the likelihood function, hyperparameters, and the \strong{initial penalty strength} vector
#'
#' If the initial penalty strength vector is \strong{not} called \code{lambda}, the name it has in \code{dat} needs to be specified using the \code{psname} argument below.
#' Its length needs to match the to the total number of random effects.
#' @param random vector of names of the random effects/ penalised parameters in \code{par}
#' 
#' \strong{Caution:} The ordering of \code{random} needs to match the order of the random effects passed to \code{\link{penalty}} inside the likelihood function.
#' @param map optional map argument, containing factor vectors to indicate parameter sharing or fixing.
#' 
#' Needs to be a named list for a subset of fixed effect parameters or penalty strength parameters. 
#' For example, if the model has four penalty strength parameters, \code{map[[psname]]} could be \code{factor(c(NA, 1, 1, 2))} to fix the first penalty strength parameter, estimate the second and third jointly, and estimate the fourth separately.
#' @param psname optional name given to the penalty strength parameter in \code{dat}. Defaults to \code{"lambda"}.
#' @param alpha optional hyperparamater for exponential smoothing of the penalty strengths.
#'
#' For larger values smoother convergence is to be expected but the algorithm may need more iterations.
#' @param smoothing optional scaling factor for the final penalty strength parameters
#' 
#' Increasing this beyond one will lead to a smoother final model. Can be an integer or a vector of length equal to the length of the penalty strength parameter.
#' @param maxiter maximum number of iterations in the outer optimisation over the penalty strength parameters.
#' @param tol Convergence tolerance for the penalty strength parameters.
#' @param control list of control parameters for \code{\link[stats:optim]{optim}} to use in the inner optimisation. Here, \code{optim} uses the \code{BFGS} method which cannot be changed.
#' 
#' We advise against changing the default values of \code{reltol} and \code{maxit} as this can decrease the accuracy of the Laplace approximation.
#' @param silent integer silencing level: 0 corresponds to full printing of inner and outer iterations, 1 to printing of outer iterations only, and 2 to no printing.
#' @param joint_unc logical, if \code{TRUE}, joint \code{RTMB} object is returned allowing for joint uncertainty quantification
#' @param saveall logical, if \code{TRUE}, then all model objects from each iteration are saved in the final model object.
#' # @param epsilon vector of two values specifying the cycling detection parameters. If the relative change of the new penalty strength to the previous one is larger than \code{epsilon[1]} but the change to the one before is smaller than \code{epsilon[2]}, the algorithm will average the two last values to prevent cycling.
#'
#' @return model object of class 'qremlModel'. This is a list containing:
#' \item{...}{everything that is reported inside \code{pnll} using \code{RTMB::REPORT()}. When using \code{forward}, \code{tpm_g}, etc., this may involve automatically reported objects.}
#' \item{obj}{\code{RTMB} AD object containing the final conditional model fit}
#' \item{psname}{final penalty strength parameter vector}
#' \item{all_psname}{list of all penalty strength parameter vectors over the iterations}
#' \item{par}{named estimated parameter list in the same structure as the initial \code{par}. Note that the name \code{par} is not fixed but depends on the original name of your \code{par} list.}
#' \item{relist_par}{function to convert the estimated parameter vector to the estimated parameter list. This is useful for uncertainty quantification based on sampling from a multivariate normal distribution.}
#' \item{par_vec}{estimated parameter vector}
#' \item{llk}{unpenalised log-likelihood at the optimum}
#' \item{n_fixpar}{number of fixed, i.e. unpenalised, parameters}
#' \item{edf}{overall effective number of parameters}
#' \item{all_edf}{list of effective number of parameters for each smooth}
#' \item{Hessian_condtional}{final Hessian of the conditional penalised fit}
#' \item{obj_joint}{if \code{joint_unc = TRUE}, joint \code{RTMB} object for joint uncertainty quantification in model and penalty parameters.}
#' 
#' @export
#'
#' @import RTMB
#'
#' @examples
#' # no example
qreml_old = function(pnll, # penalized negative log-likelihood function
                 par, # initial parameter list
                 dat, # initial dat object, currently needs to be called dat!
                 random, # names of parameters in par that are random effects/ penalized
                 map = NULL, # map for fixed effects
                 psname = "lambda", # name given to the psname parameter in dat
                 alpha = 0.25, # exponential smoothing parameter
                 smoothing = 1,
                 maxiter = 100, # maximum number of iterations
                 tol = 1e-4, # tolerance for convergence
                 control = list(reltol = 1e-10, maxit = 1000), # control list for inner optimization
                 silent = 1, # print level
                 joint_unc = TRUE, # should joint object be returned?
                 saveall = FALSE)# , # save all intermediate models?
                 #epsilon = c(1e-2, 1e-1)) # cycling detection parameters 
{
  
  # setting the argument name for par because later updated par is returned
  argname_par = as.character(substitute(par))
  argname_dat = as.character(substitute(dat))
  
  # number of random effects, each one can be a matrix where each row is a random effect, but then they have the same penalty structure
  n_re = length(random) 
  
  # list to save all model objects
  allmods = list() 
  
  # initial lambda locally
  # Define a global environment to store the captured names
  # penalty_metadata <- new.env(parent = emptyenv())
  # if(is.null(psname)){
  #   pnll(par) # call once to get the name of the lambda parameter
  #   psname = get("argname_lambda", envir = penalty_metadata)
  # }
  
  lambda = dat[[psname]]
  lambda0 = lambda # saving initial lambda so that fixed pars can always be refilled
  
  # creating the objective function as wrapper around pnll to pull lambda from local
  f = function(par){
    environment(pnll) = environment()
    
    # overloading assignment operators, currently necessary
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    getLambda = function(x) lambda
    
    dat[[psname]] = DataEval(getLambda, rep(advector(1), 0))
    
    # assigning dat to whatever it is called in pnll() (hopefully)
    assign(argname_dat, dat, envir = environment())
    
    pnll(par)
  }
  
  # creating the RTMB objective function
  if(silent %in% 0:1){
    cat("Creating AD function\n")
  } 
  
  ## mapping
  # map can contain fixed effects -> just passed to MakeADFun
  # and it can contain penalty strength parameters
  if(!is.null(map)){
    # check that no random effects are fixed
    if(any(names(map) %in% random)){
      msg <- "'map' cannot contain random effects or spline parameters"
      stop(msg)
    }
    # make factor if not
    if(!all(sapply(map, is.factor))){
      message("Converting map to factor")
    }
    map <- lapply(map, factor)
    # if there is mapping but no psname map, add psname map
    if(is.null(map[[psname]])){
      map[[psname]] = factor(seq_along(lambda))
    }
  } else {
    # initialises a list only having named element psname
    map[[psname]] = factor(seq_along(lambda))
  }
  # separate out psname map
  lambda_map <- map[[psname]]
  if(length(lambda_map) != length(lambda)){
    msg <- paste0("Length of map argument for ", psname, " has wrong length.")
    stop(msg)
  }
  
  # pop lambda_map from map list
  map <- map[names(map) != psname]
  # if the remaining map is now an empty list, set to NULL to work with MakeADFun
  if(length(map) == 0) map = NULL
  
  # deal with mapping of penalty strength parameters
  lambda_mapped = map_lambda(lambda, lambda_map)
  if(length(lambda_mapped) < 1){
    stop("No penalty paramters estimated as all are fixed.")
  }
  
  obj = MakeADFun(func = f, 
                  parameters = par, 
                  silent = TRUE,
                  map = map) # silent and replacing with own prints
  
  newpar = obj$par # saving initial parameter value as vector to initialize optimization in loop
  
  # own printing of maximum gradient component if silent = 0
  gradcounter = 1
  if(silent == 0){
    newgrad = function(par){
      gr = obj$gr(par)
      if(gradcounter %% 10 == 0){
        iter = gradcounter / 10
        cat("  iter:", iter, "-", "inner mgc:", max(abs(gr)), "\n")
      }
      gradcounter <<- gradcounter + 1
      gr
    }
  } else{
    newgrad = obj$gr
  }
  
  # prepwork
  mod0 = obj$report() # getting all necessary information from penalty report
  S = mod0$S # penalty matrix/ matrices in list format
  # S_dims = sapply(S, nrow)
  
  # finding the indices of the random effects to later index Hessian
  re_inds = list() 
  for(i in seq_len(n_re)){
    if(is.vector(par[[random[i]]])){
      re_dim = c(1, length(par[[random[i]]]))
    } else if(is.matrix(par[[random[i]]])){
      re_dim = dim(par[[random[i]]])
    } else{
      stop(paste0(random[i], " must be a vector or matrix"))
    }
    
    byrow = FALSE
    if(re_dim[1] == nrow(S[[i]])){
      byrow = TRUE
    }
    #re_dim = dim(as.matrix(par[[random[i]]]))
    
    # which(re_dim == nrow(S[[i]])) - 1
    
    # if(re_dim[2] == S_dims[i]){
    #   byrow = FALSE
    # } else{
    #   byrow = TRUE
    # }
    re_inds[[i]] = matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])#, byrow = byrow)
    if(byrow) {
      re_inds[[i]] = t(re_inds[[i]]) # if byrow, then transpose
    }
    # if(dim(re_inds[[i]])[2] == 1) re_inds[[i]] = t(re_inds[[i]]) # if only one column, then transpose
  }
  
  # get number of similar random effects for each distinct random effect (of same structure)
  re_lengths = sapply(re_inds, function(x) if (is.vector(x)) 1 else nrow(x))
  
  # initialize list of penalty strength parameters
  Lambdas = list()
  Lambdas[[1]] = reshape_lambda(re_lengths, lambda) # reshaping to match structure of random effects
  
  if(silent < 2){
    cat("Initialising with", paste0(psname, ":"), round(lambda, 3), "\n")
  }
  
  # computing rank deficiency for each penalty matrix to use in correction term
  m = numeric(length(S)) 
  for(i in seq_len(length(m))) {
    m[i] = nrow(S[[i]]) - Matrix::rankMatrix(S[[i]])
  } 
  
  # initialising convergence check index (initially for all lambdas)
  convInd <- seq_along(unlist(Lambdas[[1]]))
  
  ### updating algorithm
  # loop over outer iterations until convergence or maxiter
  for(k in seq_len(maxiter)){
    
    # fitting the model conditional on lambda: current local lambda will be pulled by f
    gradcounter = 1
    opt = stats::optim(newpar, obj$fn, newgrad, 
                       method = "BFGS", hessian = TRUE, # return hessian in the end
                       control = control)
    if(silent == 0){
      gr = obj$gr(opt$par)
      cat("  final inner mgc:", max(abs(gr)), "\n")
    }
    
    # setting new optimum par for next iteration
    newpar = opt$par 
    
    # reporting to extract penalties
    mod = obj$report() 
    
    # evaluating current Hessian
    # J = obj$he()
    J = opt$hessian
    
    # computing inverse Hessian
    J_inv = MASS::ginv(J) 
    
    # saving entire model object
    if(saveall){
      allmods[[k]] = mod
    }
    
    ## updating all lambdas
    lambdas_k = list() # temporary lambda list
    
    # looping over distinct random effects (matrices)
    edoF = rep(NA, length(lambda0)) # initialise edoF vector
    pens = rep(NA, length(lambda0)) # initialise penalty vector
    l = 1 # counter for lambda
    
    for(i in 1:n_re){
      # initializing lambda vector for i-th random effect
      lambdas_k[[i]] = numeric(nrow(re_inds[[i]]))
      
      # looping over similar random effects
      for(j in 1:nrow(re_inds[[i]])){
        idx = re_inds[[i]][j,] # indices of this random effect
        
        # effective degrees of freedom for this random effect: J^-1_p J
        edoF[l] = nrow(S[[i]]) - m[i] - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
        
        # penalty
        pens[l] = mod$Pen[[i]][j]
        
        l = l+1
        
      #   # effective degrees of freedom for this random effect: J^-1_p J
      #   edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
      #   
      #   # calculating new lambda based on updating rule
      #   lambda_new = as.numeric((edoF - m[i]) / mod$Pen[[i]][j]) # m is correction if S_i does not have full rank
      #   
      #   # potentially smoothing new lambda
      #   lambdas_k[[i]][j] = (1-alpha) * lambda_new + alpha * Lambdas[[k]][[i]][j]
      #   
      #   # check for cycling behaviour
      #   if(k > 2){
      #     if(abs((lambdas_k[[i]][j] - Lambdas[[k-1]][[i]][j]) / Lambdas[[k-1]][[i]][j]) < epsilon[1] & # change to lambda_t-2 is small
      #        abs((lambdas_k[[i]][j] - Lambdas[[k]][[i]][j]) / Lambdas[[k]][[i]][j]) > epsilon[2]) # but change to lambda_t-1 is large
      #     {
      #       cat("Cycling detected - averaging for faster convergence\n")
      #       # replacing with mean to prevent cycling
      #       lambdas_k[[i]][j] = (lambdas_k[[i]][j] + Lambdas[[k]][[i]][j]) / 2 
      #     }
      #   }
      }
    }
    
    # now loop over actual lambda_mapped to update
    outer_gr = numeric(length(lambda_mapped))
    for(i in seq_along(lambda_mapped)){
      this_level = levels(lambda_map)[i]
      this_ind = which(lambda_map == this_level)
      
      this_edoF = sum(edoF[this_ind])
      this_pen = sum(pens[this_ind])
      
      lambda_new = this_edoF / this_pen
      
      # smoothing lambda
      lambda_mapped[i] = (1-alpha) * lambda_new + alpha * lambda_mapped[i]
      
      # gradient
      outer_gr[i] = -0.5 * this_pen + 1/(2*lambda_mapped[i]) * this_edoF
    }
    
    # potentially set lambdas to "working infinity"
    lambda_mapped[which(lambda_mapped > 1e8)] <- 1e8
    
    # unmap lambda
    lambdas_k = unmap_lambda(lambda_mapped, lambda_map, lambda0)
    
    # minimum of zero for penalty strengths
    lambdas_k[which(lambdas_k < 0)] = 0
    
    # assigning new lambda to global list
    Lambdas[[k+1]] = utils::relist(lambdas_k, Lambdas[[1]])
    
    # updating lambda vector locally for next iteration
    lambda = lambdas_k
    
    # old length of convergence check indices
    oldlength <- length(convInd)
    
    if(k > 2){ # after 2 iterations, check whether any lambda > 1e5 and exclude from check
      convInd = which(lambda <= 1e6)
    }
    
    if(silent < 2){
      cat("outer", k, "-", paste0(psname, ":"), round(lambda, 3), "\n")
      if(silent == 0){
        cat("outer mgc:", max(abs(outer_gr)), "\n")
      }
      
      # print only if something changes
      if(length(convInd) != oldlength & length(seq_along(lambda)[-convInd]) > 0){
        cat(psname, seq_along(lambda)[-convInd], "excluded from convergence check (> 1e6)", "\n")
      }
    }
    
    # convergence check
    # if(all(abs(lambda - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]])) < tol)){
    if(max(abs(
      (lambda - unlist(Lambdas[[k]]))[convInd] / unlist(Lambdas[[k]])[convInd]
    )) < tol){
      
      if(silent < 2){
        cat("Converged\n")
      }
      break
    }
    
    if(k == maxiter){
      cat("No convergence\n")
      warning("No convergence\n")
    } 
  }
  
  # final model fit
  lambda = lambda * smoothing # scaling lambda by smoothing factor
  
  if(silent < 2){
    if(any(smoothing != 1)){
      cat("Smoothing factor:", smoothing, "\n")
    }
    cat("Final model fit with", paste0(psname, ":"), round(lambda, 3), "\n")
  }
  
  # fitting the model conditional on lambda: current local lambda will be pulled by f
  gradcounter = 1
  opt = stats::optim(newpar, obj$fn, newgrad, 
                     method = "BFGS", hessian = TRUE, # return hessian in the end
                     control = control)
  if(silent == 0){
    gr = obj$gr(opt$par)
    cat("  final inner mgc:", max(abs(gr)), "\n")
  }
  
  # setting new optimum par for next iteration
  newpar = opt$par 
  
  # reporting to extract penalties
  mod = obj$report() 
  
  # evaluating current Hessian
  # J = obj$he()
  J = opt$hessian
  
  # computing inverse Hessian
  J_inv = MASS::ginv(J) 
  
  # saving entire model object
  if(saveall){
    allmods[[k+1]] = mod
  }
  
  #############################################
  
  # assign RTMB obj to return object
  mod$obj <- obj
  
  # if all intermediate models should be returned, assign
  if(saveall) {
    mod$allmods = allmods
  }
  
  # assign final lambda to return object
  mod[[psname]] = lambda
  
  # assigning all lambdas to return object
  mod[[paste0("all_", psname)]] = Lambdas
  
  # calculating unpenalized log-likelihood at final parameter values
  lambda = rep(0, length(lambda))
  dat[[psname]] = lambda
  
  # format parameter to list
  # skeleton = utils::as.relistable(par)
  # parlist = utils::relist(opt$par, skeleton)
  parlist = obj$env$parList(opt$par)
  mod[[argname_par]] = parlist # and assing to return object
  
  mod[[paste0("relist_", argname_par)]] = obj$env$parList
  
  # assign estimated parameter as vector
  mod[[paste0(argname_par, "_vec")]] = opt$par
  
  # assign log-likelihood at optimum to return object
  mod$llk = -pnll(parlist)
  
  ## calculating effective degrees of freedom for final model
  mod$edf = list()
  for(i in 1:n_re){
    edoF_i = numeric(nrow(re_inds[[i]]))
    for(j in 1:nrow(re_inds[[i]])){
      idx = re_inds[[i]][j,]
      edoF_i[j] = edoF = nrow(S[[i]]) - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]]))
    }
    mod$edf[[i]] = edoF_i
  }
  
  # number of fixed parameters
  # mod$n_fixpar = length(unlist(par[!(names(par) %in% random)]))
  mod$n_fixpar = length(opt$par)
  
  # assing conditinoal Hessian
  mod$Hessian_conditional = J
  
  # removing penalty list from model object
  mod = mod[names(mod) != "Pen"] 
  
  
  if(joint_unc){
    ### constructing joint object
    parlist$loglambda = log(mod[[psname]])
    
    # finding the number of similar random effects for each random effect
    # indvec = rep(1:n_re, times = re_lengths)
    
    # computing log determinants
    logdetS = numeric(length(S))
    for(i in 1:length(S)){
      logdetS[i] = gdeterminant(S[[i]])
    }
    
    ## defining joint negative log-likelihood
    jnll = function(par) {
      
      environment(pnll) = environment()
      
      # overloading assignment operators, currently necessary
      "[<-" <- ADoverload("[<-") 
      "c" <- ADoverload("c")
      "diag<-" <- ADoverload("diag<-")
      
      dat[[psname]] = exp(par$loglambda)
      
      l_p = -pnll(par[names(par) != "loglambda"])
      
      ## computing additive constants (missing from only penalized likelihood)
      const = 0
      for(i in 1:n_re){
        for(j in 1:nrow(re_inds[[i]])){
          k = length(re_inds[[i]][j,])
          
          if(i == 1){
            loglam = par$loglambda[j]
          } else{
            loglam = par$loglambda[re_lengths[i-1] + j]
          }
          
          const = const - k * log(2*pi) + k * loglam + logdetS[i]
        }
      }
      
      l_joint = l_p + 0.5 * const
      -l_joint
    }
    
    if(is.null(map)){
      map = list(loglambda = lambda_map)
    } else{
      map$loglambda = lambda_map
    }
    
    # creating joint AD object
    obj_joint = MakeADFun(jnll, parlist,
                          random = names(par)[names(par) != "loglambda"], # REML, everything random except lambda
                          map = map)
    
    # assigning object to return object
    mod$obj_joint = obj_joint
  }
  
  class(mod) = "qremlModel"
  return(mod)
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
#' par = list(logmu = log(c(0.3, 1)), # step mean
#'            logsigma = log(c(0.2, 0.7)), # step sd
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
  
  # RTMB::REPORT(re_coef)
  
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
        # compute quadratic form
        subthispen[k] = thislambda[k] * (t(this_tp[j,]) %*% thisS[[k]] %*% this_tp[j,])
      }
      # add to overall penalty
      pen = pen + sum(subthispen)
      
      thispen[[j]] <- subthispen # in this case we have multiple penalties we need for the update
    }
    
    Pen[[i]] <- thispen
  }
  
  RTMB::REPORT(Pen) # Report the penalty list for qreml update
  pen = 0.5 * pen
  RTMB::REPORT(pen)
  pen
}

#' Quasi restricted maximum likelihood (qREML) algorithm for models with penalised splines or simple i.i.d. random effects
#'
#' @description
#' This algorithm can be used very flexibly to fit statistical models that involve \strong{penalised splines} or simple \strong{i.i.d. random effects}, i.e. that have penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
#' with smoothing parameters \eqn{\lambda_i}, coefficient vectors \eqn{b_i}, and fixed penalty matrices \eqn{S_i}.
#'
#' The \strong{qREML} algorithm is typically much faster than REML or marginal ML using the full Laplace approximation method, but may be slightly less accurate regarding the estimation of the penalty strength parameters.
#'
#' Under the hood, \code{qreml} uses the R package \code{RTMB} for automatic differentiation in the inner optimisation.
#' The user has to specify the \strong{penalised negative log-likelihood function} \code{pnll} structured as dictated by \code{RTMB} and use the \code{\link{penalty}} function to compute the quadratic-form penalty inside the likelihood.
#' 
#' @seealso \code{\link{penalty}} and \code{\link{penalty2}} to compute the penalty inside the likelihood function
#' 
#' @references Koslik, J. O. (2024). Efficient smoothness selection for nonparametric Markov-switching models via quasi restricted maximum likelihood. arXiv preprint arXiv:2411.11498.
#'
#' @param pnll penalised negative log-likelihood function that is structured as dictated by \code{RTMB} and uses the \code{\link{penalty}} function from \code{LaMa} to compute the penalty
#'
#' Needs to be a function of the named list of initial parameters \code{par} only.
#' @param par named list of initial parameters
#'
#' The random effects/ spline coefficients can be vectors or matrices, the latter summarising several random effects of the same structure, each one being a row in the matrix.
#' @param dat initial data list that contains the data used in the likelihood function, hyperparameters, and the \strong{initial penalty strength} vector
#'
#' If the initial penalty strength vector is \strong{not} called \code{lambda}, the name it has in \code{dat} needs to be specified using the \code{psname} argument below.
#' Its length needs to match the to the total number of random effects.
#' @param random vector of names of the random effects/ penalised parameters in \code{par}
#' 
#' \strong{Caution:} The ordering of \code{random} needs to match the order of the random effects passed to \code{\link{penalty}} inside the likelihood function.
#' @param map optional map argument, containing factor vectors to indicate parameter sharing or fixing.
#' 
#' Needs to be a named list for a subset of fixed effect parameters or penalty strength parameters. 
#' For example, if the model has four penalty strength parameters, \code{map[[psname]]} could be \code{factor(c(NA, 1, 1, 2))} to fix the first penalty strength parameter, estimate the second and third jointly, and estimate the fourth separately.
#' @param silent integer silencing level: 0 corresponds to full printing of inner and outer iterations, 1 to printing of outer iterations only, and 2 to no printing.
#' @param psname optional name given to the penalty strength parameter in \code{dat}. Defaults to \code{"lambda"}.
#' @param alpha optional hyperparamater for exponential smoothing of the penalty strengths.
#'
#' For larger values smoother convergence is to be expected but the algorithm may need more iterations.
#' @param smoothing optional scaling factor for the final penalty strength parameters
#' 
#' Increasing this beyond one will lead to a smoother final model. Can be an integer or a vector of length equal to the length of the penalty strength parameter.
#' @param maxiter maximum number of iterations in the outer optimisation over the penalty strength parameters.
#' @param tol Convergence tolerance for the penalty strength parameters.
#' @param method optimisation method to be used by \code{\link[stats:optim]{optim}}. Defaults to \code{"BFGS"}, but might be changed to \code{"L-BFGS-B"} for high-dimensional settings.
#' @param control list of control parameters for \code{\link[stats:optim]{optim}} to use in the inner optimisation. Here, \code{optim} uses the \code{BFGS} method which cannot be changed.
#' 
#' We advise against changing the default values of \code{reltol} and \code{maxit} as this can decrease the accuracy of the Laplace approximation.
#' @param method optimisation method to be used by \code{\link[stats:optim]{optim}}. Defaults to \code{"BFGS"}, but might be changed to \code{"L-BFGS-B"} for high-dimensional settings.
#' @param conv_crit character, convergence criterion for the penalty strength parameters. Can be \code{"relchange"} (default) or \code{"gradient"}.
#' @param joint_unc logical, if \code{TRUE}, joint \code{RTMB} object is returned allowing for joint uncertainty quantification
#' @param spHess logical, if \code{TRUE}, sparse AD Hessian is used in each outer iteration. If your Hessian is large and sparse (many cross derivatives are 0), this will speed up the computations a lot. 
#' If your Hessian is dense, this will slow down the computations slightly and might require significantly more memory.
#' @param saveall logical, if \code{TRUE}, then all model objects from each iteration are saved in the final model object.
#'
#' @return model object of class 'qremlModel'. This is a list containing:
#' \item{...}{everything that is reported inside \code{pnll} using \code{RTMB::REPORT()}. When using \code{forward}, \code{tpm_g}, etc., this may involve automatically reported objects.}
#' \item{obj}{\code{RTMB} AD object containing the final conditional model fit}
#' \item{psname}{final penalty strength parameter vector}
#' \item{all_psname}{list of all penalty strength parameter vectors over the iterations}
#' \item{par}{named estimated parameter list in the same structure as the initial \code{par}. Note that the name \code{par} is not fixed but depends on the original name of your \code{par} list.}
#' \item{relist_par}{function to convert the estimated parameter vector to the estimated parameter list. This is useful for uncertainty quantification based on sampling from a multivariate normal distribution.}
#' \item{par_vec}{estimated parameter vector}
#' \item{llk}{unpenalised log-likelihood at the optimum}
#' \item{n_fixpar}{number of fixed, i.e. unpenalised, parameters}
#' \item{edf}{overall effective number of parameters}
#' \item{all_edf}{list of effective number of parameters for each smooth}
#' \item{Hessian_condtional}{final Hessian of the conditional penalised fit}
#' \item{obj_joint}{if \code{joint_unc = TRUE}, joint \code{RTMB} object for joint uncertainty quantification in model and penalty parameters.}
#'
#' @export
#'
#' @import RTMB
#'
#' @examples
#' data = trex[1:1000,] # subset
#'
#' # initial parameter list
#' par = list(logmu = log(c(0.3, 1)), # step mean
#'            logsigma = log(c(0.2, 0.7)), # step sd
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
#'   Gamma = tpm_g(Z, cbind(beta0, betaspline), ad = TRUE) # transition probabilities
#'   delta = stationary_p(Gamma, t = 1, ad = TRUE) # initial distribution
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
qreml <- function(pnll, # penalized negative log-likelihood function
                  par, # initial parameter list
                  dat, # initial dat object, currently needs to be called dat!
                  random, # names of parameters in par that are random effects/ penalized
                  map = NULL, # map for fixed effects
                  silent = 1, # print level
                  psname = "lambda", # name given to the psname parameter in dat
                  alpha = 0.3, # exponential smoothing parameter
                  smoothing = 1,
                  maxiter = 100, # maximum number of iterations
                  tol = 1e-4, # tolerance for convergence
                  method = "BFGS", # optimization method used by optim
                  control = list(), # control list for inner optimization
                  conv_crit = "relchange",
                  spHess = FALSE, # use sparse Hessian instead of finite diff gradient
                  joint_unc = FALSE, # should joint object be returned?
                  saveall = FALSE # save all intermediate models?
                  )
{
  ### input checking arguments
  if(!is.function(pnll)){
    stop("pnll needs to be a function")
  }
  if(!conv_crit %in% c("gradient", "relchange")){
    stop("'conv_crit' needs to be either 'gradient' or 'relchange'")
  }
  if(!is.list(par)){
    stop("'par' needs to be a named list")
  }
  if(!is.list(dat)){
    stop("'dat' needs to be a named list")
  }
  if(!psname %in% names(dat)){
    stop(paste0("'dat' needs to contain a vector called '", psname, "' with initial penalty strengths"))
  }
  if(!is.character(random) || length(random) < 1){
    stop("'random' needs to be a character vector of names of random effects in 'par'")
  }
  if(!is.null(map) && !is.list(map)){
    stop("'map' needs to be a named list of factors for fixed effects or penalty strength parameters")
  }
  
  # setting the argument name for par because later updated par is returned
  argname_par <- as.character(substitute(par))
  argname_dat <- as.character(substitute(dat))
  
  # number of random effects, each one can be a matrix where each row is a random effect, but then they have the same penalty structure
  n_re <- length(random) 
  
  # list to save all model objects
  allmods <- list() 
  
  # initialising penalty strength lambda
  lambda <- dat[[psname]]
  lambda0 <- lambda # saving initial lambda so that fixed parts can always be refilled even when 'lambda' is changed
  
  # creating the objective function as wrapper around pnll to pull lambda from local
  f <- function(par){
    # setting the environment
    environment(pnll) = environment()
    
    # overloading assignment operators, currently necessary
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    # defining function that grabs lambda
    getLambda <- function(x) lambda
    # grab lambda from outside
    dat[[psname]] <- DataEval(getLambda, rep(advector(1), 0))
    
    # assigning dat to whatever it is called in pnll() (hopefully)
    assign(argname_dat, dat, envir = environment())
    
    pnll(par)
  }
  
  ## mapping
  # map can contain fixed effects -> just passed to MakeADFun
  # and it can contain penalty strength parameters
  if(!is.null(map)){
    # check that no random effects are fixed
    if(any(names(map) %in% random)){
      msg <- "'map' cannot contain random effects or spline parameters"
      stop(msg)
    }
    map <- lapply(map, factor)
  }
  # if there is mapping but no psname map, add psname map
  if(is.null(map[[psname]])){
    map[[psname]] = factor(seq_along(lambda))
  }
  # separate out psname map
  lambda_map <- map[[psname]]
  if(length(lambda_map) != length(lambda)){
    msg <- paste0("Length of map argument for ", psname, " has wrong length.")
    stop(msg)
  }
  
  # pop lambda_map from map list
  map <- map[names(map) != psname]
  # if the remaining map is now an empty list, set to NULL to work with MakeADFun
  if(length(map) == 0) map = NULL
  
  # deal with mapping of penalty strength parameters
  lambda_mapped = map_lambda(lambda, lambda_map)
  if(length(lambda_mapped) < 1){
    message("No penalty parameters will be estimated as all are fixed.")
    maxiter <- 1
  }
  Lambda_mapped <- matrix(lambda_mapped, nrow = 1, ncol = length(lambda_mapped))
  
  # creating the RTMB objective function
  if(silent %in% 0:1) {
    message("Creating AD function")
  }
  obj <- MakeADFun(func = f, 
                   parameters = par, 
                   silent = TRUE,
                   map = map) # silent and replacing with own prints
  
  newpar <- obj$par # saving initial parameter value as vector to initialize optimization in loop
  
  # use sparse Hessian?
  if(spHess) {
    Tape <- RTMB::GetTape(obj, name = "ADFun") # get the Tape
    if(silent < 2) message("Constructing sparse Hessian")
    obj$spHess <- Tape$jacfun(sparse = TRUE)$jacfun(sparse = TRUE) # construct sparse Hessian function from Tape
    rm(Tape) # removing Tape to save memory
    # obj$gr <- function(p) as.matrix(spGrad(p))
  }
  
  # gradient printing
  counter_env <- new.env() # create environment to hold a counter
  counter_env$count <- 0 # initialise with zero
  if(silent == 0){
    ctREPORT <- 10 # by default, report every 10 calls
    if(!is.null(control$REPORT)){
      ctREPORT <- control$REPORT # if report is changed, use that
      control$REPORT <- NULL # remove REPORT from control to avoid problems with optim
    }
    
    newgrad <- function(par){
      counter_env$count <- counter_env$count + 1
      ct <- counter_env$count
      gr <- obj$gr(par)
      if(ct %% ctREPORT == 0) cat("iter", ct, "- inner mgc:", round(max(abs(gr)), 5), "\n")
      gr
    }
  } else{
    newgrad <- obj$gr
  }
  
  # prepwork -> running reporting to get necessary quantities
  mod0 <- obj$report() # getting all necessary information from penalty report
  S <- mod0$S # penalty matrix/ matrices in list format
  
  # finding the indices of the random effects to later index Hessian
  re_inds <- list() 
  for(i in seq_len(n_re)){
    if(is.vector(par[[random[i]]])){
      re_dim <- c(1, length(par[[random[i]]]))
    } else if(is.matrix(par[[random[i]]])){
      re_dim <- dim(par[[random[i]]])
    } else{
      stop(paste0(random[i], " must be a vector or matrix"))
    }
    
    byrow <- FALSE
    if(is.matrix(S[[i]])){ # one penalty matrix
      if(re_dim[1] == nrow(S[[i]])){
        byrow <- TRUE
      }
    } else if(is.list(S[[i]])){ # multiple penalty matrices
      if(re_dim[1] == nrow(S[[i]][[1]])){
        byrow <- TRUE
      }
    }
    
    re_inds[[i]] <- matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])
    if(byrow) re_inds[[i]] <- t(re_inds[[i]]) # if byrow, then transpose
  }
  
  ## find how many penalty strength pars are needed for ecah random effect
  # 1: univariate smooth
  # >1: tensorproduct
  n_penalties <- sapply(S, function(x) if(is.matrix(x)) 1 else length(x))
  
  ## Compte indices of simple univariate smooths and of tensorproduct smooths
  simple_ind <- which(n_penalties == 1)
  tp_ind <- which(n_penalties > 1)
  
  # get number of similar random effects for each distinct random effect (of same structure)
  re_lengths = sapply(re_inds, function(x) if (is.vector(x)) 1 else nrow(x))
  
  ## total number of lambdas for each random effect with one penalty matrix/list
  lambda_lengths <- n_penalties * re_lengths
  if(length(lambda) != sum(lambda_lengths)){
    msg <- paste0("Length of '", psname, "' does not match the number of penalty strength parameters needed")
    stop(msg)
  }
  
  # initialize list of penalty strength parameters
  Lambdas <- list()
  Lambdas[[1]] <- reshape_lambda(lambda_lengths, lambda) # reshaping to match structure of random effects
  
  # naming lambdas better:
  # simple smooths: smooth_name.1, ..., smooth_name.re_lengths[[i]]
  for(ind in seq_along(simple_ind)){
    names(Lambdas[[1]][simple_ind][[ind]]) <- seq_along(Lambdas[[1]][simple_ind][[ind]])
  }
  # tensorproducts: same but additionally append margin name for clarity
  for(ind in seq_along(tp_ind)){
    margin_names <- names(S[[tp_ind[ind]]])
    names(Lambdas[[1]][tp_ind][[ind]]) <- paste0(rep(1:re_lengths[tp_ind[ind]], each = length(margin_names)),".",
                                                 rep(margin_names, re_lengths[tp_ind[ind]]))
  }
  lambda_names <- names(unlist(Lambdas[[1]]))
  
  if(silent < 2) message("Initialising with ", psname, ": ", paste(round(lambda, 3), collapse = " "))
  
  # Computing ranks of penalty matrices for simple_ind
  ranks <- sapply(S, function(x) if(is.matrix(x)) Matrix::rankMatrix(x) else NA)
  
  # locally define function to construct full penalty matrix from lambdas
  build_bigS <- function(lambdas) {
    bigS <- matrix(0, length(newpar), length(newpar))
    for(i in seq_len(n_re)){
      for(j in seq_len(nrow(re_inds[[i]]))){
        idx <- re_inds[[i]][j,]
        if(i %in% simple_ind){ # if simple smooth: just lambda_i * S_i
          bigS[idx, idx] <- lambdas[[i]][j] * S[[i]]
        } else { # if tensor product, we have a sum at these indices
          n_pen <- length(S[[i]])
          for(pen in 1:n_pen){ 
            bigS[idx, idx] <- bigS[idx, idx] + lambdas[[i]][(j-1) * n_pen + pen] * S[[i]][[pen]]
          }
        }
      }
    }
    bigS
  }
  
  # define restricted likelihood function
  restr_llk <- function(lp_opt, bigS, J) {
    lp_opt + gdeterminant(bigS) / 2 - gdeterminant(J) / 2
  }
  llk_r <- numeric(maxiter) # restricted likelihood vector
  
  # initialising convergence check index (initially for all lambdas)
  convInd <- seq_along(lambda_mapped)
  convInd_unmapped <- seq_along(lambda) # for unmapped lambdas
  
  # controlling optim printing
  ctl <- list(maxit = 1000)
  ctl[names(control)] <- control # overwriting with user-provided control parameters
  # if(silent == 0) ctl$trace = 1 else ctl$trace = 0 # setting trace to 1 if silent == 0, otherwise 0
  if(method == "BFGS") ctl$reltol <- 1e-10
  if(method == "L-BFGS-B") ctl$maxit <- 5000 # L-BFGS-B takes smaller steps
  
  ### updating algorithm
  # loop over outer iterations until convergence or maxiter
  for(k in seq_len(maxiter)){
    
    # set inner gradient counter to zero
    counter_env$count <- 0
    
    # fitting the model conditional on lambda: current local lambda will be pulled by f
    if(silent == 0) cat("\nInner optimisation:", "\n")
    opt <- stats::optim(newpar, obj$fn, newgrad, 
                        method = method,
                        control = ctl)
    
    gr <- obj$gr(opt$par)
    if(silent == 0){
      cat("iter", counter_env$count, "- inner mgc:", round(max(abs(gr)), 5), "\n")
    }
    
    # evaluating current penalised Hessian
    if(silent == 0) cat("evaluating Hessian...\n")
    if(spHess) {
      J <- obj$spHess(opt$par)
    } else{
      J <- stats::optimHess(opt$par, obj$fn, obj$gr)
      J <- (J + t(J))/2 # force symmetric
    }
    
    # build big penalty matrix from current lambdas
    bigS <- build_bigS(Lambdas[[k]])
    bigS <- (bigS + t(bigS)) / 2 # force symmetric 
    
    H <- J + bigS # Hessian = J + S_lambda
    # H <- (H + t(H)) / 2 # force symmetric Hessian
    
    # R <- tryCatch(chol(H), error = function(e) NULL)
    # if (is.null(R)) {
    #   if(silent == 0) cat("stabilising Hessian for inversion\n")
    #   eps <- 1e-8 * mean(diag(H))
    #   H <- H + diag(eps, nrow(H))
    #   R <- chol(H)
    # }
    
    H_inv <- safe_chol_inv(H) # chol2inv(R)
    
    # rebuild penalised Hessin pd for inversion
    J_pd <- H - bigS
    
    # # check if positive definite
    # if(!is.positive.definite(H)) {
    #   if(silent == 0) cat("replacing Hessian with nearest PD\n")
    #   H <- nearPD(H)$mat # if not, find nearest PD matrix
    # }
    
    # inverting current Hessian
    # try Cholesky
    # R <- tryCatch(chol(J_pd), error = function(e) NULL)
    # 
    # if (is.null(R)) {
    #   # not PD -> jitter
    #   eps <- 1e-8 * mean(diag(J_pd))
    #   J_pd <- J_pd + diag(eps, nrow(J_pd))
    #   R <- chol(J_pd)  # now must succeed
    # }
    
    # compute inverse
    J_inv <- safe_chol_inv(J_pd) # chol2inv(R)

    # J_inv <- tryCatch(solve(J_pd), error = function(e) NULL)
    # if(is.null(J_inv)) J_inv <- MASS::ginv(J_pd) # if problem, pseudo-inverse
    
    # setting new optimum par for next iteration
    newpar <- opt$par 
    
    # reporting to extract penalties
    mod <- obj$report() 
    
    # saving entire model object
    if(saveall){
      allmods[[k]] <- mod
    }
    
    ## calculating restricted likelihood
    llk_r[k] <- restr_llk(-opt$value, bigS, J)
    
    ### Updating all lambdas ###
    
    # looping over distinct random effects (matrices)
    edoF <- rep(NA, length(lambda0)) # initialise edoF vector
    pens <- rep(NA, length(lambda0)) # initialise penalty vector
    l <- 1 # counter for lambda vector
    
    # Loop over random effects (list entries)
    for(i in 1:n_re){
      
      # simple random effects with one smoothing parameter
      if(i %in% simple_ind){
        
        # looping over similar random effects (rows of re_coefs[[i]])
        for(j in 1:nrow(re_inds[[i]])){
          # indices of this random effect
          idx <- re_inds[[i]][j,]
          # effective degrees of freedom for this random effect
          edoF[l] <- ranks[i] - Lambdas[[k]][[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
          # quadratic penalty: b^t S b
          pens[l] <- mod$Pen[[i]][j]
          l <- l+1
        }
        
      } else if(i %in% tp_ind){ # more complicated tensorproduct random effects with multiple smoothing parameters
        
        # how many penalty matrices?
        n_pen <- length(S[[i]])
        
        # looping over similar random effects (rows of re_coefs[[i]])
        for(j in 1:nrow(re_inds[[i]])){
          # indices of this random effect
          idx <- re_inds[[i]][j,]
          # extracting old penalty strengths
          oldlambda <- Lambdas[[k]][[i]][(j-1) * n_pen + 1:n_pen]
          
          # effective degrees of freedom for this random effect
          # calculate (lambda_1* S_1 + ... + lambda_{n_pen} S_{n_pen})^-1
          thisS <- bigS[idx, idx] # extract submatrix of bigS for this random effect
          # thisS <- oldlambda[1] * S[[i]][[1]]
          # for(pen in 2:n_pen) thisS <- thisS + oldlambda[pen] * S[[i]][[pen]]
          thisS_inv <- MASS::ginv(thisS) # Moore-Penrose pseudo-inverse via SVD
          
          edoFs <- numeric(n_pen)
          for(pen in 1:n_pen){
            edoFs[pen] <- oldlambda[pen] * 
              (sum(rowSums(thisS_inv * S[[i]][[pen]])) - # tr(S^-1 S_j)
                 sum(rowSums(J_inv[idx, idx] * S[[i]][[pen]]))) # tr(J^-1 S_j)
          }
          edoF[l : (l + n_pen - 1)] <- edoFs
          
          # quadratic penalty: b^t S b, this is reported by penalty2()
          pens[l : (l + n_pen - 1)] <- mod$Pen[[i]][[j]]
          
          l <- l + n_pen
        }
      }
    }
    
    # now loop over actual lambda_mapped to update
    outer_gr <- numeric(length(lambda_mapped))
    for(i in seq_along(lambda_mapped)){
      this_level <- levels(lambda_map)[i]
      this_ind <- which(lambda_map == this_level)
      
      this_edoF <- sum(edoF[this_ind])
      this_pen <- sum(pens[this_ind])
      
      # compute new proposal
      lambda_new <- this_edoF / this_pen
      
      # smooth new proposal
      lambda_mapped[i] <- (1-alpha) * lambda_new + alpha * lambda_mapped[i]
      
      # gradient
      outer_gr[i] <- -0.5 * this_pen + 1 / (2 * lambda_mapped[i]) * this_edoF
    }
    
    # potentially set lambdas to "working infinity"
    lambda_mapped[which(lambda_mapped > 1e8)] <- 1e8
    
    # save current lambda_mapped
    Lambda_mapped <- rbind(Lambda_mapped, lambda_mapped)
    
    # unmap lambda
    lambdas_k <- unmap_lambda(lambda_mapped, lambda_map, lambda0)
    
    # minimum of zero for penalty strengths
    lambdas_k[which(lambdas_k < 0)] <- 0
    
    # assigning new lambda to global list
    Lambdas[[k+1]] <- Lambdas[[k]] # just to get the nice naming
    Lambdas[[k+1]] <- reshape_lambda(lambda_lengths, lambdas_k)
    
    # updating lambda vector locally for next iteration
    lambda <- lambdas_k
    
    # old length of convergence check indices
    oldlength <- length(convInd)
    
    if(k > 3){ # after 2 iterations, check whether any lambda > 1e5 and exclude from check
      convInd <- which(lambda_mapped <= 1e6)
      convInd_unmapped <- which(lambda <= 1e6) # indices of unmapped lambdas
    }
    
    mgc <- max(abs(outer_gr[convInd]))
    
    if(silent < 2){
      if(silent == 0) cat("\n")
      cat("outer", k, "-", paste0(psname, ":"), round(lambda, 3), "\n")
      if(silent == 0){
        cat("outer mgc:", mgc, "\n")
      }
      
      # print only if something changes
      if(length(convInd) != oldlength & length(lambda_mapped[-convInd]) > 0){
        cat(psname, seq_along(lambda)[-convInd], "excluded from convergence check (> 1e6)", "\n")
      }
    }
    
    #### convergence check ####
    if(conv_crit == "gradient"){
      if(k > 3 & (mgc < tol | opt$counts[2] < 3)) {
        if(silent < 2){
          message("Converged")
        }
        break
      }
    } else{
      # relative change of lambda
      rel_change <- abs((lambda - unlist(Lambdas[[k]])) / unlist(Lambdas[[k]]))
      
      if(k > 3 & (all(rel_change[convInd_unmapped] < tol) | opt$counts[2] < 3)) {
        if(silent < 2){
          message("Converged")
        }
        break
      }
    }
    
    if(k == maxiter){
      message("No convergence")
      warning("No convergence\n")
    } 
  }
  
  # final model fit
  lambda <- lambda * smoothing # scaling lambda by smoothing factor
  
  if(silent < 2){
    if(any(smoothing != 1)){
      message("Smoothing factor: ", paste(smoothing, collapse = " "))
    }
    if(silent == 0){
      message("\nFinal model fit with ", psname, ": ", paste(round(lambda, 3), collapse = " "))
    } else{
      message("Final model fit with ", psname, ": ", paste(round(lambda, 3), collapse = " "))
    }
  }
  
  # fitting the model conditional on final lambda
  opt <- stats::optim(newpar, obj$fn, newgrad, 
                      method = method, hessian = FALSE, # return hessian in the end
                      control = control)
  
  if(spHess) {
    J <- as.matrix(obj$spHess(opt$par))
  } else {
    J <- stats::optimHess(opt$par, obj$fn, obj$gr)
  }
  
  if(silent == 0){
    gr = obj$gr(opt$par)
    cat("final inner maximum gradient component:", round(max(abs(gr)), 5), "\n")
  }
  
  # reporting to extract penalties
  mod <- obj$report() 
  
  # save log likelihood at convergence
  pllk <- -opt$value # penalised
  llk <- pllk + mod$pen
  
  # evaluating current Hessian
  # J <- opt$hessian
  
  # computing inverse Hessian
  J_inv <- MASS::ginv(J) 
  
  # saving entire model object
  if(saveall){
    allmods[[k+1]] <- mod
  }
  
  #############################################
  
  # assign RTMB obj to return object
  mod$obj <- obj
  
  # if all intermediate models should be returned, assign
  if(saveall) {
    mod$allmods <- allmods
  }
  
  # assign gradient function
  mod$outer_gr <- function(x){
    lambda <- unmap_lambda(x, lambda_map, lambda0)
    Lambda <- reshape_lambda(lambda_lengths, lambda)
    
    environment(obj) = environment()
    
    # fitting the model conditional on lambda: current local lambda will be pulled by f
    if(silent == 0) cat("\nInner optimisation:", "\n")
    inner_opt <- stats::optim(newpar, obj$fn, newgrad,
                              method = method, hessian = TRUE, # return hessian in the end
                              control = control)
    thismod <- obj$report(inner_opt$par)
    J <- inner_opt$hessian
    J_inv <- MASS::ginv(J)
    # looping over distinct random effects (matrices)
    edoF <- rep(NA, length(lambda)) # initialise edoF vector
    pens <- rep(NA, length(lambda)) # initialise penalty vector
    l <- 1 # counter for lambda vector
    # Loop over random effects (list entries)
    for(i in 1:n_re){
      # simple random effects with one smoothing parameter
      if(i %in% simple_ind){
        for(j in 1:nrow(re_inds[[i]])){
          # indices of this random effect
          idx <- re_inds[[i]][j,]
          # effective degrees of freedom for this random effect: J^-1_p J
          edoF[l] <- ranks[i] - Lambda[[i]][j] * sum(rowSums(J_inv[idx, idx] * S[[i]])) # trace(J^-1 \lambda S)
          # quadratic penalty: b^t S b
          pens[l] <- thismod$Pen[[i]][j]
          l <- l+1
        }
      } else if(i %in% tp_ind){ # more complicated tensorproduct random effects with multiple smoothing parameters
        # how many penalty matrices?
        n_pen <- length(S[[i]])
        # looping over similar random effects (rows of re_coefs[[i]])
        for(j in 1:nrow(re_inds[[i]])){
          # indices of this random effect
          idx <- re_inds[[i]][j,]
          oldlambda <- Lambda[[i]][(j-1) * n_pen + 1:n_pen]
          # effective degrees of freedom for this random effect
          # calculate (lambda_1* S_1 + ... + lambda_{n_pen} S_{n_pen})^-1
          thisS <- oldlambda[1] * S[[i]][[1]]
          for(pen in 2:n_pen) thisS <- thisS + oldlambda[pen] * S[[i]][[pen]]
          thisS_inv <- MASS::ginv(thisS) # Moore-Penrose pseudo-inverse via SVD
          edoFs <- numeric(n_pen)
          for(pen in 1:n_pen){
            edoFs[pen] <- oldlambda[pen] *
              (sum(rowSums(thisS_inv * S[[i]][[pen]])) - # tr(S^-1 S_j)
                 sum(rowSums(J_inv[idx, idx] * S[[i]][[pen]]))) # tr(J^-1 S_j)
          }
          edoF[l : (l + n_pen - 1)] <- edoFs
          # quadratic penalty: b^t S b
          pens[l : (l + n_pen - 1)] <- thismod$Pen[[i]][[j]]
          l <- l + n_pen
        }
      }
    }
    # now loop over actual lambda_mapped to update
    outer_gr <- numeric(length(x))
    for(i in seq_along(x)){
      this_level <- levels(lambda_map)[i]
      this_ind <- which(lambda_map == this_level)
      this_edoF <- sum(edoF[this_ind])
      this_pen <- sum(pens[this_ind])
      # gradient
      outer_gr[i] <- -0.5 * this_pen + 1 / (2 * x[i]) * this_edoF
    }
    attr(outer_gr, "estimate") <- inner_opt$par
    outer_gr
  }
  environment(mod$outer_gr) <- environment()
  
  
  # assign final lambda to return object
  names(lambda) <- lambda_names
  mod[[psname]] <- lambda
  
  # assigning all lambdas to return object
  mod[[paste0("all_", psname)]] <- Lambdas
  
  # format parameter to list
  parlist <- obj$env$parList(opt$par)
  mod[[argname_par]] <- parlist # and assing to return object
  mod[[paste0("relist_", argname_par)]] <- obj$env$parList
  mod[[paste0("map_", psname)]] <- function(lambda) map_lambda(lambda, lambda_map)
  mod$psname <- psname
  mod$parname <- argname_par
  
  # assign estimated parameter as vector
  mod[[paste0(argname_par, "_vec")]] <- opt$par
  
  # assign log-likelihood at optimum to return object
  mod$llk <- llk
  
  # reassigning the correct lambda
  # lambda <- mod$lambda
  # dat[[psname]] <- lambda
  
  # number of fixed parameters
  mod$n_fixpar <- length(unlist(par[!(names(par) %in% random)]))
  
  ## compute effective degrees of freedom for each smooth (diag(J_p^-1 J))
  # building the entire model penalty matrix to compute J_0 = J_p - S
  # S_lambda = \sum_i lambda_i S_i padded out with zeros
  bigS <- build_bigS(Lambdas[[k+1]])
  
  leading_diag <- rowSums(J_inv * (J - bigS)) # computes diag(J_inv %*% (J - bigS)) more efficiently (only diagonal terms)
  Edfs <- Lambdas[[k+1]] # copy names from Lambdas if present
  for(i in seq_len(n_re)){
    if(i %in% tp_ind) Edfs[[i]] = numeric(nrow(re_inds[[i]])) # only one edf for each tensor product (not each margin)
    for(j in seq_len(nrow(re_inds[[i]]))){
        Edfs[[i]][j] = sum(leading_diag[re_inds[[i]][j,]]) # sum over the entries for each smooth
    }
  }
  mod$df <- mod$n_fixpar + sum(unlist(Edfs)) # total effective number of parameters
  mod$edf <- Edfs # seperated by smooth
  
  if(!is.null(mod$allprobs)){
    mod$nobs <- nrow(mod$allprobs) # number of observations
  }
  
  # assing conditinoal Hessian
  mod$Hessian_conditional <- J
  
  # assigning restriced likelihood
  mod$llk_restricted <- llk_r[1:k]
  
  # removing unnecessary elements that are only reported for qreml
  mod <- mod[!names(mod) %in% c("Pen", "pen", "S")] 
  
  if(length(tp_ind) == 0){ # only simple smooths, joint uncertainty possible
    if(joint_unc){
      ### constructing joint object
      parlist$loglambda <- log(mod[[psname]])
      
      # computing log determinants
      logdetS <- numeric(length(S))
      for(i in 1:length(S)){
        logdetS[i] <- gdeterminant(S[[i]])
      }
      
      ## defining joint negative log-likelihood
      jnll <- function(par) {
        
        environment(pnll) = environment()
        
        # overloading assignment operators, currently necessary
        "[<-" <- ADoverload("[<-") 
        "c" <- ADoverload("c")
        "diag<-" <- ADoverload("diag<-")
        
        dat[[psname]] <- exp(par$loglambda)
        
        l_p <- -pnll(par[names(par) != "loglambda"])
        
        ## computing additive constants (missing from only penalized likelihood)
        const <- 0
        for(i in 1:n_re){
          for(j in 1:nrow(re_inds[[i]])){
            k = length(re_inds[[i]][j,])
            
            if(i == 1){
              loglam <- par$loglambda[j]
            } else{
              loglam <- par$loglambda[re_lengths[i-1] + j]
            }
            
            const <- const - k * log(2*pi) + k * loglam + logdetS[i]
          }
        }
        
        l_joint <- l_p + 0.5 * const
        -l_joint
      }
      
      if(is.null(map)){
        map <- list(loglambda = lambda_map)
      } else{
        map$loglambda <- lambda_map
      }
      
      # creating joint AD object
      obj_joint <- MakeADFun(jnll, parlist,
                             random = names(par)[names(par) != "loglambda"], # REML, everything random except lambda
                             map = map)
      
      # assigning object to return object
      mod$obj_joint <- obj_joint
    }
  } 

  # clean up - may not be necessary
  gc()
  
  class(mod) = "qremlModel"
  return(mod)
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
  lambdas <- object[[object$psname]]
  lambda_names <- names(lambdas)
  if(is.null(lambda_names)){
    lambda_names <- paste0("lambda.", 1:length(lambdas))
  }
  for(i in seq_along(lambdas)){
    cat("\n", lambda_names[i], ": ", sep = "")
    cat(lambdas[[i]])
  }
  cat("\n")
  
  # Print additional user-specified objects, excluding unwanted ones
  excluded <- c("allprobs", "trackID", "type", "obj", "outer_gr", 
                paste0("all_", object$psname), "parname", object$parname, paste0("relist_", object$parname), 
                paste0("map_", object$psname), "psname", paste0(object$parname, "_vec"), 
                "edf", "Hessian_conditional", "obj_joint",
                "beta", "delta", "Gamma", "lambda", "llk", "n_fixpar", "df", "nobs",
                "llk_restricted", "allmods")
  
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
  
  psname <- mod$psname
  map_lambda <- mod[[paste0("map_", psname)]]
  outer_gr <- mod$outer_gr

  # map lambda
  lambda_mapped <- map_lambda(mod[[psname]])
  
  # map names of lambda
  lambda_names <- names(mod[[psname]])
  names(lambda_names) <- lambda_names
  mapped_names <- map_lambda(lambda_names)

  H <- - jacobian(outer_gr, lambda_mapped, method = "simple")
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
