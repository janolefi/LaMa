#' Get reported quantities from and \code{RTMB} object and return a \code{LaMaModel}
#' 
#' Having fitted a latent Markov model using automatic differentiation via \code{RTMB}, this function calls \code{RTMB}'s \code{obj$report()} and does some additional processing.
#' This then yields estimated parameters on their natural scale, allows for convenient calculation of \code{AIC} and \code{BIC}, state-decoding, and residual calculation.
#'
#' @param obj Object returned by \code{RTMB::MakeADFun()}
#'
#' @returns A model object of class "\code{LaMaModel}" containing a list with the reported quantities from the \code{RTMB} object, along with log-likelihood, number of parameters, and number of observations.
#' @export 
#'
#' @examples
#' data <- trex[1:200,]
#' 
#' # initial parameters and observations
#' par <- list(
#'   log_mu = log(c(0.3, 1)),      # initial means for step length (log-transformed)
#'   log_sigma = log(c(0.2, 0.7)), # initial sds for step length (log-transformed)
#'   eta = rep(-2, 2)              # initial t.p.m. parameters (on logit scale)
#' )    
#' dat <- list(
#'   step = data$step,             # hourly step lengths
#'   nStates = 2                   # number of hidden states
#' )
#' 
#' # likelihood function
#' nll <- function(par) {
#'   getAll(par, dat)
#'   REPORT(step)
#'   Gamma <- tpm(eta)
#'   delta <- stationary(Gamma)
#'   mu <- exp(log_mu); REPORT(mu)
#'   sigma <- exp(log_sigma); REPORT(sigma)
#'   allprobs <- matrix(1, length(step), nStates)
#'   ind <- which(!is.na(step))
#'   for(j in 1:nStates) {
#'     allprobs[ind,j] <- dgamma2(step[ind], mu[j], sigma[j]) 
#'   } 
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' # automatic differentiation and optimisation
#' obj <- MakeADFun(nll, par, silent = TRUE)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' 
#' ### reporting ###
#' mod <- report(obj)
#' 
#' # estimated quantities on natural scale
#' mod$mu
#' mod$sigma
#' mod$Gamma
#' 
#' # information criteria
#' AIC(mod)
#' BIC(mod)
#' 
#' # state decoding
#' states <- viterbi(mod = mod)   # global decoding
#' probs <- stateprobs(mod = mod) # local decoding
#' 
#' # residual calculation
#' pres <- pseudo_res(data$step, # observation sequence
#'                    "gamma2",  # distribution family
#'                    list(mean = mod$mu, sd = mod$sigma), # parameters for that family
#'                    mod = mod) # model object
report <- function(obj) {
  mod <- tryCatch(
    obj$report(),
    error = function(e) NULL
  )
  if(is.null(mod)) {
    stop("Does not seem to be an RTMB object.")
  } 
  
  # assign log-likelihood, number of parameters, and number of observations to the model object
  mod$ll <- -obj$fn()
  mod$df <- length(obj$par)
  mod$nobs <- tryCatch(
    nrow(mod$allprobs),
    error = function(e) NULL
    )
  
  mod$obs <- obj$env$obs
  
  class(mod) <- "LaMaModel"
  return(mod)
}

#' Extract log-likelihood from LaMaModel object
#' @param object A model fitted using RTMB and obtained via \code{report(obj)} of class "LaMaModel"
#' @param ... Additional arguments (not used)
#' @return An object of class "logLik"
#' @export
logLik.LaMaModel <- function(object, ...) {
  ll <- object$ll  # your stored log-likelihood
  df <- object$df # number of free parameters
  nobs <- object$nobs  # number of observations
  
  val <- as.numeric(ll)
  attr(val, "df") <- df
  attr(val, "nobs") <- nobs
  class(val) <- "logLik"
  val
}


# create environment to store information of distributions and parameters used in allprobs calculations
# forward() and forward_g() then access this environment and REPORT() the contents 
# this way we can capture all the parameters used in the observation distributions
# .obs_info <- new.env(parent = emptyenv())

# dist_reporter <- function(fun, family = NULL) {
#   if(is.null(family)) {
#     family <- substring(deparse(substitute(fun)), 2)
#   } 
#   
#   function(...) {
#     mc <- match.call(fun, expand.dots = FALSE)
#     stream <- deparse(mc[[2]])
#     
#     # match call to formal arguments
#     matched_call <- match.call(fun, call = mc)
#     matched_list <- as.list(matched_call)[-1]           # drop function name
#     matched_list <- matched_list[names(matched_list) != "x"]  # drop data argument
#     
#     # evaluate data and parameters
#     eval_args <- lapply(matched_list, eval, parent.frame())
#     obs_data <- eval(mc[[2]], parent.frame())  # get the actual observed data
#     
#     # initialize storage if first call
#     if (is.null(.obs_info$calls[[stream]])) {
#       .obs_info$calls[[stream]] <- list(
#         family = family,
#         obs    = obs_data,                   # store data once
#         pars   = lapply(eval_args, function(x) list())  # list per parameter
#       )
#     }
#     
#     # append each parameter value as a new element in the list
#     for (nm in names(eval_args)) {
#       .obs_info$calls[[stream]]$pars[[nm]][[length(.obs_info$calls[[stream]]$pars[[nm]]) + 1]] <- eval_args[[nm]]
#     }
#     
#     # call the original function
#     do.call(fun, list(...))
#   }
# }