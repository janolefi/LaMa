#' Build a \code{LaMaModel} from an optimised \code{RTMB} object
#' 
#' Having fitted a latent Markov model using automatic differentiation via \code{RTMB}, this function calls \code{RTMB}'s \code{obj$report()} and does some additional processing.
#' This then yields estimated parameters on their natural scale, allows for convenient calculation of \code{AIC} and \code{BIC}, state-decoding, and residual calculation.
#' \code{print()} gives a short overview of the model, \code{\link[=summary.LaMaModel]{summary()}} the estimates.
#' \code{report()} is the deprecated old name of \code{lama()}.
#'
#' @details
#' The number of observations (used by \code{BIC}) is the number of rows of the reported \code{allprobs} that are not all 1,
#' i.e. missing observations are not counted. This assumes that \code{\link{forward}} is called only once in the likelihood function.
#'
#' @param obj Optimised \code{RTMB} object
#'
#' @returns A model object of class "\code{LaMaModel}" containing a list with the reported quantities from the \code{RTMB} object, along estimated parameters (\code{par}),
#' the log-likelihood (\code{llk}), number of parameters (\code{df}) and observations (\code{nobs}), and the maximum gradient component (\code{mgc}) at the estimate.
#' It also contains \code{obj} itself (not a copy), so that functions like \code{\link{MCreport}} or \code{\link{scores}} can be called on the model object.
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
#' mod <- lama(obj)
#' mod          # short overview
#' summary(mod) # estimated state process and reported quantities
#'
#' # estimated parameters
#' mod$par
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
#'
#' # or from the model object alone, and plotted
#' plot(mod)
lama <- function(obj) {
  
  # get best parameter vector (potentially including random effects) from object
  p_hat <- tryCatch(
    obj$env$last.par.best,
    error = function(e) stop("Does not seem to be an RTMB object.")
  )
  mod <- tryCatch(
    obj$report(par = p_hat),
    error = function(e) stop("Does not seem to be an RTMB object.")
  )
  # Now trust that it is RTMB object
  
  mod$par <- obj$env$parList(par = p_hat)

  # fixed effects of the best parameter vector, at which fn and gr are evaluated
  random <- obj$env$random
  p_fix <- if(length(random)) p_hat[-random] else p_hat

  # assign log-likelihood, number of parameters, and number of observations to the model object
  mod$llk <- -obj$fn(p_fix)
  mod$df <- length(obj$par)
  mod$nobs <- count_nobs(mod$allprobs)
  mod$mgc <- max(abs(obj$gr(p_fix))) # maximum gradient component, to check convergence

  mod$obs <- obj$env$obs
  mod$obj <- obj # same object, not a copy
  
  class(mod) <- "LaMaModel"
  return(mod)
}

#' @rdname lama
#' @export
report <- function(obj) {
  .Deprecated("lama") # renamed, as it does much more than calling obj$report()
  lama(obj)
}

# Number of observations: rows of allprobs that are not all 1,
# as these are missing and contribute nothing to the likelihood
count_nobs <- function(allprobs){
  if(is.null(allprobs)) return(NULL)
  sum(rowSums(allprobs != 1) > 0)
}

# RTMB object from MakeADFun(), or the one stored by lama() or qreml()
get_obj <- function(x){
  if(inherits(x, c("LaMaModel", "qremlModel"))) x <- x$obj
  if(is.null(x$env)) stop("Expected an RTMB object from MakeADFun() or a model from lama() or qreml().")
  x
}

# estimated parameters, like obj$par: the final fit for qreml() models (whose objective changes
# with lambda, so last.par.best can be from an earlier iteration), else the best parameters of obj
est_par <- function(x){
  if(inherits(x, "qremlModel")) x[[paste0(x$parname, "_vec")]] else get_obj(x)$env$last.par.best
}

#' Extract log-likelihood from LaMaModel object
#' @param object A model fitted using RTMB and obtained via \code{lama(obj)} of class "LaMaModel"
#' @param ... Additional arguments (not used)
#' @return An object of class "logLik"
#' @export
logLik.LaMaModel <- function(object, ...) {
  ll <- object$llk  # your stored log-likelihood
  df <- object$df # number of free parameters
  nobs <- object$nobs  # number of observations
  
  val <- as.numeric(ll)
  attr(val, "df") <- df
  attr(val, "nobs") <- nobs
  class(val) <- "logLik"
  val
}