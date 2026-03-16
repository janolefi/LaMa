# inner helper to eval CDF for both continuous and discrete pseudo residuals
.eval_cdf_states <- function(x, par, stateprobs, cdf_func){
  
  nObs <- length(x)
  N <- ncol(stateprobs)
  ind <- which(!is.na(x))
  
  cdf_values <- matrix(0, nObs, N)
  
  for(state in 1:N){
    current_par <- lapply(par, function(param){
      if(is.matrix(param)){
        if(nrow(param) != nObs | ncol(param) != N){
          stop("Parameter matrix dimensions must match observations and states.")
        }
        param[, state]
      } else if(is.vector(param)){
        if(length(param) == 1 | length(param) == nObs){
          param
        } else{
          if(length(param) != N){
            stop("Parameter vector must have length equal to number of states.")
          }
          param[state]
        }
      } else{
        stop("Invalid parameter specification.")
      }
    })
    
    # safe CDF evaluation
    cdf_values[ind, state] <- tryCatch(
      do.call(cdf_func, c(list(x[ind]), current_par)),
      error = function(e)
        do.call(cdf_func, c(list(pmax(x[ind], 0)), current_par))
    )
  }
  
  rowSums(cdf_values * stateprobs)
}


#' Calculate pseudo-residuals
#' 
#' @description
#' For HMMs, pseudo-residuals are used to assess the overall goodness-of-fit of the model. 
#' These are based on the cumulative distribution function (CDF)
#' \deqn{F_{X_t}(x_t) = F(x_t \mid x_1, \dots, x_{t-1})}
#' and can be used to quantify whether an observation is extreme relative to its model-implied distribution.
#' 
#' This function calculates such residuals via the probability integral transform, based on the local state probabilities obtained by \code{\link{stateprobs}} or \code{\link{stateprobs_g}} and the respective parametric family.
#'
#' @details
#' Pseudo-residuals are based on the probability integral transform. 
#' If the cumulative distribution function (CDF) \eqn{F} of a random variable 
#' \eqn{X} is known, then \eqn{F(X)} is uniformly distributed on \eqn{[0,1]}. 
#' Applying the standard normal quantile function \eqn{\Phi^{-1}} then gives 
#' a residual that is approximately standard normal under the model.
#' 
#' For discrete-valued observations, the CDF is a step function, so the residual 
#' is not uniquely defined. One can compute a "lower" and "upper" pseudo-residual 
#' corresponding to \eqn{F(X-1)} and \eqn{F(X)}, respectively. 
#' If \code{randomise = TRUE}, a value is drawn uniformly 
#' between the lower and upper bounds. In this case, the resulting pseudo-residuals 
#' are again approximately standard normal after applying the standard normal quantile function.
#'
#' @seealso \code{\link{plot.LaMaResiduals}} for plotting pseudo-residuals.
#' 
#' @param obs vector of continuous-valued observations (of length \code{nObs})
#' @param dist character string specifying which parametric CDF to use (e.g., \code{"norm"} for normal or \code{"pois"} for Poisson) or CDF function to evaluate directly.
#' If a discrete CDF function is passed, the \code{discrete} argument needs to be set to \code{TRUE} because this cannot determined automatically.
#' @param par named parameter list for the parametric CDF
#' 
#' Names need to correspond to the parameter names in the specified distribution (e.g. \code{list(mean = c(1,2), sd = c(1,1))} for a normal distribution and 2 states).
#' This argument is as flexible as the parametric distribution allows. For example you can have a matrix of parameters with one row for each observation and one column for each state.
#' @param stateprobs matrix of local state probabilities for each observation (of dimension c(nObs, nStates)) as computed by \code{\link{stateprobs}}, \code{\link{stateprobs_g}} or \code{\link{stateprobs_p}}
#' @param mod optional model object containing required quantities
#' 
#' When using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward}}, \code{\link{forward_g}} or \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{report()} or from \code{\link{qreml}} directly to this function and avoid calculating local state proabilities manually.
#' In this case, a call should look like \code{pseudo_res(obs, "norm", par, mod = mod)}.
#' @param normal logical, if \code{TRUE}, returns Gaussian pseudo residuals
#'
#' These will be approximately standard normally distributed if the model is correct.
#' @param discrete logical, if \code{TRUE}, computes discrete pseudo residuals (which slightly differ in their definition)
#'
#' By default, will be determined using \code{dist} argument, but only works for standard discrete distributions.
#' When used with a special discrete distribution, set to \code{TRUE} manually.
#' @param randomise for discrete pseudo residuals only. Logical, if \code{TRUE}, return randomised pseudo residuals. Recommended for discrete observations.
#'
#' @return vector of pseudo residuals of class \code{LaMaResiduals} or list containing \code{lower}, \code{upper}, and \code{mean} if discrete residuals are not randomised.
#' 
#' @importFrom stats runif
#' @export
#'
#' @examples
#' ## continuous-valued observations
#' obs = rnorm(100)
#' stateprobs = matrix(0.5, nrow = 100, ncol = 2)
#' par = list(mean = c(1,2), sd = c(1,1))
#' pres = pseudo_res(obs, "norm", par, stateprobs)
#'
#' ## discrete-valued observations
#' obs = rpois(100, lambda = 1)
#' par = list(lambda = c(1,2))
#' pres = pseudo_res(obs, "pois", par, stateprobs)
#' 
#' ## custom CDF function
#' obs = rnbinom(100, size = 1, prob = 0.5)
#' par = list(size = c(0.5, 2), prob = c(0.4, 0.6))
#' pres = pseudo_res(obs, pnbinom, par, stateprobs, 
#'                   discrete = TRUE)
#' # if discrete CDF function is passed, 'discrete' needs to be set to TRUE
#' 
#' ## no CDF function available, only density (artificial example)
#' obs = rnorm(100)
#' par = list(mean = c(1,2), sd = c(1,1))
#' # construct CDF using numerical integration
#' cdf = function(x, mean, sd) integrate(dnorm, -Inf, x, mean = mean, sd = sd)$value
#' pres = pseudo_res(obs, cdf, par, stateprobs)
#' 
#' 
#' ### Full model fit example ###
#' step = trex$step[1:200]
#' 
#' nll = function(par){
#'   getAll(par)
#'   Gamma = tpm(logitGamma)
#'   delta = stationary(Gamma)
#'   mu = exp(logMu); REPORT(mu)
#'   sigma = exp(logSigma); REPORT(sigma)
#'   allprobs = matrix(1, length(step), 2)
#'   ind = which(!is.na(step))
#'   for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' par = list(logitGamma = c(-2,-2), 
#'            logMu = log(c(0.3, 2.5)), 
#'            logSigma = log(c(0.3, 0.5)))
#'            
#' obj = MakeADFun(nll, par, silent = TRUE)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
#' 
#' mod = obj$report()
#' 
#' pres = pseudo_res(step,      # observation sequence
#'                   "gamma2",  # parametric family that was used
#'                   list(mean = mod$mu, sd = mod$sigma), # parameters for that family
#'                   mod = mod) # model object
#'
#' plot(pres)
pseudo_res <- function(obs,
                       dist,
                       par,
                       stateprobs = NULL,
                       mod = NULL,
                       normal = TRUE,
                       discrete = NULL,
                       randomise = TRUE) {
  
  # handle model input
  if(!is.null(mod)){
    # checks
    if(is.null(mod$type)){
      stop("'mod' contains no type.")
    }
    if(!(mod$type) %in% c("homogeneous","inhomogeneous","periodic","continuous_time")){
      stop("'mod' contains invalid type.")
    }
    if(is.null(mod$delta)) stop("'mod' contains no initial distribution.")
    if(is.null(mod$Gamma)) stop("'mod' contains no transition matrix.")
    if(is.null(mod$allprobs)) stop("'mod' contains no state-dependent probabilities.")
    if(mod$type == "periodic" && is.null(mod$tod)){
      stop("'mod' contains no cyclic indexing variable.")
    }
    if(mod$type == "homogeneous"){
      stateprobs = stateprobs(mod = mod, forecast = TRUE)
    }
    if(mod$type %in% c("inhomogeneous","continuous_time")){
      stateprobs = stateprobs_g(mod = mod, forecast = TRUE)
    }
    if(mod$type == "periodic"){
      stateprobs = stateprobs_p(mod = mod, forecast = TRUE)
    }
  }
  
  nObs <- length(obs)
  
  if(nrow(stateprobs) != nObs){
    stop("Rows of 'stateprobs' must match length of 'obs'.")
  }
  
  # determine discrete automatically if possible
  if(is.null(discrete)){
    if(is.character(dist)){
      discrete <- dist %in% c("pois","binom","geom","nbinom","nbinom2",
                              "betabinom","genpois","zibinom","zibinom2",
                              "zipois","ztbinom","ztnbinom","ztnbinom2", "ztpois")
    } else if(is.function(dist)){
      discrete <- FALSE
      message("Assuming 'dist' evaluates a continuous CDF. Set 'discrete = TRUE' if not.")
    }
  }
  
  # obtain CDF
  if(is.character(dist)){
    cdf_func <- get(paste0("p", dist), mode = "function")
  } else if(is.function(dist)){
    cdf_func <- Vectorize(dist)
  } else{
    stop("'dist' must be a character string or function.")
  }
  
  if(discrete){
    message("Calculating discrete pseudo-residuals")
    
    u_upper <- .eval_cdf_states(obs, par, stateprobs, cdf_func)
    u_lower <- .eval_cdf_states(obs - 1, par, stateprobs, cdf_func)
    
    if(randomise){
      message("Randomised between lower and upper")
      u <- rep(NA_real_, nObs)
      ind <- which(!is.na(obs))
      u[ind] <- runif(length(ind), u_lower[ind], u_upper[ind])
    } else{
      if(normal){
        out <- list(
          lower = qnorm(u_lower),
          upper = qnorm(u_upper),
          mean  = qnorm((u_lower + u_upper)/2)
        )
        class(out) <- "LaMaResiduals"
        return(out)
      } else{
        out <- list(
          lower = u_lower,
          upper = u_upper,
          mean  = (u_lower + u_upper)/2
        )
        class(out) <- "LaMaResiduals"
        return(out)
      }
    }
  } else{
    u <- .eval_cdf_states(obs, par, stateprobs, cdf_func)
  }
  
  if(normal){
    u <- qnorm(u)
  }
  
  u[is.infinite(u)] <- NA
  
  class(u) <- "LaMaResiduals"
  
  return(u)
}


#' Plot pseudo-residuals
#' 
#' @description
#' Plot pseudo-residuals computed by \code{\link{pseudo_res}}.
#' 
#' @param x pseudo-residuals as returned by \code{\link{pseudo_res}}
#' @param col character, color for the QQ-line (and density curve if \code{histogram = TRUE})
#' @param lwd numeric, line width for the QQ-line (and density curve if \code{histogram = TRUE})
#' @param main optional character vector of main titles for the plots of length 2 (or 3 if \code{histogram = TRUE})
#' @param breaks \code{breaks} argument passed to hist
#' @param axis.lab labels used for the x and y axis of each plot (named list)
#' @param ... currently ignored. For method consistency
#' 
#' @returns NULL, plots the pseudo-residuals in a 3-panel layout
#' @export
#'
#' @importFrom graphics par lines hist abline
#' @importFrom stats acf na.pass qqnorm
#' @importFrom RTMB dnorm
#' @importFrom utils modifyList
#'
#' @examples
#' ## pseudo-residuals for the trex data
#' step = trex$step[1:1000]
#' angle = trex$angle[1:1000]
#' 
#' nll = function(par){
#'   getAll(par)
#'   Gamma = tpm(logitGamma)
#'   delta = stationary(Gamma)
#'   mu = exp(logMu); REPORT(mu)
#'   sigma = exp(logSigma); REPORT(sigma)
#'   kappa = exp(logKappa); REPORT(kappa)
#'   allprobs = matrix(1, length(step), 2)
#'   ind = which(!is.na(step) & !is.na(angle))
#'   for(j in 1:2) {
#'     allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j]) * 
#'          dvm(angle[ind], 0, kappa[j])
#'   }
#'   -forward(delta, Gamma, allprobs)
#' }
#' 
#' par = list(logitGamma = c(-2,-2), 
#'            logMu = log(c(0.3, 2.5)), 
#'            logSigma = log(c(0.3, 0.5)),
#'            logKappa = log(c(0.2, 1)))
#'            
#' obj = MakeADFun(nll, par)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
#' 
#' mod = report(obj)
#' 
#' pres_step = pseudo_res(step,      # observations
#'                        "gamma2",  # family that is used
#'                        list(mean = mod$mu, sd = mod$sigma), # the family's parameters
#'                        mod = mod) # model object
#' pres_angle = pseudo_res(angle,
#'                         "vm",
#'                         list(mu = 0, kappa = mod$kappa),
#'                         mod = mod)
#'                         
#' # separate plots 
#' plot(pres_step)
#' plot(pres_angle)
#' 
#' # together
#' par(mfrow = c(2,3))
#' plot(pres_step, main = c("", "Step Length", ""), 
#'      axis.lab = list(hist = c("Step residuals", "Density")))
#' plot(pres_angle, main = c("", "Turning Angle", ""),
#'      axis.lab = list(hist = c("Angle residuals", "Density")))
plot.LaMaResiduals <- function(
    x, 
    col = "darkblue", 
    lwd = 1.5,
    main = NULL,
    breaks = "Sturges",
    axis.lab = list(qq = c("Theoretical quantiles", "Sample quantiles"), hist = c("Pseudo-residuals", "Density"), acf = c("Lag", "ACF")),
    ...
    ) {
  
  ## default axis labels
  default_axis <- list(
    qq   = c("Theoretical quantiles", "Sample quantiles"),
    hist = c("Pseudo-residuals", "Density"),
    acf  = c("Lag", "ACF")
  )
  
  axis.lab <- modifyList(default_axis, axis.lab)
  
  # Extract mean if residuals is a list
  res <- if (is.list(x)) x$mean else x
  res_clean <- na.omit(res)
  
  if (is.null(main)) main <- rep("", 3)
  if (length(main) < 3) {
    main <- c(main, rep("", 3 - length(main)))
  }
 
  mf <- par("mfrow")
  if (all(mf == c(1,1))) {
    old <- par(mfrow = c(1,3))
    on.exit(par(old))
  } else if (mf[2] != 3) {
    old <- par(mfrow = c(1,3))
    on.exit(par(old))
  }
  
  # QQ Plot
  qqnorm(res_clean, main = main[1], 
         xlab = axis.lab$qq[1], 
         ylab = axis.lab$qq[2],
         bty = "n", pch = 16, col = "#00000070")
  abline(a = 0, b = 1, col = col, lwd = lwd)
  
  # Histogram
  r <- range(res_clean)
  dr <- diff(r)
  xgrid <- seq(r[1] - dr / 4, r[2] + dr / 4, length.out = 200)
  dens <- dnorm(xgrid)
  ylim <- c(0, max(dens) * 1.1)
  
  hist(res_clean, main = main[2], border = "white", ylim = ylim,
       prob = TRUE, breaks = breaks,
       xlab = axis.lab$hist[1],
       ylab = axis.lab$hist[2])
  
  lines(xgrid, dens, col = col, lwd = lwd)
  
  # ACF
  acf(res, na.action = na.pass, main = main[3],
      xlab = axis.lab$acf[1],
      ylab = axis.lab$acf[2],
      bty = "n")
}