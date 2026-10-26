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
        param[ind, state] # observed rows, like x[ind]
      } else if(is.vector(param)){
        if(length(param) == 1){
          param
        } else if(length(param) == nObs){
          param[ind]
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

# pseudo-residuals from the predictive CDF at the observations (upper) and, for discrete data,
# at the observations minus one (lower); NA where not observed
finish_res <- function(upper, lower = NULL, observed = !is.na(upper), normal = TRUE, randomise = TRUE){
  upper[!observed] <- NA
  if(!is.null(lower)) lower[!observed] <- NA
  q <- if(normal) qnorm else identity

  if(!is.null(lower)){
    message("Calculating discrete pseudo-residuals")
    if(!randomise){
      out <- list(lower = q(lower), upper = q(upper), mean = q((lower + upper) / 2))
      class(out) <- "LaMaResiduals"
      return(out)
    }
    message("Randomised between lower and upper")
    u <- rep(NA_real_, length(upper))
    u[observed] <- runif(sum(observed), lower[observed], upper[observed])
  } else{
    u <- upper
  }

  u <- q(u)
  u[is.infinite(u)] <- NA
  class(u) <- "LaMaResiduals"
  u
}


#' Calculate pseudo-residuals
#' 
#' @description
#' For HMMs, pseudo-residuals are used to assess the overall goodness-of-fit of the model. 
#' These are based on the cumulative distribution function (CDF)
#' \deqn{F_{X_t}(x_t) = F(x_t \mid x_1, \dots, x_{t-1})}
#' and can be used to quantify whether an observation is extreme relative to its model-implied distribution.
#' 
#' This function calculates such residuals via the probability integral transform, based on the local state probabilities obtained by \code{\link{stateprobs}} and the respective parametric family.
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
#' \strong{Model object only:} \code{pseudo_res(mod)} needs no observations, family or parameters.
#' It evaluates the likelihood function stored in the model again, with each density
#' (e.g. \code{dgamma2()}) replaced by its CDF (\code{pgamma2()}) with the same arguments, so
#' that the likelihood function itself builds the CDFs in place of \code{allprobs}. This works for
#' loops over states or time, \code{cbind()}, covariate-dependent parameters, densities on the log scale,
#' mixtures within a state and several data streams (multiplied, or summed on the log scale),
#' giving one set of residuals per stream, named after the data variable.
#' Before, it checks that the densities enter \code{allprobs} as a product or mixture, so that
#' swapping in CDFs is valid. If not, e.g. for densities transformed by hand (truncation, Jacobians,
#' zero inflation via indicators) or a latent variable integrated out inside \code{allprobs}, it stops, and
#' \code{pseudo_res(obs, dist, par, mod = mod)} needs to be used instead.
#' Densities in your own helper functions called by the likelihood function are found as well,
#' and \code{\link{qreml}} models work like models from \code{\link{lama}}.
#' The likelihood function needs to call \code{\link{forward}} once (use \code{trackID} for several tracks),
#' and it is evaluated with the data as they are when \code{pseudo_res()} is called.
#'
#' @seealso \code{\link{plot.LaMaResiduals}} for plotting pseudo-residuals.
#' 
#' @param obs vector of observations (of length \code{nObs}), or a model object from \code{\link{lama}}
#' (or the \code{RTMB} object), in which case nothing else is needed (see \sQuote{Details})
#' @param dist character string specifying which parametric CDF to use (e.g., \code{"norm"} for normal or \code{"pois"} for Poisson) or CDF function to evaluate directly.
#' If a discrete CDF function is passed, the \code{discrete} argument needs to be set to \code{TRUE} because this cannot determined automatically.
#' @param par named parameter list for the parametric CDF
#' 
#' Names need to correspond to the parameter names in the specified distribution (e.g. \code{list(mean = c(1,2), sd = c(1,1))} for a normal distribution and 2 states).
#' This argument is as flexible as the parametric distribution allows. For example you can have a matrix of parameters with one row for each observation and one column for each state.
#' @param stateprobs matrix of local state probabilities for each observation (of dimension c(nObs, nStates)) as computed by \code{\link{stateprobs}}
#' @param mod optional model object containing required quantities
#' 
#' When using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} directly to this function and avoid calculating local state proabilities manually.
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
#' For \code{pseudo_res(mod)} with several data streams, a named list of these, one per stream, also of class \code{LaMaResiduals}.
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
#' step = trex$step[1:500]
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
#'            logSigma = log(c(0.3, 1.5)))
#'            
#' obj = MakeADFun(nll, par, silent = TRUE)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
#' 
#' mod = lama(obj)
#' 
#' pres = pseudo_res(step,      # observation sequence
#'                   "gamma2",  # parametric family that was used
#'                   list(mean = mod$mu, sd = mod$sigma), # parameters for that family
#'                   mod = mod) # model object
#'
#' # or, from the model object alone
#' pres = pseudo_res(mod)
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
  
  # pseudo_res(mod): everything from the model object
  if(missing(obs) || inherits(obs, c("LaMaModel", "qremlModel")) || (is.list(obs) && !is.null(obs$env))) {
    if(!missing(obs)) mod <- obs
    return(pseudo_res_mod(mod, normal, randomise))
  }
  
  # handle model input: forecast probabilities, checks in stateprobs()
  if(!is.null(mod)){
    stateprobs = stateprobs(mod = mod, forecast = TRUE)
  }
  
  nObs <- length(obs)
  
  if(nrow(stateprobs) != nObs){
    stop("Rows of 'stateprobs' must match length of 'obs'.")
  }
  
  # determine discrete automatically if possible
  if(is.null(discrete)){
    if(is.character(dist)){
      discrete <- dist %in% .discrete_fams
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
  
  observed <- !is.na(obs)
  if(discrete){
    upper <- .eval_cdf_states(obs, par, stateprobs, cdf_func)
    lower <- .eval_cdf_states(obs - 1, par, stateprobs, cdf_func)
    finish_res(upper, lower, observed, normal, randomise)
  } else{
    finish_res(.eval_cdf_states(obs, par, stateprobs, cdf_func), observed = observed, normal = normal)
  }
}


#' Plot pseudo-residuals
#' 
#' @description
#' Plot pseudo-residuals computed by \code{\link{pseudo_res}}.
#' The QQ plot shows a pointwise 95\% band (in gray) for the sample quantiles under a correct model.
#' 
#' @param x pseudo-residuals as returned by \code{\link{pseudo_res}}.
#' For several data streams (from \code{pseudo_res(mod)}), one row of plots per stream, titled by the stream name.
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
#' @importFrom graphics par lines hist abline polygon points
#' @importFrom stats acf na.pass qqnorm ppoints
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
#'            logSigma = log(c(0.3, 1.5)),
#'            logKappa = log(c(0.1, 1.7)))
#'            
#' obj = MakeADFun(nll, par, silent = TRUE)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
#' 
#' mod = lama(obj)
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
#'
#' # or from the model object alone, one row per stream
#' plot(pseudo_res(mod))
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

  # several streams: one row per stream, titled by its name
  if (is.list(x) && all(sapply(x, inherits, "LaMaResiduals"))) {
    old <- par(mfrow = c(length(x), 3))
    on.exit(par(old))
    for (s in names(x)) {
      plot(x[[s]], col = col, lwd = lwd, main = c("", s, ""), breaks = breaks, axis.lab = axis.lab)
    }
    return(invisible(NULL))
  }

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
  
  # QQ Plot with pointwise 95% band, i-th order statistic of n uniforms is Beta(i, n-i+1)
  n <- length(res_clean)
  q <- qnorm(ppoints(n))
  lo <- qnorm(qbeta(0.025, 1:n, n:1))
  hi <- qnorm(qbeta(0.975, 1:n, n:1))
  qqnorm(res_clean, main = main[1],
         xlab = axis.lab$qq[1],
         ylab = axis.lab$qq[2],
         ylim = range(res_clean, lo, hi),
         bty = "n", type = "n")
  polygon(c(q, rev(q)), c(lo, rev(hi)), col = "gray90", border = NA)
  points(qqnorm(res_clean, plot.it = FALSE), pch = 16, col = "#00000070")
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