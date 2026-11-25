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
#' states <- viterbi(mod)   # global decoding
#' probs <- stateprobs(mod) # local decoding
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
  mod$nobs <- count_nobs(mod)
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

# Number of observations of a model (or the reported quantities): rows of allprobs that are not all 1,
# as these are missing and contribute nothing to the likelihood, without the dropped ones of HSMMs
count_nobs <- function(mod){
  if(is.null(mod$allprobs)) return(NULL)
  sum(rowSums(mod$allprobs != 1) > 0 & used_obs(mod))
}

# Which observations are used? HSMMs drop the observations before first_obs in each track
used_obs <- function(mod){
  nObs <- nrow(mod$allprobs)
  first_obs <- mod[["first_obs"]]
  if(is.null(first_obs) || first_obs == 1) return(rep(TRUE, nObs))
  tracks <- prepare_tracks(mod$trackID, nObs)
  index <- seq_len(nObs) - rep(tracks$start, tracks$end - tracks$start + 1) + 1 # 1, 2, ... in each track
  index >= first_obs
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

#' Estimated parameters of a \code{LaMaModel} object
#'
#' Returns the estimated fixed parameters of a model from \code{\link{lama}} or \code{\link{qreml}},
#' on the working scale and in the order of \code{obj$par}, i.e. without parameters fixed via \code{map}
#' and without random effects.
#' Parameters with several entries are named by their position among the estimated entries,
#' e.g. \code{logmu.1} and \code{logmu.2}, matrices in column-major order as in \code{obj$par}.
#'
#' @param object model object of class \code{"LaMaModel"}
#' @param ... additional arguments (not used)
#'
#' @returns named numeric vector of the estimated parameters
#' @seealso \code{\link{vcov.LaMaModel}}
#' @importFrom stats coef
#' @export
coef.LaMaModel <- function(object, ...) {
  p <- est_par(object)
  random <- get_obj(object)$env$random
  if(length(random)) p <- p[-random]
  # unique names, e.g. logmu.1, logmu.2, so confint() and indexing by name work
  nm <- names(p)
  rep <- nm %in% nm[duplicated(nm)]
  nm[rep] <- paste0(nm[rep], ".", stats::ave(seq_along(nm), nm, FUN = seq_along)[rep])
  names(p) <- nm
  p
}

#' Variance-covariance matrix of a \code{LaMaModel} object
#'
#' Computes the approximate covariance matrix of the estimated parameters \code{coef(object)},
#' either as the inverse Hessian of the negative log-likelihood at the estimate (\code{type = "hessian"}),
#' or as a sandwich that remains valid if the model is misspecified: with the outer product of the scores
#' (\code{type = "sandwich"}), heteroskedasticity and autocorrelation consistent (\code{type = "HAC"}) or
#' cluster-robust over tracks (\code{type = "cluster"}).
#'
#' @details
#' \strong{Hessian.} The Hessian \eqn{H} is
#' \itemize{
#'   \item the exact Hessian \code{obj$he()} for models without random effects,
#'   \item for models with random effects, the Hessian of the Laplace-approximated marginal likelihood
#'   w.r.t. the fixed parameters, by finite differences of \code{obj$gr()} (as in \code{\link[RTMB]{sdreport}}),
#'   \item for models from \code{\link{qreml}}, the Hessian of the penalised negative log-likelihood,
#'   including the penalised coefficients. The covariance is then conditional on the estimated smoothing parameters.
#' }
#' If the Hessian is not positive definite, which indicates that the model has not converged,
#' a warning is issued.
#'
#' \strong{Sandwich.} With the per-observation score contributions \eqn{u_t} from \code{\link{scores}},
#' the covariance is \eqn{H^{-1} M H^{-1}} with the outer product of the scores as meat, \eqn{M = \sum_t u_t u_t^\top}
#' (the outer product of gradients, OPG). It is robust to a misspecified predictive distribution as long as the \eqn{u_t}
#' are serially uncorrelated, which they are under correct specification (martingale differences), but not to
#' dependence the model misses. It equals the HAC sandwich with the Bartlett kernel and \code{bw} at most 1. Not to be confused with the inverse of the
#' outer product alone, \eqn{M^{-1}} (the BHHH estimator), which, like \eqn{H^{-1}}, is only valid under correct specification.
#'
#' \strong{HAC sandwich.} With the per-observation score contributions \eqn{u_t} from \code{\link{scores}},
#' the covariance is \eqn{H^{-1} M H^{-1}}, with the kernel estimate of the variance of the score sum
#' \deqn{M = \Gamma_0 + \sum_{l \ge 1} k(l / b) (\Gamma_l + \Gamma_l^\top), \qquad
#' \Gamma_l = \sum_t u_t u_{t-l}^\top,}
#' where the sum in \eqn{\Gamma_l} only runs over pairs of observations from the same track (as reported by \code{\link{forward}}),
#' as different tracks are independent. Missing observations contribute zero scores, so lags remain in units of time.
#' The kernel \eqn{k} down-weights the autocovariances at longer lags, which are estimated with much noise, and the
#' bandwidth \eqn{b} (\code{bw}) scales the lag axis. Two kernels are available:
#' \itemize{
#'   \item \code{kernel = "QS"} (default), the Quadratic Spectral kernel
#'   \deqn{k(x) = \frac{25}{12 \pi^2 x^2} \Bigl(\frac{\sin(6 \pi x / 5)}{6 \pi x / 5} - \cos(6 \pi x / 5)\Bigr).}
#'   Its weights are positive up to lag \eqn{1.2 b} and then oscillate with small amplitude, so all lags contribute.
#'   It is the most efficient kernel among those that guarantee a positive semi-definite \eqn{M} (Andrews, 1991):
#'   it keeps the short lags at nearly full weight, which reduces the bias.
#'   \item \code{kernel = "Bartlett"}, \eqn{k(x) = \max(1 - |x|, 0)}, the Newey-West estimator: with \eqn{b = L + 1},
#'   lags up to \eqn{L} enter with weights \eqn{1 - l / (L + 1)}.
#' }
#' The default bandwidth is the data-driven choice of Andrews (1991), which minimises the asymptotic mean squared error
#' of \eqn{M} for an AR(1) approximation of the scores: with \eqn{\rho_a} and \eqn{\sigma_a^2} the coefficient and innovation
#' variance of an AR(1) with intercept fitted to each score column \eqn{a} (by least squares on consecutive observations
#' within tracks),
#' \deqn{b = 1.3221 \, (\alpha_2 \bar T)^{1/5} \; \text{(QS)}, \qquad b = 1.1447 \, (\alpha_1 \bar T)^{1/3} \; \text{(Bartlett)},}
#' with \eqn{\bar T} the mean track length, \eqn{\alpha_2 = \sum_a 4 \rho_a^2 \sigma_a^4 (1 - \rho_a)^{-8} / D},
#' \eqn{\alpha_1 = \sum_a 4 \rho_a^2 \sigma_a^4 (1 - \rho_a)^{-6} (1 + \rho_a)^{-2} / D} and
#' \eqn{D = \sum_a \sigma_a^4 (1 - \rho_a)^{-4}}. For a single track, this is \code{sandwich::bwAndrews()} without
#' prewhitening and with equal weights. All score columns are weighted equally, so the bandwidth is driven by the most persistent and most variable ones.
#' Stronger serial dependence of the scores gives a larger bandwidth.
#'
#' All lags are computed exactly, with the fast Fourier transform within each track, which takes
#' \eqn{O(p T \log T + T p^2)} operations for \eqn{T} observations and \eqn{p} parameters instead of \eqn{O(T^2 p^2)} for
#' summing all lagged products directly.
#'
#' If the model is correctly specified, the \eqn{u_t} are martingale differences, so \eqn{\Gamma_l \approx 0} for \eqn{l \ge 1}
#' and \eqn{M \approx H}: both types then agree approximately, and marked differences point to misspecification.
#'
#' \strong{Cluster-robust sandwich.} With the score sums \eqn{s_k = \sum_{t \in k} u_t} of the clusters \eqn{k = 1, \dots, G}
#' (by default the tracks), the meat is
#' \deqn{M = c \sum_{k=1}^G s_k s_k^\top,}
#' with \eqn{c = 1} (\code{cr = "CR0"}) or \eqn{c = G / (G - 1)} (\code{cr = "CR1"}, default).
#' It makes no assumption about the dependence within clusters, so no bandwidth is needed, but requires
#' independent clusters and is only reliable for many of them. As the \eqn{s_k} sum to the gradient,
#' which is zero at the estimate, \eqn{M} has rank at most \eqn{G - 1}, so the covariance is singular if \eqn{G - 1} is less than
#' the number of parameters. For parameters informed by a single cluster only, such as track-specific initial distributions,
#' it is not meaningful. Clusters may combine several tracks (e.g. several tracks of one individual), but a track must not
#' be split between clusters.
#'
#' \code{cr = "CR3"} gives the jackknife estimator (MacKinnon, Nielsen and Webb, 2023)
#' \deqn{V = \frac{G - 1}{G} \sum_{k=1}^G d_k d_k^\top, \qquad d_k = \hat\theta - \hat\theta^{(k)},}
#' where \eqn{\hat\theta^{(k)}} is the estimate without cluster \eqn{k}. It is found by minimising the weighted negative
#' log-likelihood with weight zero on the observations of cluster \eqn{k}, started at \eqn{\hat\theta}, so the model is refitted
#' once per cluster, and a message says so. Compared to CR1, it inflates the contributions of clusters with high leverage
#' (e.g. long tracks), whose score sums are shrunk the most by fitting. Refits that do not converge give a warning.
#' If the Hessian at the estimate of a refit is not positive definite, as for parameters not or only weakly identified without
#' that cluster, an error is issued. That a refit has moved to another mode, for example with a different order of the states,
#' is not checked; standard errors much larger than those of the one-step version can be a sign of this.
#'
#' \code{cr = "CR3onestep"} approximates the change \eqn{d_k} by one Newton step from the estimate,
#' \deqn{d_k = (H - H_k)^{-1} s_k,}
#' where \eqn{H - H_k} is the Hessian without cluster \eqn{k}, from the weighted likelihood with weight zero on its observations.
#' It costs one Hessian per cluster and no refits. For linear least squares, the step is exact; otherwise, the approximation is
#' poorest if there are few clusters or if one of them carries much of the information. With equal leverage
#' \eqn{H_k = H / G}, it coincides with CR1. The one-step deletion is a pure Newton step, so it requires \eqn{H - H_k} to be
#' positive definite: if leaving out a cluster makes it singular or indefinite, an error is issued.
#'
#' All sandwiches require the likelihood function to call \code{\link{forward}} once (with \code{trackID} for several tracks)
#' and are not available for models with random effects or from \code{\link{qreml}}.
#'
#' @param object model object of class \code{"LaMaModel"}
#' @param type \code{"hessian"} (default) for the inverse Hessian, \code{"sandwich"} for the sandwich with the outer
#' product of the scores, \code{"HAC"} for the HAC sandwich, or \code{"cluster"} for the cluster-robust sandwich
#' @param kernel kernel of the HAC sandwich: \code{"QS"} (Quadratic Spectral, default) or \code{"Bartlett"} (Newey-West),
#' see Details
#' @param bw bandwidth of the HAC sandwich, a non-negative number. Defaults to the data-driven bandwidth of Andrews (1991)
#' for the chosen kernel, see Details. For the Bartlett kernel, \code{bw = L + 1} includes the lags up to \eqn{L},
#' and \code{bw} at most 1 gives the same as \code{type = "sandwich"}.
#' @param cluster cluster identifier for \code{type = "cluster"}, one value per observation (of any type, in the order of the
#' observations in \code{\link{forward}}). Defaults to the tracks, i.e. the \code{trackID} reported by \code{\link{forward}}.
#' @param cr small-sample adjustment for \code{type = "cluster"}: \code{"CR1"} (default) multiplies the meat by \eqn{G / (G - 1)},
#' \code{"CR0"} does not adjust, \code{"CR3"} is the jackknife, which refits the model without each cluster, and
#' \code{"CR3onestep"} its approximation by one Newton step per cluster, see Details.
#' @param ... additional arguments (not used)
#'
#' @returns covariance matrix with rows and columns named like \code{coef(object)}
#' @references Andrews, D. W. K. (1991). Heteroskedasticity and autocorrelation consistent covariance matrix estimation.
#' \emph{Econometrica}, 59(3), 817--858.
#'
#' Newey, W. K. and West, K. D. (1987). A simple, positive semi-definite, heteroskedasticity
#' and autocorrelation consistent covariance matrix. \emph{Econometrica}, 55(3), 703--708.
#'
#' MacKinnon, J. G., Nielsen, M. O. and Webb, M. D. (2023). Cluster-robust inference: A guide to empirical practice.
#' \emph{Journal of Econometrics}, 232(2), 272--299.
#'
#' Cameron, A. C. and Miller, D. L. (2015). A practitioner's guide to cluster-robust inference.
#' \emph{Journal of Human Resources}, 50(2), 317--372.
#' @seealso \code{\link{coef.LaMaModel}}, \code{\link{scores}}
#' @importFrom stats vcov optimHess
#' @export
vcov.LaMaModel <- function(object, type = c("hessian", "sandwich", "HAC", "cluster"), kernel = c("QS", "Bartlett"),
                           bw = NULL, cluster = NULL, cr = c("CR1", "CR0", "CR3", "CR3onestep"), ...) {
  type <- match.arg(type)
  kernel <- match.arg(kernel)
  cr <- match.arg(cr)
  obj <- get_obj(object)
  if(type != "hessian") {
    if(inherits(object, "qremlModel")) stop("type = \"", type, "\" is not available for models from qreml().")
    if(length(obj$env$random)) stop("type = \"", type, "\" is not available for models with random effects.")
  }
  p <- coef(object)
  H <- if(inherits(object, "qremlModel")) object$Hessian_conditional
       else if(length(obj$env$random)) optimHess(p, obj$fn, obj$gr)
       else obj$he(p)
  H <- (H + t(H)) / 2

  V <- tryCatch(chol2inv(chol(H)), error = function(e) NULL)
  if(is.null(V)) {
    warning("Hessian is not positive definite, the model may not have converged ",
            "(maximum gradient component ", signif(object$mgc, 3), ").", call. = FALSE)
    V <- tryCatch(solve(H), error = function(e) matrix(NA_real_, nrow(H), ncol(H)))
  }

  if(type != "hessian") {
    tp <- weighted_tape(object)
    U <- tape_scores(tp)
    if(type == "cluster" && cr %in% c("CR3", "CR3onestep")) {
      V <- cr3_vcov(tp, U, cluster_ids(U, object$trackID, cluster), onestep = (cr == "CR3onestep"))
    } else {
      M <- switch(type,
                  sandwich = crossprod(U),
                  HAC = hac_meat(U, object$trackID, kernel, bw),
                  cluster = cluster_meat(U, object$trackID, cluster, cr, npar = length(p)))
      V <- V %*% M %*% V
    }
    V <- (V + t(V)) / 2
  }
  dimnames(V) <- list(names(p), names(p))
  V
}

#' Takeuchi information criterion
#'
#' Computes the Takeuchi information criterion (TIC), a version of AIC whose penalty remains valid if the model is
#' misspecified, for one or several fitted models.
#'
#' @details
#' With \eqn{H} the Hessian of the negative log-likelihood and \eqn{J} the variance of the score sum (the meat of the
#' sandwich, see \code{\link{vcov.LaMaModel}}),
#' \deqn{\mathrm{TIC} = -2 \log L(\hat\theta) + 2 \, \mathrm{edf}, \qquad \mathrm{edf} = \mathrm{tr}(J H^{-1}).}
#' The penalty of AIC, the number of parameters \eqn{p}, assumes \eqn{J = H}, which holds if the model is correctly
#' specified. TIC estimates it instead: the effective number of parameters \code{edf} approximates the optimism
#' \eqn{E[(\hat\theta - \theta^*)^\top H (\hat\theta - \theta^*)] \approx \mathrm{tr}(H \, \mathrm{Var}(\hat\theta))}
#' of the maximised log-likelihood, so models whose estimates are more variable than the model itself suggests are
#' penalised more.
#'
#' It is computed from the robust covariance \eqn{V = H^{-1} J H^{-1}} as \eqn{\mathrm{edf} = \mathrm{tr}(H V)}, so
#' every robust \code{type} of \code{\link{vcov.LaMaModel}} can be used, including the jackknives \code{cr = "CR3"}
#' and \code{cr = "CR3onestep"}.
#' The default is the HAC sandwich, as the penalty needs the correct variance of the estimator, which for serially
#' dependent scores is not given by the outer product of the scores (\code{type = "sandwich"}). For several
#' independent tracks, \code{type = "cluster"} is an alternative.
#'
#' Under correct specification, \code{edf} is close to \eqn{p}, and TIC close to AIC. Values well above \eqn{p} indicate
#' misspecification, e.g. serial dependence the model misses. Values below \eqn{p} arise from noise in the
#' estimated \eqn{J}, from a model that overstates the variability of the data, or from parameters that are nearly not
#' identified, e.g. a superfluous state. In the last case, neither AIC nor TIC is justified, so neither should be
#' used to choose the number of states. \code{edf} is itself estimated and, especially for short series, noisy.
#'
#' @param object model object of class \code{"LaMaModel"}
#' @param ... further models of class \code{"LaMaModel"}, compared as in \code{\link[stats]{AIC}}
#' @param type robust covariance used for the penalty: \code{"HAC"} (default), \code{"sandwich"} or \code{"cluster"},
#' see \code{\link{vcov.LaMaModel}}
#' @param kernel,bw,cluster,cr passed to \code{\link{vcov.LaMaModel}}
#'
#' @returns For one model, the TIC, with the effective number of parameters as attribute \code{"edf"}.
#' For several, a data frame with one row per model and columns \code{df} (number of parameters), \code{edf} and \code{TIC}.
#' @references Takeuchi, K. (1976). Distribution of informational statistics and a criterion of model fitting.
#' \emph{Suri-Kagaku (Mathematical Sciences)}, 153, 12--18 (in Japanese).
#'
#' Konishi, S. and Kitagawa, G. (2008). \emph{Information Criteria and Statistical Modeling}. Springer, New York.
#' @seealso \code{\link{vcov.LaMaModel}}, \code{\link[stats]{AIC}}
#' @export
TIC <- function(object, ...) UseMethod("TIC")

#' @rdname TIC
#' @export
TIC.LaMaModel <- function(object, ..., type = c("HAC", "sandwich", "cluster"), kernel = c("QS", "Bartlett"),
                          bw = NULL, cluster = NULL, cr = c("CR1", "CR0", "CR3", "CR3onestep")) {
  type <- match.arg(type)
  kernel <- match.arg(kernel)
  cr <- match.arg(cr)
  models <- c(list(object), list(...))
  if(!all(vapply(models, inherits, TRUE, what = "LaMaModel"))) stop("All models must be of class \"LaMaModel\".")
  out <- data.frame(df = numeric(0), edf = numeric(0), TIC = numeric(0))
  for(m in models) {
    H <- solve(vcov(m)) # Hessian of the negative log-likelihood
    V <- vcov(m, type = type, kernel = kernel, bw = bw, cluster = cluster, cr = cr)
    edf <- sum(H * V)   # tr(H V), as H and V are symmetric
    out <- rbind(out, data.frame(df = length(coef(m)), edf = edf, TIC = -2 * as.numeric(logLik(m)) + 2 * edf))
  }
  if(length(models) == 1) return(structure(out$TIC, edf = out$edf))
  # row names as in AIC(): the model arguments as written in the call
  cl <- match.call(expand.dots = FALSE)
  rownames(out) <- c(deparse(cl$object), vapply(cl$..., deparse, ""))
  out
}

# HAC meat from score contributions U (one row per observation): M = sum_l k(l / bw) Gamma_l over all lags,
# with autocovariances only within tracks. Tracks are contiguous, as forward() requires.
hac_meat <- function(U, trackID = NULL, kernel = "QS", bw = NULL){
  n <- nrow(U)
  if(is.null(trackID)) trackID <- rep(1L, n)
  if(length(trackID) != n) {
    stop("Length of the reported trackID (", length(trackID), ") does not match the number of score contributions (", n,
         "). The likelihood function should call forward() once, with trackID for several tracks.")
  }
  tracks <- split(seq_len(n), factor(trackID, levels = unique(trackID))) # row indices of each track
  if(is.null(bw)) bw <- bw_andrews(U, tracks, kernel)
  if(!is.numeric(bw) || length(bw) != 1 || !is.finite(bw) || bw < 0) stop("bw must be a non-negative number.")

  M <- matrix(0, ncol(U), ncol(U))
  for(i in tracks) {
    lags <- seq_along(i) - 1
    w <- if(bw > 0) hac_kernel(lags / bw, kernel) else as.numeric(lags == 0) # bw = 0: lag 0 only
    M <- M + fft_meat(U[i, , drop = FALSE], w)
  }
  M
}

# kernel weights k(x), x = lag / bandwidth
hac_kernel <- function(x, kernel = "QS"){
  x <- abs(x)
  if(kernel == "Bartlett") return(pmax(1 - x, 0))
  z <- 6 * pi * x / 5
  w <- 25 / (12 * pi^2 * x^2) * (sin(z) / z - cos(z))
  w[x == 0] <- 1
  w
}

# Gamma_0 + sum_{l >= 1} w[l + 1] (Gamma_l + Gamma_l') for one track, over all lags l = 1, ..., n - 1, by FFT:
# with F the Fourier transforms of the score columns and W that of the weights, the weighted sum of all lagged
# products is sum_j W_j F_j F_j^H / (2n). Zero-padding to length 2n avoids wrap-around, so this is exact.
fft_meat <- function(U, w){
  n <- nrow(U)
  if(n == 1) return(w[1] * crossprod(U))
  F <- stats::mvfft(rbind(U, matrix(0, n, ncol(U))))  # Fourier transform of each zero-padded score column
  W <- Re(stats::fft(c(w, 0, rev(w[-1]))))     # weights of lags 0, ..., n - 1, (n), -(n - 1), ..., -1
  Re(crossprod(F, W * Conj(F))) / (2 * n)
}

# Andrews (1991) bandwidth from AR(1) fits to each score column (least squares with intercept on the consecutive
# pairs within tracks, as stats::ar.ols()), all columns weighted equally; T is the mean track length.
# For a single track, as sandwich::bwAndrews(prewhite = 0, weights = 1).
bw_andrews <- function(U, tracks, kernel = "QS"){
  t0 <- unlist(lapply(tracks, function(i) i[-length(i)])) # pairs (t0, t1 = t0 + 1) within a track
  t1 <- unlist(lapply(tracks, function(i) i[-1]))
  num <- den <- 0
  for(a in seq_len(ncol(U))) {
    y <- U[t1, a] - mean(U[t1, a])
    z <- U[t0, a] - mean(U[t0, a])
    if(sum(z^2) == 0) next # constant score column, e.g. a parameter with no information
    rho <- sum(y * z) / sum(z^2)  # least squares slope of u_t on u_{t-1}
    s2 <- mean((y - rho * z)^2)   # innovation variance
    den <- den + s2^2 / (1 - rho)^4
    num <- num + 4 * rho^2 * s2^2 / (if(kernel == "QS") (1 - rho)^8 else (1 - rho)^6 * (1 + rho)^2)
  }
  if(den == 0) return(0)
  nT <- nrow(U) / length(tracks)
  if(kernel == "QS") 1.3221 * (num / den * nT)^(1/5) else 1.1447 * (num / den * nT)^(1/3)
}

# Cluster identifiers, one per observation (default: the tracks), checked against the scores and tracks
cluster_ids <- function(U, trackID = NULL, cluster = NULL){
  n <- nrow(U)
  if(!is.null(trackID) && length(trackID) != n) {
    stop("Length of the reported trackID (", length(trackID), ") does not match the number of score contributions (", n,
         "). The likelihood function should call forward() once, with trackID for several tracks.")
  }
  if(is.null(cluster)) {
    if(is.null(trackID)) stop("type = \"cluster\" needs several tracks (trackID in forward()) or a cluster argument.")
    cluster <- trackID
  } else {
    if(length(cluster) != n) stop("cluster must have one value per observation (", n, "), but has length ", length(cluster), ".")
    if(anyNA(cluster)) stop("cluster must not contain missing values.")
    # tracks are dependent internally, so splitting one between clusters is invalid
    if(!is.null(trackID) && any(tapply(cluster, trackID, function(x) length(unique(x))) > 1)) {
      stop("A track is split between clusters: each track must lie within a single cluster.")
    }
  }
  if(length(unique(cluster)) < 2) stop("type = \"cluster\" needs at least two clusters.")
  cluster
}

# Cluster-robust meat c * sum_k s_k s_k' from score contributions U, with s_k the score sums of
# the clusters (default: the tracks) and c = G / (G - 1) for CR1, c = 1 for CR0.
cluster_meat <- function(U, trackID = NULL, cluster = NULL, cr = "CR1", npar = ncol(U)){
  cluster <- cluster_ids(U, trackID, cluster)
  S <- rowsum(U, cluster, reorder = FALSE) # G x p score sums
  G <- nrow(S)
  if(G - 1 < npar) {
    warning("Only ", G, " clusters for ", npar, " parameters: the cluster-robust covariance has rank at most ",
            G - 1, " and is singular.", call. = FALSE)
  }
  adj <- if(cr == "CR1") G / (G - 1) else 1
  adj * crossprod(S)
}

# CR3 (jackknife): (G - 1) / G * sum_k d_k d_k', with d_k the change in the estimate (up to its sign)
# when leaving out cluster k, exactly by refitting (default) or approximated by one Newton step.
cr3_vcov <- function(tp, U, cluster, onestep = FALSE){
  grp <- match(cluster, unique(cluster)) # 1, ..., G in order of appearance
  S <- rowsum(U, grp, reorder = FALSE)
  G <- nrow(S)
  if(G < tp$p) {
    warning("Only ", G, " clusters for ", tp$p, " parameters: the cluster-robust covariance has rank at most ",
            G, " and is singular.", call. = FALSE)
  }
  if(!onestep) {
    message("CR3: refitting the model without each of the ", G, " clusters.")
    R <- chol_or_null(weighted_hessian(tp)) # scales the refits, one Hessian for all clusters
  }
  D <- matrix(0, G, tp$p)
  converged <- rep(TRUE, G)
  for(k in seq_len(G)) {
    label <- format(unique(cluster)[k])
    if(onestep) {
      D[k, ] <- delete_onestep(tp, grp != k, S[k, ], label)
    } else {
      del <- delete_refit(tp, grp != k, R, label)
      D[k, ] <- del$d
      converged[k] <- del$converged
    }
  }
  if(any(!converged)) {
    warning("CR3: the refits without cluster(s) ", paste(unique(cluster)[!converged], collapse = ", "),
            " did not converge.", call. = FALSE)
  }
  (G - 1) / G * crossprod(D)
}

# One-step deletion: (H - H_k)^{-1} s_k, one Newton step from the estimate on the likelihood without
# cluster k (exact for least squares). H - H_k is the Hessian with weights w, zero on cluster k.
delete_onestep <- function(tp, w, s_k, label){
  R <- chol_or_null(weighted_hessian(tp, w))
  # a pure Newton step is only meaningful for a positive definite Hessian
  if(is.null(R)) {
    stop("The Hessian without cluster ", label, " is singular or not positive definite: ",
         "some parameters are not or only weakly identified without this cluster (e.g. track-specific parameters), ",
         "so CR3onestep is not available.")
  }
  backsolve(R, forwardsolve(t(R), s_k))
}

# Exact deletion: refit with weights w, zero on cluster k, started at the estimate. The weights are
# fixed on the tape, so neither the data nor the likelihood function are touched. nlminb() uses the
# gradient only, so the refit does not depend on the Hessian without cluster k being positive definite.
# The refit is in the parameters phi = R (theta - theta_hat), with R the Cholesky factor of the full Hessian
# H = R'R: the curvature in phi is then close to the identity, from which the quasi-Newton updates
# start, so a few iterations suffice. Without R (H not positive definite) theta is not rescaled.
# Returns the change in the estimate and whether nlminb() converged.
delete_refit <- function(tp, w, R, label){
  w <- as.numeric(w)
  th <- unname(tp$theta)
  if(is.null(R)) R <- diag(tp$p)
  theta <- function(phi) th + backsolve(R, phi)
  fn <- function(phi) tp$F(c(theta(phi), w))
  gr <- function(phi) backsolve(R, tp$G(c(theta(phi), w))[seq_len(tp$p)], transpose = TRUE)
  opt <- stats::nlminb(rep(0, tp$p), fn, gr, control = list(eval.max = 1000, iter.max = 1000))
  th_k <- theta(opt$par)
  # A parameter only informed by cluster k has gradient zero and would stay at its estimate, as if
  # leaving out the cluster made no difference: check the Hessian at the new estimate
  if(opt$convergence == 0 && is.null(chol_or_null(weighted_hessian(tp, w, th_k)))) {
    stop("The Hessian of the refit without cluster ", label, " is singular or not positive definite: ",
         "some parameters are not or only weakly identified without this cluster (e.g. track-specific parameters), ",
         "so CR3 is not available.")
  }
  list(d = th_k - th, converged = (opt$convergence == 0))
}

# Cholesky factor, or NULL for a (numerically) singular or not positive definite matrix
chol_or_null <- function(H){
  if(anyNA(H) || rcond(H) < 1e-12) return(NULL)
  tryCatch(chol(H), error = function(e) NULL)
}
