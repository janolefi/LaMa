# Internal helpers for qreml() ------------------------------------------------

## Cholesky with a ridge repair, as in mgcv's gam.fit5(). Pivoting reports a
## rank, so indefiniteness is detected rather than just causing a failure.
pd_chol <- function(J, silent = 1) {
  p <- ncol(J)
  D <- diag(J)
  if(any(!is.finite(D))) stop("non-finite values in Hessian")

  repaired <- FALSE
  if(min(D) <= 0){
    small <- max(D) * sqrt(.Machine$double.eps)
    if(-min(D) < small) D[D < small] <- small # only numerically zero, not indefinite
    else repaired <- TRUE
  }

  if(repaired){
    if(silent == 0) cat("Hessian indefinite; adding ridge\n")
    bump <- abs(max(D)) * sqrt(.Machine$double.eps)
    J <- J + diag(abs(min(D)) + bump, p) # sized by the most negative diagonal, never a mean
    d <- rep(1, p) # no preconditioning once a ridge is in
  } else {
    d <- D^-0.5 # unit diagonal, for a stable factorisation
    J <- d * t(d * J)
    bump <- sqrt(.Machine$double.eps)
  }

  R <- suppressWarnings(chol(J, pivot = TRUE))
  while(attr(R, "rank") < p){ # escalate until it factorises at full rank
    repaired <- TRUE
    J <- J + diag(bump, p)
    bump <- bump * 100
    R <- suppressWarnings(chol(J, pivot = TRUE))
  }

  piv <- attr(R, "pivot")
  ipiv <- piv; ipiv[piv] <- seq_len(p)
  logdet <- 2 * sum(log(diag(R))) - 2 * sum(log(d)) # log|J|, free from the factor

  list(R = R, piv = piv, ipiv = ipiv, d = d, logdet = logdet, repaired = repaired)
}

## One diagonal block of J^-1, without forming the full inverse: solve the factor
## against that smooth's indicator columns. Needs p x length(idx), not p x p.
block_inv <- function(fac, idx) {
  p <- length(fac$piv); q <- length(idx)
  E <- matrix(0, p, q)
  E[cbind(idx, seq_len(q))] <- 1
  B <- fac$d * E # J^-1 = Dm Js^-1 Dm, with Dm = diag(d)
  Z <- backsolve(fac$R, forwardsolve(t(fac$R), B[fac$piv, , drop = FALSE]))[fac$ipiv, , drop = FALSE]
  (fac$d * Z)[idx, , drop = FALSE]
}

## Rank, generalised log-determinant and Moore-Penrose inverse of one penalty
## block. Always q x q, never the full p x p penalty.
gen_inverse <- function(S) {
  ev <- eigen(S, symmetric = TRUE)
  keep <- ev$values > max(ev$values) * .Machine$double.eps^0.75
  V <- ev$vectors[, keep, drop = FALSE]
  list(rank = sum(keep), logdet = sum(log(ev$values[keep])),
       inv = V %*% (t(V) / ev$values[keep]))
}


#' Automatic smoothness selection for arbitrary 
#' statistical models with penalised splines, simple random effects, or tensorproducts
#'
#' @description
#' Efficiently fits models involving quadratic penalties of the form
#' \deqn{\sum_{i} \lambda_i b^T S_i b,}
#' by combining automatic differentiation via \code{RTMB} with a custom implementation of the extended Fellner-Schall update.
#' 
#' Users only need to supply a custom penalised log-likelihood function that calls \code{penalty()} or \code{penalty2()}.
#' 
#' @details
#' Step size control follows the extended Fellner-Schall implementation in \code{mgcv} while convergence is judged on the effective degrees of freedom (differing from \code{mgcv}).
#'
#' \strong{The criterion.} The penalty strengths are chosen to maximise the Laplace-approximate restricted log-likelihood
#' \deqn{V(\lambda) = \ell_p(\hat{b}; \lambda) + \tfrac{1}{2} \log |S_\lambda|_+ - \tfrac{1}{2} \log |J|,}{V(lambda) = l_p(b; lambda) + 0.5 log|S_lambda|_+ - 0.5 log|J|,}
#' where \eqn{\ell_p}{l_p} is the penalised log-likelihood at its mode \eqn{\hat{b}}{b}, \eqn{S_\lambda = \sum_i \lambda_i S_i}{S_lambda = sum_i lambda_i S_i}, \eqn{|\cdot|_+}{|.|_+} denotes the product of the positive eigenvalues, and \eqn{J} is the Hessian of the penalised negative log-likelihood at \eqn{\hat{b}}{b}.
#' Every outer iteration refits the inner problem to obtain \eqn{\hat{b}}{b} and \eqn{J}, so one call is a sequence of complete model fits, each warm started at the previous one.
#'
#' \strong{The update.} The penalty strengths are updated multiplicatively, \eqn{\lambda_i \leftarrow \lambda_i r_i}{lambda_i <- lambda_i r_i}, by the extended Fellner-Schall ratio
#' \deqn{r_i = \frac{\mathrm{tr}(S_\lambda^- S_i) - \mathrm{tr}(J^{-1} S_i)}{\hat{b}^T S_i \hat{b}},}{r_i = (tr(S_lambda^- S_i) - tr(J^-1 S_i)) / (b^T S_i b),}
#' which equals one exactly where \eqn{\partial V / \partial \lambda_i = 0}{dV/dlambda_i = 0} and otherwise moves \eqn{\lambda_i} in the direction of the gradient, without ever making it negative.
#' For a smooth carrying a single penalty the first trace is \eqn{\mathrm{rank}(S_i) / \lambda_i}{rank(S_i) / lambda_i}, obtained once from an eigendecomposition instead of at every iteration.
#' Only the diagonal blocks of \eqn{J^{-1}}{J^-1} enter, and these come from the Cholesky factor by triangular solves rather than from a full inverse.
#'
#' \strong{Step control.} The step is taken in \eqn{\log \lambda}{log(lambda)} and scaled by a multiplier that persists across iterations: it doubles while the criterion keeps improving and is halved, up to \code{max_halve} times, while it does not.
#' A step that still worsens the criterion after that many halvings is accepted, which is the only way the iteration can move downhill.
#' In addition no penalty strength may fall by more than a factor \code{alpha} in one iteration, and \eqn{\log \lambda}{log(lambda)} is capped at \code{lsp_max}.
#' The fit returned is the best iterate seen, not necessarily the last one.
#'
#' \strong{Convergence.} The iteration stops once no smooth's effective degrees of freedom
#' \deqn{\mathrm{edf}_i = K_i - \mathrm{tr}(J^{-1} \lambda_i S_i)}{edf_i = K_i - tr(J^-1 lambda_i S_i)}
#' has moved by more than \code{tol_edf} over the last four iterations, provided the step in \eqn{\log \lambda}{log(lambda)} is also small; see \code{tol_edf}.
#' Being eigenvalue sums, the effective degrees of freedom must lie in \eqn{[0, K_i]}, and they leave that range only if the data Hessian is indefinite, in which case the traces driving the update are not the quantities the method assumes.
#' That is reported rather than repaired, through the returned \code{edf_valid} and a warning, and the criterion then falls back to \code{tol} on the restricted log-likelihood.
#'
#' \strong{Numerics.} \eqn{J} is the exact automatic differentiation Hessian, which matters because the finite differencing error of \code{\link[stats:optimHess]{optimHess}} is of the same order as the criterion differences the convergence test has to resolve.
#' It is factorised by a pivoted Cholesky decomposition after scaling to a unit diagonal; should it not be positive definite, a ridge is added and reported through the returned \code{hessian_repaired}.
#'
#' \strong{Uncertainty.} The returned \code{outer_hessian()} gives the Hessian of the negative restricted log-likelihood in the penalty strengths, which \code{\link{sdreport_outer}} turns into standard errors.
#' Note that these describe the curvature of \eqn{V} only, and treat the approximation to it as exact.
#'
#' @seealso \code{\link{penalty}} and \code{\link{penalty2}} to compute the penalty inside the likelihood function, and \code{\link{qreml_old}} for the original algorithm
#'
#' @references Wood, S. N., & Fasiolo, M. (2017). A generalized Fellner-Schall method for smoothing parameter optimization with application to Tweedie location, scale and shape models. Biometrics, 73(4), 1071-1081.
#' @references Koslik, J. O. (2024). Efficient smoothness selection for nonparametric Markov-switching models via quasi restricted maximum likelihood. arXiv preprint arXiv:2411.11498.
#'
#' @param pnll penalised negative log-likelihood function that is structured as dictated by \code{RTMB} and uses the \code{\link{penalty}} or \code{\link{penalty2}} function to compute the penalty
#' @param par named list of initial parameters
#' @param dat initial data list that contains the data used in the likelihood function, hyperparameters, and the \strong{initial penalty strength} vector
#' @param random vector of names of the random effects/ penalised parameters in \code{par}
#'
#' \strong{Caution:} The ordering of \code{random} needs to match the order of the random effects passed to \code{penalty}.
#' @param map optional map argument, containing factor vectors to indicate parameter sharing or fixing
#' @param silent integer silencing level: 0 corresponds to full printing of inner and outer iterations, 1 to printing of outer iterations only, and 2 to no printing
#' @param spname optional name given to the penalty strength parameter in \code{dat}. Defaults to \code{"lambda"}
#' @param tol_edf convergence tolerance on the \strong{effective degrees of freedom}. Defaults to 0.001.
#' 
#' This is the primary convergence criterion. The iteration stops once no smooth's effective degrees of freedom has changed by more than \code{tol_edf} across the last four outer iterations, provided the step in \code{log(lambda)} is also small.
#' The effective degrees of freedom are what say whether the \strong{fitted smooth} is still changing, and unlike a tolerance on the restricted likelihood they mean the same thing on every model, being measured in effective parameters rather than in nats.
#' The penalty strengths themselves are a poor proxy: \code{lambda} can slide along a flat ridge for many iterations, moving substantially in relative terms, while neither the criterion nor the fit changes appreciably.
#' @param maxiter maximum number of outer iterations
#' @param alpha smallest factor by which a penalty strength may \strong{decrease} in one outer iteration, a number in [0, 1). Defaults to 0.1.
#'
#' Penalty strengths are free to increase as fast as the update proposes, but cannot collapse faster than this per iteration.
#' Reducing a penalty strength too quickly can push the inner optimisation into a local optimum or a numerically awkward region, which matters more here than in a GAM because the likelihood is user-written.
#' Set to zero to remove the floor entirely. Step length is handled separately, by \code{max_halve} and the adaptive multiplier.
#' @param tol \strong{fallback} convergence tolerance, on the restricted log-likelihood. Defaults to 0.01.
#'
#' Used only when the effective degrees of freedom are not trustworthy, i.e.\ when they fall outside \eqn{[0, K_i]}, which happens when the data Hessian is indefinite (see \code{tol_edf} and the returned \code{edf_valid}).
#' The iteration then stops once the restricted log-likelihood has changed by less than \code{tol} over the last four outer iterations.
#' A tolerance in nats is not comparable across models, which is why it is the fallback rather than the primary criterion.
#' @param smoothing optional scaling factor for the final penalty strength parameters. Increasing this beyond one leads to a smoother final model
#' @param method inner optimisation method to be used by \code{\link[stats:optim]{optim}}. Defaults to \code{"BFGS"}
#' @param lsp_max largest value allowed for \code{log(lambda)}. Defaults to 15, as in \code{mgcv}, i.e. penalty strengths saturate at roughly 3.3e6
#' @param step_small size of a step in \code{log(lambda)} below which the step multiplier is allowed to double. Defaults to 0.05, as in \code{mgcv}
#' @param max_halve maximum number of times a step that decreases the restricted likelihood is halved before it is accepted anyway. Defaults to 6.
#'
#' \code{mgcv} never shortens below the full Fellner-Schall step and accepts a worse one instead, which for a user-written likelihood can drift downhill for tens of iterations.
#' @param control list of control parameters for \code{\link[stats:optim]{optim}} to use in the inner optimisation
#' @param spHess logical, if \code{TRUE}, the sparse automatic differentiation Hessian is used for evaluation. The factorisation is dense either way
#' @param joint_unc logical, if \code{TRUE}, joint \code{RTMB} object is returned allowing for joint uncertainty quantification
#' @param saveall logical, if \code{TRUE}, then all model objects from each iteration are saved in the final model object
#'
#' @return model object of class \code{"qremlModel"}, carrying the fitted quantities, the smoothness selection diagnostics shown by \code{\link{summary.qremlModel}}, and \code{outer_hessian()}, the outer Hessian that \code{\link{sdreport_outer}} turns into standard errors
#'
#' @export
#'
#' @import RTMB
#'
#' @examples
#' data = trex[1:2000,] # subset of the data
#'
#' # initial parameter list
#' par = list(logmu = log(c(0.3, 2.5)), # step mean
#'            logsigma = log(c(0.3, 1.5)), # step sd
#'            beta0 = c(-2,-2), # state process intercept
#'            beta_spline = matrix(rep(0, 18), nrow = 2)) # state process spline coefs
#'           
#' # data object with initial penalty strength lambda
#' dat = list(step = data$step, # step length
#'            tod = data$tod, # time of day covariate
#'            N = 2, # number of states
#'            lambda = rep(20,2)) # initial penalty strength
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
#'   Gamma = tpm(cbind(beta0, beta_spline), Z) # transition probabilities
#'   delta = stationary_p(Gamma, t = 1) # initial distribution
#'   mu = exp(logmu) # step mean
#'   sigma = exp(logsigma) # step sd
#'   # calculating all state-dependent densities
#'   allprobs = matrix(1, nrow = length(step), ncol = N)
#'   ind = which(!is.na(step)) # only for non-NA obs.
#'   for(j in 1:N) allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j])
#'   -forward(delta, Gamma[,,tod], allprobs) +
#'       penalty(beta_spline, S, lambda) # this does all the penalization work
#' }
#'
#' # model fitting
#' mod = qreml(pnll, par, dat, random = "beta_spline", silent = 2)
qreml <- function(pnll, # penalised negative log-likelihood function
                  par, # initial parameter list
                  dat, # initial dat object, currently needs to be called dat!
                  random, # names of parameters in par that are random effects/ penalised
                  map = NULL, # map for fixed effects
                  silent = 1, # print level
                  spname = "lambda", # name given to the smoothing parameter parameter in dat
                  tol_edf = 0.001, # convergence tolerance on the effective degrees of freedom
                  maxiter = 100, # maximum number of outer iterations
                  alpha = 0.1, # smallest factor by which lambda may decrease per iteration
                  tol = 0.01, # fallback tolerance on the restricted log-likelihood
                  smoothing = 1,
                  method = "BFGS", # optimisation method used by optim
                  lsp_max = 15, # largest allowed log(lambda)
                  step_small = 0.05, # step size below which the multiplier may double
                  max_halve = 6, # how often a worsening step may be halved
                  control = list(), # control list for inner optimisation
                  spHess = FALSE, # evaluate the Hessian sparsely
                  joint_unc = FALSE, # should joint object be returned?
                  saveall = FALSE) # save all intermediate models?
{
  ### Input checking
  if(!is.function(pnll)) stop("'pnll' needs to be a function")
  if(!is.list(par)) stop("'par' needs to be a named list")
  if(!is.list(dat)) stop("'dat' needs to be a named list")
  if(!spname %in% names(dat)){
    stop(paste0("'dat' needs to contain a vector called '", spname, "' with initial penalty strengths"))
  }
  if(!is.character(random) || length(random) < 1){
    stop("'random' needs to be a character vector of names of random effects in 'par'")
  }
  if(!is.null(map) && !is.list(map)){
    stop("'map' needs to be a named list of factors for fixed effects or penalty strength parameters")
  }
  if(!is.numeric(alpha) || length(alpha) != 1 || alpha < 0 || alpha >= 1){
    stop("'alpha' needs to be a single number in [0, 1)")
  }
  if(!is.numeric(max_halve) || length(max_halve) != 1 || max_halve < 0){
    stop("'max_halve' needs to be a single non-negative number")
  }
  if(!is.numeric(tol_edf) || length(tol_edf) != 1 || tol_edf < 0){
    stop("'tol_edf' needs to be a single non-negative number")
  }

  # setting the argument names because later updated par is returned
  argname_par <- as.character(substitute(par))
  argname_dat <- as.character(substitute(dat))

  n_re <- length(random) # number of distinct random effects
  allmods <- list() # list to save all model objects

  # initialising penalty strength lambda
  lambda <- dat[[spname]]
  lambda0 <- lambda # so that fixed parts can be refilled when 'lambda' is changed

  # creating the objective function as wrapper around pnll to pull lambda from local
  f <- function(par){
    environment(pnll) = environment()

    # overloading assignment operators, currently necessary
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")

    # defining function that grabs lambda
    getLambda <- function(x) lambda
    dat[[spname]] <- DataEval(getLambda, rep(advector(1), 0))

    # assigning dat to whatever it is called in pnll()
    assign(argname_dat, dat, envir = environment())
    pnll(par)
  }

  ## mapping
  if(!is.null(map)){
    if(any(names(map) %in% random)){
      stop("'map' cannot contain random effects or spline parameters")
    }
    map <- lapply(map, factor)
  }
  if(is.null(map[[spname]])) map[[spname]] <- factor(seq_along(lambda))
  lambda_map <- map[[spname]]
  if(length(lambda_map) != length(lambda)){
    stop(paste0("Length of map argument for ", spname, " has wrong length."))
  }
  map <- map[names(map) != spname]
  if(length(map) == 0) map <- NULL

  # local versions of the mapping helpers, so only the vector needs passing
  map_lambda <- function(lambda){ # tied entries collapse to one value per level
    grp <- split(lambda, lambda_map)
    if(is.character(lambda)) vapply(grp, paste, character(1), collapse = "&")
    else vapply(grp, mean, numeric(1))
  }
  unmap_lambda <- function(lambda_mapped){ # entries mapped to NA stay fixed
    lambda <- lambda0
    free <- !is.na(lambda_map)
    lambda[free] <- lambda_mapped[as.integer(lambda_map[free])]
    lambda
  }

  lambda_mapped <- map_lambda(lambda)
  # with every penalty strength fixed there is nothing for the outer iteration to
  # do, and it must be skipped rather than run once: the step is over an empty
  # vector, which is not a convergence failure
  estimate_ps <- length(lambda_mapped) > 0
  if(!estimate_ps) message("No penalty parameters will be estimated as all are fixed.")

  ## creating the RTMB objective function
  if(silent %in% 0:1) message("Creating AD function")
  obj <- MakeADFun(func = f, parameters = par, silent = TRUE, map = map)
  newpar <- obj$par

  ## Choosing how the Hessian is evaluated (3 potential routes)
  # 1) default: obj$he(); the exact AD Hessian for accuracy
  # 2) fallback only: optimHess()
  # 3) sparse Hessian: spHess via RTMB's sparse Hessian route
  # -> fairly expensive to build but well worth it if model has a sparse Hessian
  if(spHess){ # sparse option
    Tape <- RTMB::GetTape(obj, name = "ADFun")
    if(silent < 2) message("Constructing sparse Hessian")
    spH <- Tape$jacfun(sparse = TRUE)$jacfun(sparse = TRUE)
    rm(Tape)
    gc(verbose = FALSE)
    hessian_at <- function(p) as.matrix(spH(p))
  } else if(!is.null(obj$he)){ # regular option
    hessian_at <- function(p) obj$he(p)
  } else { # fallback
    if(silent < 2) message("obj$he() unavailable, falling back to finite differences")
    hessian_at <- function(p) stats::optimHess(p, obj$fn, obj$gr)
  }

  ## Gradient printing (only when silent == 0)
  counter_env <- new.env()
  counter_env$count <- 0
  if(silent == 0){
    ctREPORT <- 10
    if(!is.null(control$REPORT)){
      ctREPORT <- control$REPORT
      control$REPORT <- NULL
    }
    newgrad <- function(par){
      counter_env$count <- counter_env$count + 1
      ct <- counter_env$count
      gr <- obj$gr(par)
      if(ct %% ctREPORT == 0) cat("iter", ct, "- inner mgc:", round(max(abs(gr)), 5), "\n")
      gr
    }
  } else newgrad <- obj$gr

  ## prepwork -> running reporting to get necessary quantities
  mod0 <- obj$report()
  S <- mod0$S

  # index of each random effect's coefficients in obj$par, one smooth per row:
  # re_inds[[i]][j, ] are the coefficients of the j-th smooth
  re_inds <- lapply(seq_len(n_re), function(i){
    b <- par[[random[i]]]
    dims <- if(is.null(dim(b))) c(1, length(b)) else dim(b) # a vector is one smooth
    if(length(dims) != 2) stop(random[i], " must be a vector or matrix")
    inds <- matrix(which(names(obj$par) == random[i]), dims[1], dims[2]) # RTMB flattens by column
    # a smooth has as many coefficients as its penalty matrix has rows; when those
    # sit in the columns of b, transpose so that rows index smooths either way
    S_i <- if(is.matrix(S[[i]])) S[[i]] else S[[i]][[1]] # a list of S = tensor product
    if(dims[1] == nrow(S_i)) t(inds) else inds
  })

  ## how many penalty strengths per random effect: 1 = simple smooth, >1 = tensor product
  n_penalties <- sapply(S, function(x) if(is.matrix(x)) 1 else length(x))
  simple_ind <- which(n_penalties == 1)
  tp_ind <- which(n_penalties > 1)

  re_lengths <- sapply(re_inds, function(x) if(is.vector(x)) 1 else nrow(x))
  lambda_lengths <- n_penalties * re_lengths
  if(length(lambda) != sum(lambda_lengths)){
    stop(paste0("Length of '", spname, "' does not match the number of penalty strength parameters needed"))
  }

  # splits a lambda vector into one entry per smooth
  reshape_lambda <- function(lambda){
    to <- cumsum(lambda_lengths)
    Map(function(from, to) lambda[from:to], to - lambda_lengths + 1, to)
  }

  ## rank and generalised log-determinant of each simple penalty matrix
  # both come from one eigen-decomposition, computed once: for a simple smooth
  # S_lambda = lambda * S, so log|S_lambda|_+ = rank * log(lambda) + log|S|_+ and
  # tr(S_lambda^- S) = rank / lambda, neither of which needs redoing per iteration
  ranks <- rep(NA_real_, n_re)
  logdetS0 <- rep(NA_real_, n_re)
  for(i in simple_ind){
    gi <- gen_inverse(S[[i]])
    ranks[i] <- gi$rank
    logdetS0[i] <- gi$logdet
  }

  ## naming the penalty strengths
  Lambda0 <- reshape_lambda(lambda)
  for(ind in seq_along(simple_ind)){
    names(Lambda0[simple_ind][[ind]]) <- seq_along(Lambda0[simple_ind][[ind]])
  }
  for(ind in seq_along(tp_ind)){
    margin_names <- names(S[[tp_ind[ind]]])
    names(Lambda0[tp_ind][[ind]]) <- paste0(rep(1:re_lengths[tp_ind[ind]], each = length(margin_names)), ".",
                                            rep(margin_names, re_lengths[tp_ind[ind]]))
  }
  lambda_names <- names(unlist(Lambda0))

  ## penalty block of one smooth, i.e. lambda * S or sum_k lambda_k S_k
  block_S <- function(i, j, Lambda) {
    if(i %in% simple_ind) return(Lambda[[i]][j] * S[[i]])
    n_pen <- length(S[[i]])
    out <- Lambda[[i]][(j-1) * n_pen + 1] * S[[i]][[1]]
    for(pen in 2:n_pen) out <- out + Lambda[[i]][(j-1) * n_pen + pen] * S[[i]][[pen]]
    out
  }

  ## log|S_lambda|_+ summed over all smooths, block by block
  logdet_Slambda <- function(Lambda) {
    out <- 0
    for(i in seq_len(n_re)){
      for(j in seq_len(nrow(re_inds[[i]]))){
        if(i %in% simple_ind){
          out <- out + ranks[i] * log(Lambda[[i]][j]) + logdetS0[i]
        } else {
          out <- out + gen_inverse(block_S(i, j, Lambda))$logdet
        }
      }
    }
    out
  }

  # number of coefficients per smooth, the upper limit for that smooth's edf
  block_dim <- unlist(lapply(seq_len(n_re),
                             function(i) rep(ncol(re_inds[[i]]), nrow(re_inds[[i]]))))
  edf_hist <- matrix(NA_real_, length(block_dim), maxiter) # one column per iteration

  ## controlling optim
  ctl <- list(maxit = 1000)
  ctl[names(control)] <- control
  if(method == "BFGS") ctl$reltol <- 1e-10
  if(method == "L-BFGS-B") ctl$maxit <- 5000

  ## one complete evaluation at a given log penalty strength:
  ## inner fit, Hessian, factorisation and criterion
  fit_at <- function(lsp_try, start) {
    lambda <<- unmap_lambda(exp(lsp_try))
    Lambda <- reshape_lambda(lambda)

    if(silent == 0) cat("\nInner optimisation at", spname, "=", round(exp(lsp_try), 3), "\n")
    counter_env$count <- 0

    # RTMB caches the objective value at the last parameter vector it saw, which is
    # exactly 'start'. Without this nudge optim's first evaluation returns the value
    # belonging to the PREVIOUS lambda and it stops there.
    invisible(try(obj$fn(start + 1e-8 * (1 + abs(start))), silent = TRUE))

    opt <- stats::optim(start, obj$fn, newgrad, method = method, control = ctl)
    gr <- obj$gr(opt$par) # forces evaluation at opt$par so that report() matches
    if(silent == 0) cat("iter", counter_env$count, "- inner mgc:", round(max(abs(gr)), 5), "\n")

    mod <- obj$report()
    if(silent == 0) cat("evaluating Hessian...\n")
    J <- hessian_at(opt$par)
    fac <- pd_chol((J + t(J)) / 2, silent) # force symmetric
    llk_r <- -opt$value + logdet_Slambda(Lambda) / 2 - fac$logdet / 2 # criterion is -llk_r

    # J is not kept: everything downstream works from the factor, and cur, trial
    # and trial2 can be alive at once
    list(opt = opt, mod = mod, fac = fac, Lambda = Lambda,
         lsp = lsp_try, llk_r = llk_r, crit = -llk_r, mgc = max(abs(gr)))
  }

  ## Fellner-Schall ratio on the mapped scale: lambda_new = lambda * r, so the step
  ## in log(lambda) is log(r). Also returns the pieces of the outer gradient and the
  ## effective degrees of freedom, which the same inverse blocks give for free.
  efs_ratio <- function(state) {
    a <- bSb <- rep(NA_real_, length(lambda0))
    edf <- numeric(length(block_dim))
    l <- 1; bl <- 1
    for(i in seq_len(n_re)) for(j in seq_len(nrow(re_inds[[i]]))){
      idx <- re_inds[[i]][j, ]
      Jinv <- block_inv(state$fac, idx) # only the block we need
      edf[bl] <- length(idx) - sum(Jinv * block_S(i, j, state$Lambda)); bl <- bl + 1

      if(i %in% simple_ind){
        a[l] <- ranks[i] / state$Lambda[[i]][j] - sum(Jinv * S[[i]]) # tr(S_lambda^- S) = rank/lambda
        bSb[l] <- state$mod$Pen[[i]][j]
        l <- l + 1
      } else {
        n_pen <- length(S[[i]])
        Sinv <- gen_inverse(block_S(i, j, state$Lambda))$inv
        for(k in seq_len(n_pen)) a[l+k-1] <- sum(Sinv * S[[i]][[k]]) - sum(Jinv * S[[i]][[k]])
        bSb[l:(l+n_pen-1)] <- state$mod$Pen[[i]][[j]] # penalty2() reports b'S_k b, without lambda_k
        l <- l + n_pen
      }
    }

    # only positivity is enforced: mgcv's sqrt(eps) floor would clamp both for a
    # strongly penalised smooth and freeze its lambda at r = 1
    a <- pmax(.Machine$double.xmin, a)
    bSb <- pmax(.Machine$double.xmin, bSb)

    grp <- split(seq_along(a), lambda_map) # indices of each tied group
    a_m <- vapply(grp, function(ii) sum(a[ii]), 0)
    bSb_m <- vapply(grp, function(ii) sum(bSb[ii]), 0)

    r <- a_m / bSb_m
    r[!is.finite(r)] <- 1e6
    r <- pmin(pmax(r, 1e-6), 1e6) # no step may move lambda by more than a factor 1e6
    list(r = unname(r), a = unname(a_m), bSb = unname(bSb_m), edf = edf)
  }

  ### updating algorithm
  lsp <- log(lambda_mapped)
  mult <- 1 # step multiplier, persistent across iterations, as in mgcv
  n_success <- 0 # consecutive improving steps
  crit_hist <- rep(NA_real_, maxiter)

  if(silent < 2) message("Initialising with ", spname, ": ", paste(round(lambda, 3), collapse = " "))
  if(silent == 0) cat("\nouter 0 - initial fit\n")

  # one fit at the starting lambda is unavoidable: the first update needs the
  # penalties and traces at a converged inner solution
  cur <- fit_at(lsp, newpar)
  converged <- FALSE

  best_lsp <- lsp # best iterate seen, in case a run ends worse than it once was
  best_crit <- cur$crit
  best_iter <- 0
  edf_ok_run <- TRUE # flips once the edf leave [0, block dimension]

  if(!estimate_ps){
    crit_hist[1] <- cur$crit
    converged <- TRUE
  }

  for(iter in seq_len(if(estimate_ps) maxiter else 0)){

    if(saveall) allmods[[iter]] <- cur$mod

    ef <- efs_ratio(cur) # compute edf ratio quantities
    step <- log(ef$r) # compute proposed step on log-scale

    # the floor applies to the step actually taken, mult * step: applied to the raw
    # step it would become alpha^mult, i.e. no guarantee when moving fastest
    lo <- if(alpha > 0) log(alpha) else -Inf # lambda may not fall faster than alpha
    take <- function(m) pmin(lsp + pmax(m * step, lo), lsp_max)

    lsp1 <- take(mult) # updated smoothing parameters
    max_step <- max(abs(lsp1 - lsp))
    n_halve <- 0
    if(silent == 0) cat("\nouter", iter, "- proposed", paste0(spname, ":"), round(exp(lsp1), 3), "\n")

    trial <- fit_at(lsp1, cur$opt$par) # trial fit at new smoothing parameters

    if(trial$crit <= cur$crit){ ## improved
      n_success <- n_success + 1
      accelerated <- FALSE
      # mgcv gates acceleration on the largest step across ALL lambda, so one busy
      # parameter throttles the rest; two improving steps in a row opens it too
      if(max_step < step_small || n_success >= 2){
        lsp2 <- take(2 * mult)
        trial2 <- fit_at(lsp2, cur$opt$par)
        if(trial2$crit < trial$crit){
          trial <- trial2; lsp1 <- lsp2; mult <- 2 * mult; accelerated <- TRUE
          if(silent == 0) cat("accelerating: step multiplier now", mult, "\n")
        }
      }
      # must sit outside the branch above: once mult is small the steps are too, so
      # a rejected doubling would leave it collapsed for the rest of the run
      if(!accelerated && mult < 1) mult <- min(2 * mult, 1)
    } else { ## worsened
      n_success <- 0
      # mgcv accepts a worse step at mult = 1; here that can drift downhill for
      # tens of iterations, so it is genuinely backtracked
      while(trial$crit > cur$crit && n_halve < max_halve){
        mult <- mult / 2
        n_halve <- n_halve + 1
        if(silent == 0) cat("criterion increased; step multiplier now", signif(mult, 4), "\n")
        lsp1 <- take(mult)
        trial <- fit_at(lsp1, cur$opt$par)
      }
    }

    lsp <- lsp1
    cur <- trial
    crit_hist[iter] <- cur$crit
    if(cur$crit < best_crit){
      best_crit <- cur$crit
      best_lsp <- lsp
      best_iter <- iter
    }

    if(silent == 1){
      cat("outer", iter, "-", paste0(spname, ":"), round(exp(lsp), 3), "\n")
    } else if(silent == 0){
      if(n_halve > 0) cat("outer", iter, "- accepted", paste0(spname, ":"), round(exp(lsp), 3), "\n")
      cat("outer", iter, "- restricted llk:", round(cur$llk_r, 5),
          "- max step:", round(max_step, 5), "\n")
    }

    #### Convergence check ####
    # Neither the smoothing parameters nor the restricted log-likelihood are good 
    # proxies for "the fit has stopped changing"; lambda can slide along a flat 
    # ridge for a hundred iterations while the criterion barely moves.
    #
    # The effective degrees of freedom directly tell us whether the fitted model
    # is still changing and come free from the inverse blocks above. Converge
    # once no smooth's edf has moved by more than tol_edf.
    #
    # The edf are eigenvalue sums of J^-1 H and must lie in [0, block dimension].
    # Outside that range the data Hessian is indefinite, which means the traces
    # driving the update are not the quantities the algorithm assumes, so it is
    # worth saying so at the iteration where it starts rather than only at the
    # end, after minutes of iterating on them.
    if(edf_ok_run && any(ef$edf < -1e-3 | ef$edf > block_dim + 1e-3)){
      edf_ok_run <- FALSE
      if(silent < 2){
        message("Iteration ", iter, ": effective degrees of freedom left [0, block ",
                "dimension] (min ", signif(min(ef$edf), 4), "); the data Hessian is ",
                "indefinite and the update is running on unreliable traces")
      }
    }

    edf_hist[, iter] <- ef$edf
    gc(verbose = FALSE) # each iteration allocates several p x p matrices

    # Stop when the fitted smooths have settled.
    # Gates: four iterations, a small step, and a window rather than one 
    # consecutive pair.
    if(iter > 3 && max_step < step_small){
      if(edf_ok_run){
        w <- t(edf_hist[, (iter-3):iter, drop = FALSE]) # iterations down the rows
        if(max(abs(diff(w))) < tol_edf){
          converged <- TRUE
          if(silent == 0) cat("effective degrees of freedom settled\n")
        }
      } else if(max(abs(diff(crit_hist[(iter-3):iter]))) < tol){
        converged <- TRUE # fallback: edf unreliable, so use the restricted likelihood
        if(silent == 0) cat("restricted likelihood settled (edf unreliable)\n")
      }
    }
    if(converged){
      if(silent < 2) message("Converged")
      break
    }

    if(iter == maxiter){
      message("No convergence")
      warning("No convergence\n")
    }
  }

  # a zero-length for() sets its loop variable to NULL, so restore it
  if(!estimate_ps) iter <- 1

  ## final model fit, at the best penalty strengths found
  # once converged the criterion wobbles, so the last iterate is rarely the very
  # best one and saying so every time is noise. 'tol' is the scale on which a
  # change in the restricted likelihood counts as negligible, so use it here too
  shortfall <- crit_hist[iter] - best_crit
  if(shortfall > tol && silent < 2){
    message("Returning the ", spname, " from iteration ", best_iter,
            "; the run ended ", signif(shortfall, 4), " worse than that")
  }
  lsp <- best_lsp
  lambda <- unmap_lambda(exp(lsp)) * smoothing
  if(silent < 2){
    if(any(smoothing != 1)) message("Smoothing factor: ", paste(smoothing, collapse = " "))
    message("Final model fit with ", spname, ": ", paste(round(lambda, 3), collapse = " "))
  }

  final_lsp <- log(map_lambda(lambda))
  # cur is already the fit at this lambda unless smoothing was applied or the last
  # iterate was not the best one, so refitting would repeat a fit for nothing
  final <- if(isTRUE(all.equal(final_lsp, cur$lsp))) cur else fit_at(final_lsp, cur$opt$par)
  mod <- final$mod
  opt <- final$opt
  Lambda <- final$Lambda
  if(saveall) allmods[[iter + 1]] <- mod

  pllk <- -opt$value # penalised log-likelihood
  llk <- pllk + mod$pen

  mod$obj <- obj
  if(saveall) mod$allmods <- allmods

  names(lambda) <- lambda_names
  mod[[spname]] <- lambda
  mod[[paste0("all_", spname)]] <- Lambda

  parlist <- obj$env$parList(opt$par)
  mod[[argname_par]] <- parlist
  mod[[paste0("relist_", argname_par)]] <- obj$env$parList
  mod[[paste0("map_", spname)]] <- map_lambda
  mod$spname <- spname
  mod$parname <- argname_par
  mod[[paste0(argname_par, "_vec")]] <- opt$par
  mod$llk <- llk
  mod$n_fixpar <- length(unlist(par[!(names(par) %in% random)]))

  ## Effective degrees of freedom, from the same inverse blocks: summing
  ## diag(J^-1 H) over a smooth gives q - sum(J^-1[idx,idx] * S_lambda[idx,idx]),
  ## since S_lambda is block diagonal and no off-block entry is needed
  Edfs <- Lambda # just to copy the nested list structure
  for(i in seq_len(n_re)){
    if(i %in% tp_ind) Edfs[[i]] <- numeric(nrow(re_inds[[i]]))
    for(j in seq_len(nrow(re_inds[[i]]))){
      idx <- re_inds[[i]][j, ]
      Jinv <- block_inv(final$fac, idx)
      Edfs[[i]][j] <- length(idx) - sum(Jinv * block_S(i, j, Lambda))
    }
  }
  mod$df <- mod$n_fixpar + sum(unlist(Edfs)) # one number for the full model
  mod$edf <- Edfs # list of individual edfs

  ## sanity check: each edf must lie in [0, K_i], being a sum of eigenvalues of
  ## J^-1 H. Outside that the data Hessian is indefinite and the inverse is not
  ## trustworthy. Reported rather than repaired: a ridge big enough to fix an
  ## indefinite H also smooths away the weakly identified directions.
  edf_ok <- TRUE
  for(i in seq_len(n_re)){
    q_i <- ncol(re_inds[[i]])
    if(any(Edfs[[i]] < -1e-3 | Edfs[[i]] > q_i + 1e-3)) {
      edf_ok <- FALSE
    } 
  }
  mod$edf_valid <- edf_ok

  ## lambda pinned at the upper bound: the smooth is shrunk to its null space
  ## because the bound was binding, not because the data asked for it
  at_bound <- which(log(lambda) >= lsp_max - 1e-8)
  mod$lambda_at_bound <- at_bound
  if(length(at_bound) > 0 && silent < 2){
    message(length(at_bound), " of ", length(lambda), " ", spname,
            " ended at the upper bound exp(lsp_max) = ", signif(exp(lsp_max), 4),
            "; those smooths are penalised to their null space")
  }
  if(!edf_ok){
    lam_span <- range(lambda[lambda > 0])
    warning("Effective degrees of freedom outside [0, block dimension]: the inverse ",
            "Hessian is unreliable and the fit should not be trusted. The penalty ",
            "strengths span ", signif(lam_span[1], 3), " to ", signif(lam_span[2], 3),
            "; refitting with a smaller 'lsp_max' usually resolves it.", call. = FALSE)
  }

  if(!is.null(mod$allprobs)) mod$nobs <- nrow(mod$allprobs)

  mod$Hessian_conditional <- hessian_at(final$opt$par)
  mod$llk_restricted <- -crit_hist[seq_len(iter)]

  ## the iteration is meant to be monotone -- only an exhausted backtrack is
  ## accepted -- and this is the number that says whether that held
  mod$max_drop <- if(iter > 1) max(c(0, diff(crit_hist[seq_len(iter)]))) else 0
  if(mod$max_drop > tol && silent < 2){
    message("Restricted likelihood decreased by up to ", signif(mod$max_drop, 4),
            " during the iteration; the step control could not keep it monotone")
  }
  mod$tol <- tol # the scale on which a change in the criterion is negligible
  mod$converged <- converged
  mod$iter <- iter
  mod$best_iter <- best_iter
  mod$hessian_repaired <- final$fac$repaired

  # dV/dlog(lambda) = 0.5 * lambda * (a - b'Sb), zero exactly where the update has
  # a fixed point (a = b'Sb, i.e. r = 1)
  fr <- efs_ratio(final)
  mod$outer_grad <- 0.5 * map_lambda(lambda) * (fr$a - fr$bSb)
  names(mod$outer_grad) <- levels(lambda_map)

  # the same gradient, w.r.t. lambda itself and as a function of it. Kept local:
  # only the Hessian below is reported. Every call refits the inner problem, since
  # a and b'Sb are only defined at the conditional mode.
  outer_gr <- function(x){
    # started cold, not at opt$par: from the fitted mode optim's relative tolerance
    # fires before the mode has moved, and the differences taken below are far too
    # small to survive that, which inflates the outer Hessian
    fr <- efs_ratio(fit_at(log(x), newpar))
    0.5 * (fr$a - fr$bSb)
  }

  # Hessian of the negative restricted log-likelihood, as Hessian_conditional is
  # for the inner problem. lambda_hat is captured here because fit_at() overwrites
  # 'lambda' as it differences.
  lambda_hat <- map_lambda(lambda)
  mod$outer_hessian <- function(x = lambda_hat){
    H <- -numDeriv::jacobian(outer_gr, x, method = "simple")
    (H + t(H)) / 2 # one-sided differences leave it slightly asymmetric
  }

  # removing elements only reported for the update
  mod <- mod[!names(mod) %in% c("Pen", "pen", "S")]

  
  ## If the model is simple (i.e. contains no tensor-products), set up full REML
  # marginal likelihood which can be used for *joined* uncertainty quantification 
  # (as opposed to conditional on the fitted smoothing parameters) via sdreport()
  if(length(tp_ind) == 0 && joint_unc){
    ### constructing joint object
    parlist$loglambda <- log(mod[[spname]])
    logdetS <- numeric(length(S))
    for(i in seq_along(S)) logdetS[i] <- gdeterminant(S[[i]])

    # joint negative log-likelihood -> pnll with log normalisation constants 
    jnll <- function(par) {
      environment(pnll) = environment()
      "[<-" <- ADoverload("[<-")
      "c" <- ADoverload("c")
      "diag<-" <- ADoverload("diag<-")
      dat[[spname]] <- exp(par$loglambda)
      l_p <- -pnll(par[names(par) != "loglambda"])
      const <- 0
      for(i in 1:n_re){
        for(j in 1:nrow(re_inds[[i]])){
          k <- length(re_inds[[i]][j, ])
          loglam <- if(i == 1) par$loglambda[j] else par$loglambda[re_lengths[i-1] + j]
          const <- const - k * log(2*pi) + k * loglam + logdetS[i]
        }
      }
      -(l_p + 0.5 * const)
    }

    if(is.null(map)) map <- list(loglambda = lambda_map) else map$loglambda <- lambda_map
    mod$obj_joint <- MakeADFun(jnll, parlist,
                               random = names(par)[names(par) != "loglambda"],
                               map = map)
  }

  class(mod) <- "qremlModel"
  return(mod)
}
