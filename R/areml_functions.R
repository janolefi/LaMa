# Internal helpers for areml() ------------------------------------------------

## Positive definite Cholesky factorisation with repair
#
# Follows the strategy used by mgcv's gam.fit5(): the matrix is first checked
# for indefiniteness via its diagonal, repaired by a ridge whose size is set by
# the most negative diagonal element, preconditioned to unit diagonal, and then
# factorised with a *pivoted* Cholesky. Pivoting reports a rank, so a rank
# deficient matrix is detected rather than merely failing, and the ridge is
# escalated until the factorisation has full rank.
#
# Returns the factor together with everything needed to solve with it later.
pd_chol <- function(J, silent = 1) {
  p <- ncol(J)
  D <- diag(J)
  if(any(!is.finite(D))) stop("non-finite values in Hessian")

  repaired <- FALSE

  if(min(D) <= 0) {
    # a diagonal entry that is only numerically zero is not indefiniteness
    Dthresh <- max(D) * sqrt(.Machine$double.eps)
    if(-min(D) < Dthresh) {
      D[D < Dthresh] <- Dthresh
    } else {
      repaired <- TRUE
    }
  }

  if(repaired) {
    if(silent == 0) cat("Hessian indefinite; adding ridge\n")
    # ridge = size of the most negative diagonal, plus a small one set by the
    # largest. Never a mean, which can itself be negative.
    Ip <- abs(max(D)) * sqrt(.Machine$double.eps)
    J <- J + diag(abs(min(D)) + Ip, p)
    d <- rep(1, p) # no preconditioning once a ridge has been added
  } else {
    d <- D^-0.5 # scale to unit diagonal: this is what keeps the factorisation stable
    J <- d * t(d * J)
    Ip <- sqrt(.Machine$double.eps)
  }

  R <- suppressWarnings(chol(J, pivot = TRUE))
  while(attr(R, "rank") < p) { # escalate until it factorises at full rank
    repaired <- TRUE
    J <- J + diag(Ip, p)
    Ip <- Ip * 100
    R <- suppressWarnings(chol(J, pivot = TRUE))
  }

  piv <- attr(R, "pivot")
  ipiv <- piv
  ipiv[piv] <- seq_len(p)

  # log|J|, free from the factor: pivoting does not change the determinant, and
  # the unit diagonal scaling contributes a known factor
  logdet <- 2 * sum(log(diag(R))) - 2 * sum(log(d))

  list(R = R, piv = piv, ipiv = ipiv, d = d, logdet = logdet, repaired = repaired)
}

## Diagonal block of the inverse, without ever forming the full inverse
#
# Solves the (already computed) factor against the indicator columns of one
# smooth's coefficients, which costs one triangular solve pair per penalised
# coefficient and needs p x length(idx) memory rather than p x p. Everything the
# update and the effective degrees of freedom need lives in this block.
block_inv <- function(fac, idx) {
  p <- length(fac$piv)
  q <- length(idx)
  E <- matrix(0, p, q)
  E[cbind(idx, seq_len(q))] <- 1

  # J = Dm^-1 Js Dm^-1 with Dm = diag(d), so J^-1 = Dm Js^-1 Dm
  B <- fac$d * E
  Z <- backsolve(fac$R, forwardsolve(t(fac$R), B[fac$piv, , drop = FALSE]))[fac$ipiv, , drop = FALSE]
  Z <- fac$d * Z
  Z[idx, , drop = FALSE]
}

## Generalised eigen-decomposition of one penalty block
#
# Returns the rank, the generalised log-determinant (product of the non-zero
# eigenvalues) and the Moore-Penrose inverse. Used once per smooth at setup for
# simple smooths, and once per iteration for tensor products, always on the
# small q x q block rather than the full p x p penalty matrix.
gen_inverse <- function(S) {
  e <- eigen(S, symmetric = TRUE)
  ev <- e$values
  keep <- ev > max(ev) * .Machine$double.eps^0.75
  V <- e$vectors[, keep, drop = FALSE]
  list(rank = sum(keep),
       logdet = sum(log(ev[keep])),
       inv = V %*% (t(V) / ev[keep]))
}


#' Accelerated restricted maximum likelihood (aREML) for models with penalised splines or simple i.i.d. random effects
#'
#' @description
#' Fits statistical models involving \strong{penalised splines} or simple \strong{i.i.d. random effects}, i.e. that have penalties of the form
#' \deqn{0.5 \sum_{i} \lambda_i b_i^T S_i b_i,}
#' by an accelerated version of the extended Fellner-Schall update.
#'
#' \code{areml} is a drop-in replacement for \code{\link{qreml}} with the same interface. It differs in how the outer iteration over the penalty strengths is controlled:
#' \itemize{
#'   \item the update is taken on the log scale with an \strong{adaptive step multiplier} that doubles whenever a small step keeps improving the criterion, which is what makes penalty strengths approaching zero or infinity converge in a sensible number of iterations,
#'   \item \strong{convergence is judged on the restricted likelihood} rather than on the relative change of the penalty strengths, so a penalty strength drifting towards a boundary no longer keeps the iteration alive,
#'   \item the Hessian is obtained by \strong{automatic differentiation} (\code{obj$he}) rather than by finite differencing the gradient, and is repaired with a rank-revealing pivoted Cholesky,
#'   \item only the diagonal blocks of the inverse Hessian that are actually needed are computed, instead of the full inverse.
#' }
#' The step control follows the extended Fellner-Schall implementation in \code{mgcv}.
#'
#' @seealso \code{\link{penalty}} and \code{\link{penalty2}} to compute the penalty inside the likelihood function, and \code{\link{qreml}} for the original algorithm
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
#' @param psname optional name given to the penalty strength parameter in \code{dat}. Defaults to \code{"lambda"}
#' @param alpha smallest factor by which a penalty strength may \strong{decrease} in one outer iteration, a number in [0, 1). Defaults to 0.3.
#'
#' Penalty strengths are free to increase as fast as the update proposes, but cannot collapse faster than this per iteration.
#' Reducing a penalty strength too quickly can push the inner optimisation into a local optimum or a numerically awkward region, which matters more here than in a GAM because the likelihood is user-written.
#' Set to zero to remove the floor entirely. Step length is handled separately, by \code{max_halve} and the adaptive multiplier.
#' @param smoothing optional scaling factor for the final penalty strength parameters. Increasing this beyond one leads to a smoother final model
#' @param maxiter maximum number of outer iterations
#' @param tol convergence tolerance: the iteration stops once the restricted log-likelihood has changed by less than \code{tol} over the last four outer iterations and the step is small. Defaults to 0.1, as in \code{mgcv}
#' @param lsp_max largest value allowed for \code{log(lambda)}. Defaults to 15, as in \code{mgcv}, i.e. penalty strengths saturate at roughly 3.3e6
#' @param step_small size of a step in \code{log(lambda)} below which the step multiplier is allowed to double. Defaults to 0.05, as in \code{mgcv}
#' @param max_halve maximum number of times a step that decreases the restricted likelihood is halved before it is accepted anyway. Defaults to 6.
#'
#' \code{mgcv} never shortens below the full Fellner-Schall step and accepts a worse one instead, which for a user-written likelihood can drift downhill for tens of iterations.
#' @param method optimisation method to be used by \code{\link[stats:optim]{optim}}. Defaults to \code{"BFGS"}
#' @param control list of control parameters for \code{\link[stats:optim]{optim}} to use in the inner optimisation
#' @param spHess logical, if \code{TRUE}, the sparse automatic differentiation Hessian is used for evaluation. The factorisation is dense either way
#' @param joint_unc logical, if \code{TRUE}, joint \code{RTMB} object is returned allowing for joint uncertainty quantification
#' @param saveall logical, if \code{TRUE}, then all model objects from each iteration are saved in the final model object
#'
#' @return model object of class \code{c("aremlModel", "qremlModel")}, so that all methods defined for \code{\link{qreml}} objects apply
#'
#' @export
#'
#' @import RTMB
#'
#' @examples
#' # see ?qreml for a full model-fitting example: areml() is called the same way
#' # mod = areml(pnll, par, dat, random = "betaspline")
areml <- function(pnll, # penalised negative log-likelihood function
                  par, # initial parameter list
                  dat, # initial dat object, currently needs to be called dat!
                  random, # names of parameters in par that are random effects/ penalised
                  map = NULL, # map for fixed effects
                  silent = 1, # print level
                  psname = "lambda", # name given to the psname parameter in dat
                  alpha = 0.3, # smallest factor by which lambda may decrease per iteration
                  smoothing = 1,
                  maxiter = 200, # maximum number of outer iterations
                  tol = 0.1, # convergence tolerance on the restricted log-likelihood
                  lsp_max = 15, # largest allowed log(lambda)
                  step_small = 0.05, # step size below which the multiplier may double
                  max_halve = 6, # how often a worsening step may be halved
                  method = "BFGS", # optimisation method used by optim
                  control = list(), # control list for inner optimisation
                  spHess = FALSE, # evaluate the Hessian sparsely
                  joint_unc = FALSE, # should joint object be returned?
                  saveall = FALSE) # save all intermediate models?
{
  ### input checking
  if(!is.function(pnll)) stop("'pnll' needs to be a function")
  if(!is.list(par)) stop("'par' needs to be a named list")
  if(!is.list(dat)) stop("'dat' needs to be a named list")
  if(!psname %in% names(dat)){
    stop(paste0("'dat' needs to contain a vector called '", psname, "' with initial penalty strengths"))
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

  # setting the argument names because later updated par is returned
  argname_par <- as.character(substitute(par))
  argname_dat <- as.character(substitute(dat))

  n_re <- length(random) # number of distinct random effects
  allmods <- list() # list to save all model objects

  # initialising penalty strength lambda
  lambda <- dat[[psname]]
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
    dat[[psname]] <- DataEval(getLambda, rep(advector(1), 0))

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
  if(is.null(map[[psname]])) map[[psname]] <- factor(seq_along(lambda))
  lambda_map <- map[[psname]]
  if(length(lambda_map) != length(lambda)){
    stop(paste0("Length of map argument for ", psname, " has wrong length."))
  }
  map <- map[names(map) != psname]
  if(length(map) == 0) map <- NULL

  lambda_mapped <- map_lambda(lambda, lambda_map)
  # with every penalty strength fixed there is nothing for the outer iteration to
  # do, and it must be skipped rather than run once: the step is over an empty
  # vector, which is not a convergence failure
  estimate_ps <- length(lambda_mapped) > 0
  if(!estimate_ps) message("No penalty parameters will be estimated as all are fixed.")

  ## creating the RTMB objective function
  if(silent %in% 0:1) message("Creating AD function")
  obj <- MakeADFun(func = f, parameters = par, silent = TRUE, map = map)
  newpar <- obj$par

  ## choosing how the Hessian is evaluated
  # obj$he() is the exact AD Hessian and avoids the finite differencing error in
  # optimHess(), which otherwise sits right on top of the criterion differences
  # the convergence test has to resolve
  if(spHess){
    Tape <- RTMB::GetTape(obj, name = "ADFun")
    if(silent < 2) message("Constructing sparse Hessian")
    spH <- Tape$jacfun(sparse = TRUE)$jacfun(sparse = TRUE)
    rm(Tape)
    hessian_at <- function(p) as.matrix(spH(p))
  } else if(!is.null(obj$he)){
    hessian_at <- function(p) obj$he(p)
  } else {
    if(silent < 2) message("obj$he() unavailable, falling back to finite differences")
    hessian_at <- function(p) stats::optimHess(p, obj$fn, obj$gr)
  }

  ## gradient printing
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

  # finding the indices of the random effects to later index the Hessian
  re_inds <- list()
  for(i in seq_len(n_re)){
    if(is.vector(par[[random[i]]])){
      re_dim <- c(1, length(par[[random[i]]]))
    } else if(is.matrix(par[[random[i]]])){
      re_dim <- dim(par[[random[i]]])
    } else stop(paste0(random[i], " must be a vector or matrix"))

    byrow <- FALSE
    if(is.matrix(S[[i]])){
      if(re_dim[1] == nrow(S[[i]])) byrow <- TRUE
    } else if(is.list(S[[i]])){
      if(re_dim[1] == nrow(S[[i]][[1]])) byrow <- TRUE
    }

    re_inds[[i]] <- matrix(which(names(obj$par) == random[i]), nrow = re_dim[1], ncol = re_dim[2])
    if(byrow) re_inds[[i]] <- t(re_inds[[i]])
  }

  ## how many penalty strengths per random effect: 1 = simple smooth, >1 = tensor product
  n_penalties <- sapply(S, function(x) if(is.matrix(x)) 1 else length(x))
  simple_ind <- which(n_penalties == 1)
  tp_ind <- which(n_penalties > 1)

  re_lengths <- sapply(re_inds, function(x) if(is.vector(x)) 1 else nrow(x))
  lambda_lengths <- n_penalties * re_lengths
  if(length(lambda) != sum(lambda_lengths)){
    stop(paste0("Length of '", psname, "' does not match the number of penalty strength parameters needed"))
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
  Lambda0 <- reshape_lambda(lambda_lengths, lambda)
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

  ## controlling optim
  ctl <- list(maxit = 1000)
  ctl[names(control)] <- control
  if(method == "BFGS") ctl$reltol <- 1e-10
  if(method == "L-BFGS-B") ctl$maxit <- 5000

  ## one complete evaluation at a given log penalty strength:
  ## inner fit, Hessian, factorisation and criterion
  fit_at <- function(lsp_try, start) {
    lambda <<- unmap_lambda(exp(lsp_try), lambda_map, lambda0)
    Lambda <- reshape_lambda(lambda_lengths, lambda)

    # naming the penalty strengths here makes the trace unambiguous: the outer
    # line below reports lambda AFTER the update, so without this the initial fit
    # and the first updated fit look like the same fit done twice
    if(silent == 0){
      cat("\nInner optimisation at", psname, "=", round(exp(lsp_try), 3), "\n")
    }
    counter_env$count <- 0

    # RTMB/TMB caches the objective value at the last parameter vector it was
    # evaluated at. optim's very first evaluation is at 'start', which is exactly
    # where the previous inner fit left off, so without this the first function
    # value is the one belonging to the PREVIOUS lambda. Being the optimum of a
    # smaller penalty it is too low, no trial point can beat it, and optim returns
    # immediately with that stale value: the criterion then looks far better than
    # it is and the step is accepted on false evidence.
    # Evaluating once at a nudged parameter vector invalidates the cache, so the
    # evaluation optim then makes at 'start' is a real one. The value here is
    # discarded, and the nudge is never used as a starting value.
    # qreml() never hit this because optimHess() finite differences obj$fn at many
    # perturbed points and so invalidates the cache as a side effect; obj$he()
    # never touches obj$fn, which is what exposed it.
    invisible(try(obj$fn(start + 1e-8 * (1 + abs(start))), silent = TRUE))

    opt <- stats::optim(start, obj$fn, newgrad, method = method, control = ctl)

    gr <- obj$gr(opt$par) # forces evaluation at opt$par so that report() matches
    if(silent == 0) cat("iter", counter_env$count, "- inner mgc:", round(max(abs(gr)), 5), "\n")

    mod <- obj$report()
    if(silent == 0) cat("evaluating Hessian...\n")
    J <- hessian_at(opt$par)
    J <- (J + t(J)) / 2 # force symmetric
    fac <- pd_chol(J, silent)

    # restricted log-likelihood; the criterion is its negative, so it is minimised
    llk_r <- -opt$value + logdet_Slambda(Lambda) / 2 - fac$logdet / 2

    list(opt = opt, mod = mod, J = J, fac = fac, Lambda = Lambda,
         lsp = lsp_try, llk_r = llk_r, crit = -llk_r, mgc = max(abs(gr)))
  }

  ## the extended Fellner-Schall ratio r, on the mapped scale
  # lambda_new = lambda * r, i.e. the step in log(lambda) is log(r)
  efs_ratio <- function(state) {
    a <- rep(NA_real_, length(lambda0))
    bSb <- rep(NA_real_, length(lambda0))
    l <- 1
    for(i in seq_len(n_re)){
      for(j in seq_len(nrow(re_inds[[i]]))){
        idx <- re_inds[[i]][j, ]
        Jinv <- block_inv(state$fac, idx) # only the block we need

        if(i %in% simple_ind){
          # tr(S_lambda^- S) = rank / lambda for a single penalty matrix
          a[l] <- ranks[i] / state$Lambda[[i]][j] - sum(Jinv * S[[i]])
          bSb[l] <- state$mod$Pen[[i]][j]
          l <- l + 1
        } else {
          n_pen <- length(S[[i]])
          Sinv <- gen_inverse(block_S(i, j, state$Lambda))$inv
          for(pen in seq_len(n_pen)){
            a[l + pen - 1] <- sum(Sinv * S[[i]][[pen]]) - sum(Jinv * S[[i]][[pen]])
          }
          # penalty2() reports the bare b^t S_k b, without lambda_k
          bSb[l:(l + n_pen - 1)] <- state$mod$Pen[[i]][[j]]
          l <- l + n_pen
        }
      }
    }

    # Flooring both keeps the ratio positive, which is what guarantees a positive
    # penalty strength and means no test on the data Hessian is needed for that.
    # The floor has to be far below anything meaningful, though: mgcv uses
    # sqrt(.Machine$double.eps), which here would clamp BOTH quantities for a
    # strongly penalised smooth whose coefficients have been shrunk to nothing,
    # giving r = 1 and a penalty strength frozen for good. The ratio of two very
    # small numbers is still informative, so only positivity is enforced.
    tiny <- .Machine$double.xmin
    a <- pmax(tiny, a)
    bSb <- pmax(tiny, bSb)

    # summing within groups of tied penalty strengths
    a_m <- bSb_m <- numeric(length(lambda_mapped))
    for(m in seq_along(lambda_mapped)){
      ind <- which(lambda_map == levels(lambda_map)[m])
      a_m[m] <- sum(a[ind])
      bSb_m[m] <- sum(bSb[ind])
    }
    r <- a_m / bSb_m
    r[!is.finite(r)] <- 1e6
    # no single step may move a penalty strength by more than a factor 1e6,
    # mirroring the bound mgcv puts on a non-finite ratio
    r <- pmin(pmax(r, 1e-6), 1e6)
    list(r = r, a = a_m, bSb = bSb_m)
  }

  ### updating algorithm
  lsp <- log(lambda_mapped)
  mult <- 1 # step multiplier, persistent across iterations, as in mgcv
  crit_hist <- rep(NA_real_, maxiter)
  llk_hist <- rep(NA_real_, maxiter)

  if(silent < 2) message("Initialising with ", psname, ": ", paste(round(lambda, 3), collapse = " "))

  cur <- fit_at(lsp, newpar)
  converged <- FALSE

  # Best iterate seen so far. Once the multiplier is at its floor the step can no
  # longer be shortened, so a step that worsens the criterion is accepted anyway
  # (this is mgcv's rule, and the iteration often climbs back out of such a dip).
  # Keeping the best point means a run that does not climb back out still returns
  # the best penalty strengths it found rather than wherever it happened to stop.
  best_lsp <- lsp
  best_crit <- cur$crit
  best_iter <- 0

  if(!estimate_ps){
    crit_hist[1] <- cur$crit
    converged <- TRUE
  }

  for(iter in seq_len(if(estimate_ps) maxiter else 0)){

    if(saveall) allmods[[iter]] <- cur$mod

    step <- log(efs_ratio(cur)$r)
    # downward floor: a penalty strength may not fall by more than a factor alpha
    # in one iteration, which protects the inner optimisation early on
    if(alpha > 0) step <- pmax(step, log(alpha))

    lsp1 <- pmin(lsp + mult * step, lsp_max)
    max_step <- max(abs(lsp1 - lsp))

    trial <- fit_at(lsp1, cur$opt$par)

    if(trial$crit <= cur$crit){
      ## improved
      if(max_step < step_small){
        # mgcv's acceleration: the step is small and still paying, so try twice
        # as far and keep the doubling if it pays again
        lsp2 <- pmin(lsp + 2 * mult * step, lsp_max)
        trial2 <- fit_at(lsp2, cur$opt$par)
        if(trial2$crit < trial$crit){
          trial <- trial2
          lsp1 <- lsp2
          mult <- mult * 2
          if(silent == 0) cat("accelerating: step multiplier now", mult, "\n")
        }
      } else if(mult < 1){
        # a shortening from an earlier iteration must not become permanent, or
        # the iteration crawls for the rest of the run: walk the multiplier back
        # towards the full step now that it is paying again
        mult <- min(2 * mult, 1)
      }
    } else {
      ## worsened: shorten the step until it pays. mgcv stops at the full step and
      ## accepts a worse one, which here can drift downhill for tens of iterations
      ## at a time, so the step is genuinely backtracked instead.
      n_halve <- 0
      while(trial$crit > cur$crit && n_halve < max_halve){
        mult <- mult / 2
        n_halve <- n_halve + 1
        if(silent == 0) cat("criterion increased; step multiplier now", signif(mult, 4), "\n")
        lsp1 <- pmin(lsp + mult * step, lsp_max)
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
    llk_hist[iter] <- -cur$opt$value + cur$mod$pen # unpenalised log-likelihood

    if(silent < 2){
      if(silent == 0) cat("\n")
      cat("outer", iter, "-", paste0(psname, ":"), round(exp(lsp), 3), "\n")
      if(silent == 0) cat("restricted llk:", round(cur$llk_r, 5), "- max step:", round(max_step, 5), "\n")
    }

    #### convergence check ####
    # on the criterion, not on the penalty strengths: a penalty strength drifting
    # towards a boundary changes the criterion by nothing, so it no longer keeps
    # the iteration alive. The window of four absorbs numerical wobble.
    if(iter > 3 && max_step < step_small &&
       max(abs(diff(crit_hist[(iter-3):iter]))) < tol){
      converged <- TRUE
    }
    # secondary criterion: the likelihood itself has stopped changing. Unlike in
    # mgcv this also requires the step to be small: the primary tolerance here is
    # tight enough that the secondary would otherwise stop the iteration while
    # the penalty strengths are still moving quickly.
    if(iter > 1 && max_step < step_small &&
       abs(llk_hist[iter] - llk_hist[iter-1]) < 1e-5 * abs(llk_hist[iter])){
      converged <- TRUE
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

  ## final model fit
  # the best penalty strengths found, which is the last iterate whenever the
  # criterion improved to the end
  if(best_iter < iter && silent < 2){
    message("Restricted likelihood did not improve after iteration ", best_iter,
            "; returning the best penalty strengths found")
  }
  lsp <- best_lsp
  lambda <- unmap_lambda(exp(lsp), lambda_map, lambda0) * smoothing
  if(silent < 2){
    if(any(smoothing != 1)) message("Smoothing factor: ", paste(smoothing, collapse = " "))
    message("Final model fit with ", psname, ": ", paste(round(lambda, 3), collapse = " "))
  }

  final_lsp <- log(map_lambda(lambda, lambda_map))
  # cur is already the converged fit at this lambda whenever no smoothing factor
  # was applied and the last iterate was the best one, so refitting would repeat
  # an inner optimisation, a Hessian and a factorisation for nothing
  if(isTRUE(all.equal(final_lsp, cur$lsp))){
    final <- cur
  } else {
    final <- fit_at(final_lsp, cur$opt$par)
  }
  mod <- final$mod
  opt <- final$opt
  Lambda <- final$Lambda
  if(saveall) allmods[[iter + 1]] <- mod

  pllk <- -opt$value # penalised log-likelihood
  llk <- pllk + mod$pen

  #############################################

  mod$obj <- obj
  if(saveall) mod$allmods <- allmods

  names(lambda) <- lambda_names
  mod[[psname]] <- lambda
  mod[[paste0("all_", psname)]] <- Lambda

  parlist <- obj$env$parList(opt$par)
  mod[[argname_par]] <- parlist
  mod[[paste0("relist_", argname_par)]] <- obj$env$parList
  mod[[paste0("map_", psname)]] <- function(lambda) map_lambda(lambda, lambda_map)
  mod$psname <- psname
  mod$parname <- argname_par
  mod[[paste0(argname_par, "_vec")]] <- opt$par
  mod$llk <- llk
  mod$n_fixpar <- length(unlist(par[!(names(par) %in% random)]))

  ## effective degrees of freedom, from the same blocks of the inverse
  # sum_j diag(J^-1 H)_jj over a smooth equals q - sum(J^-1[idx,idx] * S_lambda[idx,idx]),
  # because S_lambda is block diagonal, so no off-block entry is ever needed
  Edfs <- Lambda
  for(i in seq_len(n_re)){
    if(i %in% tp_ind) Edfs[[i]] <- numeric(nrow(re_inds[[i]]))
    for(j in seq_len(nrow(re_inds[[i]]))){
      idx <- re_inds[[i]][j, ]
      Jinv <- block_inv(final$fac, idx)
      Edfs[[i]][j] <- length(idx) - sum(Jinv * block_S(i, j, Lambda))
    }
  }
  mod$df <- mod$n_fixpar + sum(unlist(Edfs))
  mod$edf <- Edfs

  if(!is.null(mod$allprobs)) mod$nobs <- nrow(mod$allprobs)

  mod$Hessian_conditional <- final$J
  mod$llk_restricted <- -crit_hist[seq_len(iter)]
  mod$converged <- converged
  mod$iter <- iter
  mod$best_iter <- best_iter
  mod$hessian_repaired <- final$fac$repaired

  # gradient of the restricted log-likelihood with respect to log(lambda):
  # dV/dlog(lambda) = 0.5 * lambda * (a - b^t S b), which is zero exactly where
  # the multiplicative update has a fixed point (a = b^t S b, i.e. r = 1)
  fr <- efs_ratio(final)
  mod$outer_grad <- 0.5 * map_lambda(lambda, lambda_map) * (fr$a - fr$bSb)
  names(mod$outer_grad) <- levels(lambda_map)

  # removing elements only reported for the update
  mod <- mod[!names(mod) %in% c("Pen", "pen", "S")]

  if(length(tp_ind) == 0 && joint_unc){
    ### constructing joint object
    parlist$loglambda <- log(mod[[psname]])
    logdetS <- numeric(length(S))
    for(i in seq_along(S)) logdetS[i] <- gdeterminant(S[[i]])

    jnll <- function(par) {
      environment(pnll) = environment()
      "[<-" <- ADoverload("[<-")
      "c" <- ADoverload("c")
      "diag<-" <- ADoverload("diag<-")

      dat[[psname]] <- exp(par$loglambda)
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

  class(mod) <- c("aremlModel", "qremlModel")
  mod
}
