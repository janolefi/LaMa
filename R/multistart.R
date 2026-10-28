#' Optimise an \code{RTMB} object from several starting values
#'
#' Latent Markov models often have several local optima of the likelihood. \code{multistart()} optimises
#' \code{obj} from several starting values: the first is \code{obj$par} (the initial parameters), the others
#' are random perturbations of it on the working scale. Afterwards, \code{\link{lama}(obj)} builds the model
#' from the best fit, as \code{RTMB} keeps the parameters with the lowest negative log-likelihood it evaluated.
#'
#' @details
#' \strong{Stopping.} Two runs reach the same optimum if their negative log-likelihoods differ by less than
#' \code{tol} (relative), so permuted states (label switching) count as the same optimum.
#' Only converged runs (convergence code 0 and maximum gradient component below 0.01) count for stopping.
#' In any case, at most \code{nstarts} starts are run.
#' \itemize{
#'   \item \code{stop = "bayes"} uses the Bayesian stopping rule of Boender and Rinnooy Kan (1987):
#'   after \eqn{k} runs with \eqn{w} distinct optima, the expected number of optima is \eqn{w(k-1)/(k-w-2)},
#'   and \code{multistart()} stops when this is below \eqn{w + 0.5}, i.e. probably no optimum is missed.
#'   The rule assumes that the starts are independent draws from one distribution, so only the random starts count,
#'   not the first run from \code{obj$par}, and it is not used for given \code{starts}.
#'   It needs at least 8 converged random starts.
#'   \item \code{stop = "hits"} stops once the best optimum was reached \code{hits} times.
#' }
#'
#' \strong{Caution.} Stopping does not mean that the global optimum was found. "No optimum missed" only refers to the
#' distribution of the starts: optima far away from it, or with a small basin of attraction, can be missed.
#' For example, with a very small \code{sd}, all starts converge to the same optimum and the rule stops quickly,
#' although other starts might find a better one. So choose \code{sd} wide enough and still check the fit.
#'
#' \strong{Scatter.}
#' \itemize{
#'   \item With \code{scale = "fixed"}, the starts are \code{obj$par} plus normal noise with standard deviation \code{sd}.
#'   As most working parameters are on log or logit scale, \code{sd = 1} is a relative change, so one value usually fits all parameters.
#'   Parameters on identity scale (e.g. means of normal distributions) may need their own \code{sd},
#'   given as a named list like \code{list(eta = 1, mu = 5)}.
#'   \item With \code{scale = "unit_info"}, the starts are drawn around the first fit from a normal distribution with
#'   covariance \eqn{n H^{-1}}, where \eqn{H} is the Hessian of the negative log-likelihood and \eqn{n} the number of observations.
#'   This is the information of a single observation (unit information), which does not shrink with \eqn{n}
#'   and takes correlations of the parameters into account. Standard deviations are capped at \code{sd_max},
#'   for weakly identified parameters. If the Hessian is not positive definite (or \code{allprobs} is not reported,
#'   so \eqn{n} is unknown), the fixed \code{sd} is used instead.
#' }
#'
#' @param obj \code{RTMB} object from \code{MakeADFun()}
#' @param nstarts maximum number of starts, including \code{obj$par}
#' @param sd standard deviation of the perturbations on the working scale: a number, a vector of the length of \code{obj$par},
#' or a named list with one value per parameter (names as in the parameter list)
#' @param starts optional matrix of starting values, one per row, with columns as \code{obj$par}; then \code{nstarts} and \code{sd} are ignored
#' @param scale \code{"fixed"} for perturbations of \code{obj$par} with \code{sd}, or \code{"unit_info"} for perturbations of the first fit
#' scaled by the unit information, see \sQuote{Details}
#' @param stop stopping rule: \code{"bayes"} (not used for given \code{starts}), \code{"hits"} or \code{"none"}, see \sQuote{Details}
#' @param hits number of times the best optimum needs to be reached for \code{stop = "hits"}
#' @param sd_max largest standard deviation for \code{scale = "unit_info"}
#' @param tol relative tolerance for two runs to reach the same optimum
#' @param seed optional seed for the random number generator
#' @param optimizer optimisation function called as \code{optimizer(par, fn, gr, ...)}, returning \code{par} and
#' \code{objective} (like \code{\link[stats:nlminb]{nlminb}}) or \code{value} (like \code{\link[stats:optim]{optim}})
#' @param silent logical; if \code{TRUE}, nothing is printed
#' @param ... further arguments passed to \code{optimizer}
#'
#' @return output of \code{optimizer} for the best run, with two more elements: \code{runs}, a data frame
#' with the negative log-likelihood \code{nll}, \code{convergence} code and maximum gradient component \code{mgc} of each run,
#' and \code{starts}, the matrix of starting values used.
#'
#' @references Boender, C. G. E., & Rinnooy Kan, A. H. G. (1987). Bayesian stopping rules for multistart global
#' optimization methods. \emph{Mathematical Programming}, 37(1), 59-80.
#' @export
#'
#' @examples
#' nll <- function(par) {
#'   getAll(par, dat)
#'   Gamma <- tpm(eta)
#'   mu <- exp(log_mu); REPORT(mu)
#'   allprobs <- cbind(dgamma2(step, mu[1], exp(log_sigma[1])),
#'                    dgamma2(step, mu[2], exp(log_sigma[2])))
#'   -forward(stationary(Gamma), Gamma, allprobs)
#' }
#' dat <- list(step = trex$step[1:500])
#' par <- list(eta = c(-2, -2), log_mu = log(c(0.3, 2.5)), log_sigma = log(c(0.2, 1.5)))
#' obj <- MakeADFun(nll, par, silent = TRUE)
#'
#' opt <- multistart(obj, seed = 1)
#' opt$runs
#' mod <- lama(obj) # best fit
#'
#' # starts scaled by the first fit
#' opt <- multistart(obj, scale = "unit_info", seed = 1)
multistart <- function(obj, nstarts = 20, sd = 1, starts = NULL, scale = c("fixed", "unit_info"),
                       stop = c("bayes", "hits", "none"), hits = 3, sd_max = 2, tol = 1e-6,
                       seed = NULL, optimizer = stats::nlminb, silent = FALSE, ...) {
  obj <- get_obj(obj)
  scale <- match.arg(scale)
  stop <- match.arg(stop)
  p0 <- obj$par
  if(is.list(sd)) sd <- unlist(sd)[names(p0)] # one value per parameter block, as names repeat in obj$par
  if(anyNA(sd)) stop("sd needs a value for each parameter in obj$par.")
  if(!is.null(starts)) {
    starts <- as.matrix(starts)
    if(ncol(starts) != length(p0)) stop("starts needs one column per element of obj$par.")
    nstarts <- nrow(starts)
  }

  # keep the random number stream of the session, as in simulate()
  if(!is.null(seed)) {
    if(!exists(".Random.seed", envir = globalenv())) stats::runif(1)
    old <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old, envir = globalenv()), add = TRUE)
    set.seed(seed)
  }

  # starts around center with covariance Sigma (diagonal for scale = "fixed")
  center <- p0
  Sigma <- diag(rep_len(sd^2, length(p0)), length(p0))
  draw <- function() MASS::mvrnorm(1, center, Sigma)

  # one run, failures give NA, so that the others continue
  fit <- function(p) {
    r <- tryCatch(optimizer(p, obj$fn, obj$gr, ...), error = function(e) NULL)
    if(is.null(r)) return(list(opt = NULL, nll = NA, convergence = NA, mgc = NA))
    nll <- if(!is.null(r$objective)) r$objective else r$value
    list(opt = r, nll = nll, convergence = r$convergence, mgc = max(abs(obj$gr(r$par))))
  }

  opts <- list()
  P <- matrix(NA, 0, length(p0), dimnames = list(NULL, names(p0)))
  runs <- data.frame(nll = numeric(0), convergence = integer(0), mgc = numeric(0))
  stopped <- "nstarts reached"
  for(i in seq_len(nstarts)) {
    p <- if(!is.null(starts)) starts[i, ] else if(i == 1) p0 else draw()
    r <- fit(p)
    P <- rbind(P, p)
    opts[[i]] <- r$opt
    runs[i, ] <- r[c("nll", "convergence", "mgc")]
    if(!silent) cat("start", i, "- nll:", round(r$nll, 4), "\n")

    # after the first fit: scatter around it, scaled by the unit information
    if(i == 1 && scale == "unit_info" && is.null(starts) && !is.null(r$opt)) {
      Sigma_ui <- unit_info(obj, r$opt$par, sd_max)
      if(is.null(Sigma_ui)) {
        if(!silent) message("Unit information not available, using sd instead")
      } else Sigma <- Sigma_ui
      center <- r$opt$par
    }

    # stopping rules on the converged runs
    nll <- runs$nll[converged(runs)]
    if(length(nll) && stop != "none") {
      best <- min(nll)
      if(stop == "hits" && sum(nll - best <= tol * max(1, abs(best))) >= hits) {
        stopped <- "best optimum reached often enough"; break
      }
      # Bayesian rule only on random draws: not run 1 (obj$par) or given starts
      nll_r <- runs$nll[setdiff(converged(runs), 1)]
      k <- length(nll_r)
      w <- if(k) n_optima(nll_r, tol) else 0
      if(stop == "bayes" && is.null(starts) && k > w + 2 && w * (k - 1) / (k - w - 2) < w + 0.5) {
        stopped <- "Bayesian stopping rule"; break
      }
    }
  }

  ibest <- which.min(runs$nll)
  if(!length(ibest)) stop("All runs failed.")
  out <- opts[[ibest]]
  out$runs <- runs
  out$starts <- unname(P)

  if(!silent) {
    nll <- runs$nll[converged(runs)]
    best <- runs$nll[ibest]
    cat(nrow(runs), "starts,", length(nll), "converged,", if(length(nll)) n_optima(nll, tol) else 0, "distinct optima,",
        "best reached", sum(nll - best <= tol * max(1, abs(best))), "times; stopped:", stopped, "\n")
  }
  out
}

# converged runs: convergence code 0 and small gradient, as optimisers can stop on flat regions
converged <- function(runs) which(runs$convergence == 0 & runs$mgc < 0.01)

# number of distinct optima: sorted nll values further apart than tol (relative)
n_optima <- function(nll, tol) {
  nll <- sort(nll)
  sum(c(TRUE, diff(nll) > tol * pmax(1, abs(nll[-1]))))
}

# covariance of one observation n * H^-1 at p, standard deviations capped at sd_max;
# NULL if H is not positive definite or n unknown
unit_info <- function(obj, p, sd_max) {
  n <- count_nobs(obj$report(obj$env$last.par.best)$allprobs)
  if(is.null(n)) return(NULL) # allprobs not reported
  H <- tryCatch(obj$he(p), error = function(e) stats::optimHess(p, obj$fn, obj$gr)) # obj$he() has no random effects
  R <- tryCatch(chol(H), error = function(e) NULL)
  if(is.null(R)) return(NULL)
  Sigma <- n * chol2inv(R)
  s <- sqrt(diag(Sigma))
  Sigma * tcrossprod(pmin(s, sd_max) / s) # same correlations, capped standard deviations
}
