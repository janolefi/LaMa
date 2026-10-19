# Number of observations in all forward() calls of nll (counting mode),
# found by one evaluation of nll at par without AD
count_obs <- function(nll, par){
  .lama_env$count <- 0
  on.exit(.lama_env$count <- NULL)
  nll(par)
  .lama_env$count
}


#' Per-observation score contributions
#'
#' Computes the score contribution of each observation,
#' \deqn{u_t = \frac{\partial}{\partial \theta} \log p(x_t \mid x_1, \dots, x_{t-1}),}
#' for a model whose negative log-likelihood function calls \code{\link{forward}}.
#' Their column sums are the gradient of the log-likelihood. They are the building block for
#' robust (sandwich) standard errors, e.g. cluster-robust ones via \code{rowsum(U, ID)}.
#'
#' @details
#' Each log-likelihood contribution in \code{\link{forward}} gets a weight \eqn{w_t}.
#' The gradient of the weighted log-likelihood is linear in the weights, so its derivative
#' w.r.t. \eqn{w_t} is \eqn{u_t}. This needs one tape of the gradient and one reverse sweep
#' per parameter. The weights are handed to \code{\link{forward}} internally, so \code{nll}
#' does not need to change.
#'
#' The number of observations is found by one evaluation of \code{nll} without automatic differentiation.
#' If \code{nll} calls \code{\link{forward}} several times, e.g. once per track, the rows of the
#' result follow the order of these calls. Only \code{\link{forward}} (and its deprecated
#' \code{_g} and \code{_p} versions) is supported, and models with random effects are not.
#'
#' @param nll negative log-likelihood function of a named parameter list, as passed to \code{RTMB::MakeADFun}
#' @param par named list of parameters, as passed to \code{RTMB::MakeADFun}
#' @param est optional vector of parameter values in the order of \code{unlist(par)}, e.g. \code{opt$par}.
#' Defaults to the values in \code{par}.
#'
#' @return matrix of score contributions of dimension \code{c(nObs, length(unlist(par)))}
#' @export
#'
#' @examples
#' nll <- function(par) {
#'   getAll(par)
#'   Gamma <- tpm(eta)
#'   delta <- stationary(Gamma)
#'   allprobs <- cbind(dgamma2(step, mu[1], sigma[1]), dgamma2(step, mu[2], sigma[2]))
#'   -forward(delta, Gamma, allprobs)
#' }
#' step <- trex$step[1:200]
#' par <- list(eta = c(-2, -2), mu = c(0.3, 2.5), sigma = c(0.2, 1.5))
#' obj <- RTMB::MakeADFun(nll, par, silent = TRUE)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' U <- scores(nll, par, est = opt$par)
#' colSums(U) # ~ 0 at the MLE
scores <- function(nll, par, est = NULL){
  theta <- unlist(par)
  p <- length(theta)
  if(!is.null(est)) {
    if(length(est) != p) stop("est must have the same length as unlist(par).")
    theta[] <- est
  }
  par <- utils::relist(unname(theta), par) # parameter list at est

  nObs <- count_obs(nll, par)
  if(nObs == 0) stop("nll contains no call of forward().")

  # nll with weights w for the log-likelihood contributions in forward()
  used <- NA
  f <- function(pw) {
    .lama_env$weights <- pw$weights_
    .lama_env$pos <- 0
    on.exit(.lama_env$weights <- NULL)
    pw$weights_ <- NULL
    val <- nll(pw)
    used <<- .lama_env$pos
    val
  }
  obj_w <- RTMB::MakeADFun(f, c(par, list(weights_ = rep(1, nObs))), silent = TRUE)
  if(used != nObs) stop("forward() used ", used, " of ", nObs, " observations while taping.")

  # gradient w.r.t. (theta, w), restricted to theta as a function of w (linear in w)
  G <- RTMB::GetTape(obj_w)$jacfun()
  Gw <- RTMB::MakeTape(function(w) G(c(RTMB::AD(unname(theta)), w))[1:p], rep(1, nObs))

  # d/dw_t of the gradient of nll is -u_t
  U <- -t(Gw$jacobian(rep(1, nObs)))
  colnames(U) <- names(theta)
  U
}
