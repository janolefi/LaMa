# Number of observations in all forward() calls of nll (counting mode),
# found by one evaluation of nll at par without AD
count_obs <- function(nll, par){
  .lama_env$count <- 0
  on.exit(.lama_env$count <- NULL)
  nll(par)
  .lama_env$count
}


# Frame of the MakeADFun() call that created obj, with the user's function (func), map and random.
# Not part of RTMB's interface, hence checked here.
obj_frame <- function(obj){
  fr <- environment(obj$report)
  if(!is.function(fr$func)) stop("Could not find the likelihood function in obj (has RTMB changed?).")
  fr
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
#' per parameter. The weights are handed to \code{\link{forward}} internally, so the
#' likelihood function does not need to change.
#'
#' The likelihood function is taken from \code{obj} and evaluated again, with the data as they are
#' when \code{scores()} is called. Parameters fixed via \code{map} are respected, and the columns
#' correspond to \code{obj$par}, like \code{obj$gr()} and \code{obj$he()}.
#'
#' If the likelihood function calls \code{\link{forward}} several times, e.g. once per track, the rows of the
#' result follow the order of these calls. Only \code{\link{forward}} (and its deprecated
#' \code{_g} and \code{_p} versions) is supported, and models with random effects are not.
#'
#' @param obj optimised \code{RTMB} object from \code{\link[RTMB]{MakeADFun}}, or a model object from \code{\link{lama}} or \code{\link{qreml}}
#' @param par optional parameter vector like \code{obj$par}, e.g. \code{opt$par}.
#' Defaults to the estimates (for \code{\link{qreml}} models the final fit).
#'
#' @return matrix of score contributions of dimension \code{c(nObs, length(obj$par))}
#' @export
#'
#' @examples
#' nll <- function(par) {
#'   getAll(par)
#'   Gamma <- tpm(eta)
#'   delta <- stationary(Gamma)
#'   allprobs <- cbind(dgamma2(step, mu[1], sigma[1]), 
#'                     dgamma2(step, mu[2], sigma[2]))
#'   -forward(delta, Gamma, allprobs)
#' }
#' step <- trex$step[1:200]
#' par <- list(eta = c(-2, -2), mu = c(0.3, 2.5), sigma = c(0.2, 1.5))
#' obj <- MakeADFun(nll, par, silent = TRUE)
#' opt <- nlminb(obj$par, obj$fn, obj$gr)
#' mod <- lama(obj)
#' U <- scores(mod)
#' colSums(U) # ~ 0 at the MLE
scores <- function(obj, par = NULL){
  tape_scores(weighted_tape(obj, par))
}

# score contributions from a tape of weighted_tape(), so that callers needing the tape
# for more (e.g. weighted_hessian()) only tape once
tape_scores <- function(tp){
  w1 <- rep(1, tp$nObs)
  # gradient w.r.t. theta as a function of w (linear in w); d/dw_t of the gradient of the
  # negative log-likelihood is -u_t
  Gw <- RTMB::MakeTape(function(w) tp$G(c(RTMB::AD(unname(tp$theta)), w))[seq_len(tp$p)], w1)
  U <- -t(Gw$jacobian(w1))
  colnames(U) <- tp$names
  U
}


# Taped gradient of the weighted negative log-likelihood -sum_t w_t log p(x_t | x_1, ..., x_t-1),
# with weights w_t on the log-likelihood contributions in forward(). Shared by scores() and
# weighted_hessian(). The likelihood function is taken from obj and evaluated again, with the data
# as they are now; map is respected, random effects are not supported.
# Returns a list with
#   theta: parameters (like obj$par) at which the tape was made, by default the estimates
#   G:     function of c(theta, w) returning the gradient w.r.t. (theta, w); its first p entries
#          are the gradient of the weighted negative log-likelihood w.r.t. theta
#   p, nObs, names: number of parameters, number of observations, names(obj$par)
weighted_tape <- function(obj, par = NULL){
  if(is.null(par)) par <- est_par(obj)
  obj <- get_obj(obj)
  fr <- obj_frame(obj)
  if(length(fr$random)) stop("Score contributions are not available for models with random effects.")
  p <- length(obj$par)
  if(length(par) != p) stop("par must have the same length as obj$par.")

  parlist <- obj$env$parList(par) # all parameters at par, respecting map
  attr(parlist, "check.passed") <- NULL

  nObs <- count_obs(fr$func, parlist)
  if(nObs == 0) stop("The likelihood function contains no call of forward().")

  # likelihood function with weights w for the log-likelihood contributions in forward()
  used <- NA
  f <- function(pw) {
    .lama_env$weights <- pw$weights_
    .lama_env$pos <- 0
    on.exit(.lama_env$weights <- NULL)
    pw$weights_ <- NULL
    val <- fr$func(pw)
    used <<- .lama_env$pos
    val
  }
  obj_w <- RTMB::MakeADFun(f, c(parlist, list(weights_ = rep(1, nObs))), map = fr$map, silent = TRUE)
  if(used != nObs) stop("forward() used ", used, " of ", nObs, " observations while taping.")

  list(theta = obj_w$par[seq_len(p)], G = RTMB::GetTape(obj_w)$jacfun(), p = p, nObs = nObs,
       names = names(obj$par))
}


# Hessian H(w) = -sum_t w_t d^2 log p(x_t | x_1, ..., x_t-1) / d theta d theta' of the weighted negative
# log-likelihood, the full p x p matrix in the order of obj$par: the Jacobian of the taped gradient
# w.r.t. theta at fixed weights (about p gradient sweeps). Linear in w, so w = 1 gives obj$he(),
# w = (trackID == k) the contribution H_k of track k and w = (trackID != k) the leave-one-track-out
# Hessian H - H_k. w is just a numeric (or logical) vector, so any type of track identifier works.
weighted_hessian <- function(tp, w = rep(1, tp$nObs), theta = tp$theta){
  if(length(w) != tp$nObs) stop("w must have length ", tp$nObs, " (number of observations in forward()).")
  w <- as.numeric(w)
  th <- unname(theta)
  Hw <- RTMB::MakeTape(function(x) tp$G(c(x, RTMB::AD(w)))[seq_len(tp$p)], th)
  H <- Hw$jacobian(th)
  H <- (H + t(H)) / 2
  dimnames(H) <- list(tp$names, tp$names)
  H
}
