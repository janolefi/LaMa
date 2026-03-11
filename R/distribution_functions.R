#' Von Mises distribution
#' 
#' Density, distribution function and random generation for the von Mises distribution.
#' 
#' @details
#' This implementation of \code{dvm} allows for automatic differentiation with \code{RTMB}. 
#' \code{rvm} and \code{pvm} are simply wrappers of the corresponding functions from \code{circular}.
#'
#' @param x,q vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param kappa non-negative numeric value for the concentration parameter of the distribution.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param tol the precision in evaluating the distribution function
#' @param from value from which the integration for CDF starts. If \code{NULL}, is set to \code{mu - pi}.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval [-pi, pi].
#'
#' @return \code{dvm} gives the density, \code{pvm} gives the distribution function, and \code{rvm} generates random deviates.
#'
#' @examples 
#' set.seed(1)
#' x = rvm(10, 0, 1)
#' d = dvm(x, 0, 1)
#' p = pvm(x, 0, 1)
#' @name vm
NULL

#' @rdname vm
#' @export
#' @importFrom RTMBdist dvm
dvm <- function(x, mu = 0, kappa = 1, log = FALSE) {
  RTMBdist::dvm(x, mu = mu, kappa = kappa, log = log)
}

#' @rdname vm
#' @export
#' @importFrom RTMBdist pvm
pvm <- function(q, mu = 0, kappa = 1, from = NULL, tol = 1e-20) {
  RTMBdist::pvm(q, mu=mu, kappa=kappa, from=from, tol=tol)
}

#' @rdname vm
#' @export
#' @importFrom RTMBdist rvm
rvm <- function(n, mu = 0, kappa = 1, wrap = TRUE) {
  RTMBdist::rvm(n, mu=mu, kappa=kappa, wrap=wrap)
}


#' Wrapped Cauchy distribution
#' 
#' Density and random generation for the wrapped Cauchy distribution.
#' 
#' @details
#' This implementation of \code{dwrpcauchy} allows for automatic differentiation with \code{RTMB}. 
#' \code{rwrpcauchy} is simply a wrapper for \code{rwrappedcauchy} imported from \code{circular}.
#'
#' @param x vector of angles measured in radians at which to evaluate the density function.
#' @param mu mean direction of the distribution measured in radians.
#' @param rho concentration parameter of the distribution, must be in the interval from 0 to 1.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param wrap logical; if \code{TRUE}, generated angles are wrapped to the interval [-pi, pi].
#'
#' @return \code{dwrpcauchy} gives the density and \code{rwrpcauchy} generates random deviates.
#'
#' @examples 
#' set.seed(1)
#' x = rwrpcauchy(10, 0, 0.3)
#' d = dwrpcauchy(x, 0, 0.3)
#' @name wrpcauchy
NULL

#' @rdname wrpcauchy
#' @export
#' @importFrom RTMBdist dwrpcauchy
dwrpcauchy <- function(x, mu = 0, rho, log = FALSE) {
  RTMBdist::dwrpcauchy(x, mu=mu, rho=rho, log=log)
}

#' @rdname wrpcauchy
#' @export
#' @importFrom RTMBdist rwrpcauchy
rwrpcauchy <- function(n, mu = 0, rho, wrap = TRUE) {
  RTMBdist::rwrpcauchy(n, mu=mu, rho=rho, wrap=wrap)
}


#' Reparametrised gamma distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the gamma distribution reparametrised in terms of mean and standard deviation.
#' 
#' @details
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param mean mean parameter, must be positive scalar.
#' @param sd standard deviation parameter, must be positive scalar.
#' @param log,log.p logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param lower.tail logical; if \code{TRUE}, probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @return
#' \code{dgamma2} gives the density, \code{pgamma2} gives the distribution function, \code{qgamma2} gives the quantile function, and \code{rgamma2} generates random deviates.
#'
#' @examples
#' x = rgamma2(1)
#' d = dgamma2(x)
#' p = pgamma2(x)
#' q = qgamma2(p)
#' @name gamma2
NULL

#' @rdname gamma2
#' @export
#' @importFrom RTMBdist dgamma2
dgamma2 = function(x, mean = 1, sd = 1, log = FALSE) {
  RTMBdist::dgamma2(x, mean=mean, sd=sd, log=log)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMBdist pgamma2
pgamma2 = function(q, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  RTMBdist::pgamma2(q, mean=mean, sd=sd, lower.tail=lower.tail, log.p=log.p)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMBdist qgamma2
qgamma2 = function(p, mean = 1, sd = 1, lower.tail = TRUE, log.p = FALSE) {
  RTMBdist::qgamma2(p, mean=mean, sd=sd, lower.tail = lower.tail, log.p = log.p)
}

#' @rdname gamma2
#' @export
#' @importFrom RTMBdist rgamma2
rgamma2 = function(n, mean = 1, sd = 1) {
  RTMBdist::rgamma2(n, mean=mean, sd=sd)
}



#' Skew normal distribution
#' 
#' Density, distribution function, quantile function and random generation for
#' the skew normal distribution.
#' 
#' @details
#' This implementation of \code{dskewnorm} allows for automatic differentiation with \code{RTMB} while the other functions are imported from the \code{sn} package.
#'
#' @param x,q vector of quantiles
#' @param p vector of probabilities
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param xi location parameter
#' @param omega scale parameter, must be positive.
#' @param alpha skewness parameter, +/- \code{Inf} is allowed.
#' @param log logical; if \code{TRUE}, probabilities/ densities \eqn{p} are returned as \eqn{\log(p)}.
#' @param ... additional parameters to be passed to the \code{sn} package functions for \code{pskewnorm} and \code{qskewnorm}.
#'
#' @return
#' \code{dskewnorm} gives the density, \code{pskewnorm} gives the distribution function, \code{qskewnorm} gives the quantile function, and \code{rskewnorm} generates random deviates.
#'
#' @examples
#' x = rskewnorm(1)
#' d = dskewnorm(x)
#' p = pskewnorm(x)
#' q = qskewnorm(p)
#' @name skewnorm
NULL

#' @rdname skewnorm
#' @export
#' @importFrom RTMBdist dskewnorm
dskewnorm <- function(x, xi = 0, omega = 1, alpha = 0, log = FALSE) {
  RTMBdist::dskewnorm(x, xi=xi, omega=omega, alpha=alpha, log=log)
}

#' @rdname skewnorm
#' @export
#' @importFrom RTMBdist pskewnorm
pskewnorm <- function(q, xi = 0, omega = 1, alpha = 0, ...) {
  RTMBdist::pskewnorm(q, xi=xi, omega=omega, alpha=alpha, log=log, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom RTMBdist qskewnorm
qskewnorm <- function(p, xi = 0, omega = 1, alpha = 0, ...) {
  RTMBdist::qskewnorm(p, xi=xi, omega=omega, alpha=alpha, ...)
}

#' @rdname skewnorm
#' @export
#' @importFrom RTMBdist rskewnorm
rskewnorm <- function(n, xi = 0, omega = 1, alpha = 0) {
  RTMBdist::rskewnorm(n, xi=xi, omega=omega, alpha=alpha)
}


#' Zero-inflated density constructer
#' 
#' @description 
#' Constructs a zero-inflated density function from a given probability density function
#' 
#' @details
#' The definition of zero-inflation is different for discrete and continuous distributions.
#' For discrete distributions with p.m.f. \eqn{f} and zero-inflation probability \eqn{p}, we have
#' \deqn{\Pr(X = 0) = p + (1 - p) \cdot f(0),} and
#' \deqn{\Pr(X = x) = (1 - p) \cdot f(x), \quad x > 0.}
#' 
#' For continuous distributions with p.d.f. \eqn{f}, we have
#' \deqn{f_{\text{zinfl}}(x) = p \cdot \delta_0(x) + (1 - p) \cdot f(x),}
#' where \eqn{\delta_0} is the Dirac delta function at zero.
#' 
#' @param dist either a probability density function or a probability mass function
#' @param discrete logical; if \code{TRUE}, the density for \code{x = 0} will be \code{zeroprob + (1-zeroprob) * dist(0, ...)}. Otherwise it will just be \code{zeroprob}.
#' In standard cases, this will be determined automatically. For non-standard cases, set this to \code{TRUE} or \code{FALSE} depending on the type of \code{dist}. See details.
#'
#' @returns zero-inflated density function with first argument \code{x}, second argument \code{zeroprob}, and additional arguments \code{...} that will be passed to \code{dist}.
#' @export
#'
#' @examples
#' dzinorm <- zero_inflate(dnorm)
#' dzinorm(c(NA, 0, 2), 0.5, mean = 1, sd = 1)
#' 
#' zipois <- zero_inflate(dpois)
#' zipois(c(NA, 0, 1), 0.5, 1)
zero_inflate <- function(dist, discrete = NULL) {
  dist_name <- deparse(substitute(dist))
  
  # Auto-detect whether the function is discrete
  discrete_dists <- c("dpois", "dbinom", "dnbinom", "dgeom", "dmultinom", "dhyper")
  if (is.null(discrete)) {
    discrete <- dist_name %in% discrete_dists
  }
  
  # Extract function argument names, excluding `x`
  dist_args <- setdiff(names(formals(dist)), "x")
  
  # Define core function logic for both discrete and continuous cases
  core_function <- function(x, zeroprob, dist, discrete, dist_args, ...) {
    "[<-" <- RTMB::ADoverload("[<-")  # Ensure correct assignment behavior
    
    args <- list(...)
    
    # Initialize output vector
    out <- numeric(length(x))
    
    # Handle NA values properly
    is_na <- is.na(x)
    is_zero <- x == 0 & !is_na
    is_nonzero <- x != 0 & !is_na  # Exclude NA and 0
    
    # Common logic: assign probability at zero
    out[is_zero] <- zeroprob  
    
    if (discrete) {
      # Logic for discrete distributions: Handle x = 0 case with dist(0)
      out[is_zero] <- out[is_zero] + (1 - zeroprob) * do.call(dist, c(list(0), args))
    }
    
    # Logic for non-zero values
    out[is_nonzero] <- (1 - zeroprob) * do.call(dist, c(list(x[is_nonzero]), args))
    
    # Ensure NAs remain as NA
    out[is_na] <- NA
    
    return(out)
  }
  
  # Return function with core logic encapsulated
  function(x, zeroprob, ...) {
    core_function(x, zeroprob, dist, discrete, dist_args, ...)
  }
}


#' Reparametrised multivariate Gaussian distribution
#'
#' Density function of the multivariate Gaussian distribution reparametrised in terms of its precision matrix (inverse variance).
#' This implementation is particularly useful for defining the \strong{joint log-likelihood} with penalised splines or i.i.d. random effects that have a multivariate Gaussian distribution with fixed precision/ penalty matrix \eqn{\lambda S}.
#' As \eqn{S} is fixed and only scaled by \eqn{\lambda}, it is more efficient to precompute the determinant of \eqn{S} (for the normalisation constant) and only scale the quadratic form by \eqn{\lambda}
#' when multiple spline parameters/ random effects with different \eqn{\lambda}'s but the same penalty matrix \eqn{S} are evaluated.
#'
#' This implementation allows for automatic differentiation with \code{RTMB}.
#'
#' @param x density evaluation point, either a vector or a matrix
#' @param mu mean parameter. Either scalar or vector
#' @param S unscaled precision matrix
#' @param lambda precision scaling parameter
#' 
#' Can be a vector if \code{x} is a matrix. Then each row of \code{x} is evaluated with the corresponding \code{lambda}.
#' This is benefitial from an efficiency perspective because the determinant of \code{S} is only computed once.
#' @param logdetS Optional precomputed log determinant of the precision matrix \code{S}. If the precision matrix does not depend on parameters, it can be precomputed and passed to the function.
#' @param log logical; if \code{TRUE}, densities are returned on the log scale.
#'
#' @return vector of density values
#' @export
#' @import RTMB
#'
#' @examples
#' x = matrix(runif(30), nrow = 3)
#'
#' # iid random effects
#' S = diag(10)
#' sigma = c(1, 2, 3) # random effect standard deviations
#' lambda = 1 / sigma^2
#' d = dgmrf2(x, 0, S, lambda)
#'
#' # P-splines
#' L = diff(diag(10), diff = 2) # second-order difference matrix
#' S = t(L) %*% L
#' lambda = c(1,2,3)
#' d = dgmrf2(x, 0, S, lambda, log = TRUE)
dgmrf2 = function(x, 
                  mu = 0, 
                  S, 
                  lambda, 
                  logdetS = NULL,
                  log = FALSE) {
  # if logdet not specified, compute generalised determinant
  if(is.null(logdetS)){
    logdetS = gdeterminant(S)
  }
  k = nrow(S)
  
  x_centered = x - mu # center data
  
  if(is.null(dim(x_centered))){ # vector
    quadform = lambda * crossprod(x_centered, S %*% x_centered)
  } else{ # matrix
    y = S %*% t(x_centered)
    quadform = lambda * rowSums(x_centered * t(y))
  }
  
  logdens = as.numeric(0.5 * (-k * log(2*pi) + k * log(lambda) + logdetS - quadform))
  
  if(log){
    return(logdens)
  } else{
    return(exp(logdens))
  }
}



#' State dwell-time distributions of periodically inhomogeneous Markov chains
#' 
#' @description
#' Computes the dwell-time distribution of a periodically inhomogeneous Markov chain for a given transition probability matrix.
#' 
#' @details
#' For Markov chains whose transition probabilities vary only periodically, which is achieved for example by
#' expressing the transition probability matrix as a periodic function of the time of day using \code{\link{tpm_p}} or \code{\link{cosinor}}, the probability distribution of time spent in a state can be computed analytically.
#' This function computes said distribution, either for a specific time point (conditioning on transitioning into the state at that time point) or for the overall distribution (conditioning on transitioning into the state at any time point).
#'
#' @references Koslik, J. O., Feldmann, C. C., Mews, S., Michels, R., & Langrock, R. (2023). Inference on the state process of periodically inhomogeneous hidden Markov models for animal behavior. arXiv preprint arXiv:2312.14583.
#'
#' @param x vector of (non-negative) dwell times to compute the dwell-time distribution for
#' @param Gamma array of \code{L} unique transition probability matrices of a periodically inhomogeneous Markov chain, with dimensions \code{c(N,N,L)}, where \code{N} is the number of states and \code{L} is the cycle length
#' @param time integer vector of time points in \code{1:L} at which to compute the dwell-time distribution. If \code{NULL}, the overall dwell-time distribution is computed.
#' @param state integer vector of state indices for which to compute the dwell-time distribution. If \code{NULL}, dwell-time distributions for all states are returned in a named list.
#'
#' @return either time-varying dwell-time distribution(s) if \code{time} is specified, or overall dwell-time distribution if \code{time} is \code{NULL}. 
#' If more than one \code{state} is specified, a named list over states is returned.
#'
#' @examples 
#' # setting parameters for trigonometric link
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(beta = beta, degree = 1)
#' 
#' # at specific times and for specific state
#' ddwell(1:20, Gamma, time = 1:4, state = 1)
#' # results in 4x20 matrix
#' 
#' # or overall distribution for all states
#' ddwell(1:20, Gamma)
#' # results in list of length 2, each element is a vector of length 20
#' 
#' @export
ddwell <- function(
    x, 
    Gamma,
    time = NULL,
    state = NULL
    ) {
  
  # check if x is a vector
  if (!is.vector(x) | any(x < 0)) {
    stop("'x' must be a vector of non-negative integers.")
  }
  
  # check if Gamma is a 3D array
  if (!is.array(Gamma) || length(dim(Gamma)) != 3 || dim(Gamma)[1] != dim(Gamma)[2]) {
    stop("'Gamma' must be a 3D array with dimensions 'c(N, N, L)', where 'N' is the number of states and 'L' is the cycle length.")
  }
  
  # assign dimensions
  L <- dim(Gamma)[3] # cycle length
  nStates <- dim(Gamma)[1] # number of states
  
  # check if state is NULL or a valid integer
  if (is.null(state)) {
    state <- seq_len(nStates) # compute for all states
  } else if (!is.numeric(state) || any(state < 1) || any(state > nStates)) {
    stop("'state' must be an integer vector in the range [1, N], where 'N' is the number of states.")
  }
  
  # create output list
  out <- list()
  
  # locally define time-varying dwell-time distribution
  ddwell_t <- function(x, t, state){
    ind <- (t + (1:max(x)) - 1) %% L
    ind[which(ind == 0)] <- L
    gamma_ii <- Gamma[state, state, ind]
    pmf <- c(1, cumprod(gamma_ii)[-length(ind)]) * (1 - gamma_ii)
    return(pmf[x])
  }
  
  # check if t is NULL or a valid vector
  if(!is.null(time)){
    if (!is.vector(time) || any(time < 1) || any(time > L)) {
      stop("'time' must be a vector of integers in the range [1, L].")
    }
    
    for(s in state) {
      thisout <- t(sapply(time, function(t) ddwell_t(x, t, s)))
      rownames(thisout) <- paste0("t", time)
      colnames(thisout) <- x
      
      if(nrow(thisout) == 1){
        thisout <- as.numeric(thisout) # if only one time point, return vector
        names(thisout) <- x
      }
      
      out[[paste("state", s)]] <- thisout
    }
  } else { # compute overall dwell-time distribution
    # compute all periodically stationary distributions
    Delta <- stationary_p(Gamma)
    
    for(s in state) {
      weights <- numeric(L) # calculate weights
      weights[1] <- sum(Delta[L, -s] * Gamma[-s, s, L])
      for (k in 2:L) { 
        weights[k] <- sum(Delta[k-1, -s] * Gamma[-s, s, k-1])
      }
      weights <- weights / sum(weights)
      
      pmfs_weighted <- matrix(NA, L, length(x))
      for(k in 1:L) { 
        pmfs_weighted[k,] <- weights[k] * ddwell_t(x, k, s) 
      }
      
      pmf <- as.numeric(colSums(pmfs_weighted))
      names(pmf) <- x
      
      out[[paste("state", s)]] <- pmf
    }
  }
  
  # check if only one state -> if so: drop top level
  if(length(out) == 1) {
    out <- out[[1]]
  } else {
    names(out) <- paste("state", state)
  }
  
  return(out)
}
