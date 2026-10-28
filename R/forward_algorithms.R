# HMMs, ct-HMMs and MMPPs -------------------------------------------------

#' Forward algorithm to calculate the HMM log-likelihood
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model
#' using the \strong{forward algorithm} (Zucchini, MacDonald & Langrock, 2016).
#'
#' @references 
#' Zucchini, W., MacDonald, I.L., & Langrock, R. (2016). 
#' \emph{Hidden Markov Models for Time Series: An Introduction Using R} 
#' (2nd ed.). Chapman & Hall/CRC.
#'
#' Mews, Koslik, & Langrock (2025). 
#' How to build your latent Markov model: The role of time and space. 
#' Statistical Modelling, 25(6), 481--507.
#'
#' @family forward algorithms
#'
#' @param delta 
#' initial distribution; either
#'   \itemize{
#'     \item a vector of length \code{nStates}, or
#'     \item a matrix of dimension \code{c(nTracks, nStates)}, if \code{trackID} is provided, or
#'     \item \code{"stationary"}, for the (periodically) stationary distribution of \code{Gamma} at the start of each track (see \sQuote{Details}).
#'   }
#'   
#' @param Gamma 
#' transition probability matrix; either
#'   \itemize{
#'     \item a matrix of dimension \code{c(nStates, nStates)},
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} if \code{trackID} is provided,
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, or
#'     \item an array of dimension \code{c(nStates, nStates, L)} for periodic transition probabilities, if \code{tod} is provided.
#'   }
#'   
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#' 
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks (see \sQuote{Details}).
#' The observations of each track need to be consecutive.
#' 
#' @param logspace 
#' logical; if \code{TRUE}, \code{allprobs} is assumed to be on the log-scale, improving numerical stability for small probabilities.
#' 
#' @param bw 
#' optional positive integer specifying the bandwidth for a banded approximation of the forward algorithm, inducing a banded Hessian w.r.t. the observations. 
#' Defaults to \code{NULL} (exact algorithm). Approximation error decays geometrically in \code{bw}.
#' Only used with automatic differentiation; otherwise the exact likelihood is computed.
#' 
#' @param report 
#' logical; if \code{TRUE} (default), \code{delta}, \code{Gamma}, \code{allprobs}, and \code{trackID} are reported from the fitted model. Requires \code{ad = TRUE}.
#' 
#' @param ad 
#' logical; whether to use automatic differentiation. Determined automatically — for debugging only.
#' 
#' @param tod 
#' optional integer vector of length \code{nObs} with values in 1, ..., L, indexing the slices of \code{Gamma} for periodic models
#' (e.g. time of day with L = 24 for hourly data).
#' \code{Gamma[,,tod[t]]} is the transition into time t, hence the periodically stationary initial distribution is \code{stationary_p(Gamma, tod[1] \%\% L + 1)}.
#'
#' @details
#' If \code{trackID} is provided, the total log-likelihood is the sum of each track's likelihood
#' contribution. In this case, \code{Gamma} can be
#'   \itemize{
#'     \item a matrix (same transition probabilities for each track),
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} (track-specific transition probability matrix), or
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, or
#'     \item an array of dimension \code{c(nStates, nStates, L)} for periodic transition probabilities, if \code{tod} is provided.
#'   }
#' Additionally, \code{delta} can be a vector (same initial distribution for each track) or a
#' matrix of dimension \code{c(nTracks, nStates)} (track-specific initial distributions).
#' 
#' With \code{delta = "stationary"}, the initial distribution of each track is computed from \code{Gamma}:
#'   \itemize{
#'     \item the stationary distribution of \code{Gamma} (or of its slice for each track, if track-specific),
#'     \item the periodically stationary distribution at the first time point of each track, if \code{tod} is provided (see \code{\link{stationary_p}}),
#'     \item the stationary distribution of the first slice of each track, for time-varying transition probabilities, i.e. at the covariate values of the first observation.
#'   }
#' This requires \code{Gamma} to contain transition probability matrices (not, e.g., sub-stochastic matrices), which is checked without automatic differentiation.
#' The computed \code{delta} is reported as a matrix of dimension \code{c(nTracks, nStates)}.
#' 
#' \strong{Note:} \code{forward} should only be called \strong{once} in the likelihood function, with a \code{trackID} argument for multiple tracks.
#' Downstream functions like \code{\link{lama}}, \code{\link{viterbi}}, \code{\link{stateprobs}} or \code{\link{pseudo_res}} rely on the quantities it reports,
#' and only those of the last call are kept.
#'
#' @return HMM log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' ## negative log likelihood function
#' nll = function(par, step) {
#'  # parameter transformations for unconstrained optimisation
#'  Gamma = tpm(par[1:2]) # multinomial logit link
#'  delta = stationary(Gamma) # stationary HMM
#'  mu = exp(par[3:4])
#'  sigma = exp(par[5:6])
#'  # calculate all state-dependent probabilities
#'  allprobs = matrix(1, length(step), 2)
#'  ind = which(!is.na(step))
#'  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'  # simple forward algorithm to calculate log-likelihood
#'  -forward(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the trex data
#' par = c(-2,-2,            # initial tpm params (logit-scale)
#'         log(c(0.3, 2.5)), # initial means for step length (log-transformed)
#'         log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
#' mod = nlm(nll, par, step = trex$step[1:1000])
forward <- function(delta,
                    Gamma,
                    allprobs,
                    trackID = NULL,
                    logspace = FALSE,
                    bw = NULL,
                    report = TRUE,
                    ad = NULL,
                    tod = NULL
                    ){

  # streams from statedens(), cleared for the next evaluation of the likelihood
  force(allprobs) # statedens() may be called in the argument, so evaluate it before reading the streams
  streams <- .lama_env$streams
  .lama_env$streams <- NULL

  # obj$simulate() with data from OBS(): simulate the states and the data instead of evaluating the likelihood
  if(any(vapply(streams, function(s) inherits(s$x, "simref"), logical(1)))) {
    return(sim_forward(streams, delta, Gamma, allprobs, trackID, logspace, tod))
  }

  # record mode of pseudo_res(): store allprobs instead of computing the likelihood
  if(!is.null(.lama_env$record)) {
    .lama_env$record[[length(.lama_env$record) + 1]] <- if(logspace) exp(allprobs) else allprobs
    return(0)
  }

  # check inputs and bring them into the form used by the algorithm
  inp <- prepare_forward(delta, Gamma, allprobs, trackID, logspace, bw, gidx = tod)

  # report quantities for easy use later
  if(report) {
    if(logspace) allprobs <- exp(allprobs)
    if(is.character(delta)) delta <- inp$delta # computed stationary distribution(s)
    if(!is.null(tod)) REPORT(tod)
    if(!is.null(streams)) REPORT(streams)
    report_forward(delta, Gamma, allprobs, inp$trackID, inp$type)
  }

  run_forward(inp, ad)
}

# Checks inputs of the forward algorithm and brings them into a common form:
# - delta: matrix of dimension c(nTracks, nStates)
# - Gamma: array of dimension c(nStates, nStates, nSlices)
# - gidx: slice of Gamma for the transition into each time point
# - start, end: first and last index of each track
# - allprobs: if logspace, divided by its row maxima (shift), added back on the log scale
prepare_forward <- function(delta, Gamma, allprobs, trackID = NULL,
                            logspace = FALSE, bw = NULL, gidx = NULL){

  # Check allprobs
  if(!is.matrix(allprobs)) {
    stop("allprobs needs to be a matrix of dimension c(n, N).")
  }
  nObs <- nrow(allprobs)
  nStates <- ncol(allprobs)

  # check bw
  if (!is.null(bw)) {
    if (length(bw) != 1 || is.na(bw) || bw != as.integer(bw) || bw < 2) {
      stop("bw must be a single integer >= 2 (bandwidth of 1 is degenerate).")
    }
    bw <- as.integer(bw)
  }

  # check trackID
  if (!is.null(trackID)) {
    if (length(trackID) != nObs) {
      stop("Length of trackID must equal number of rows in allprobs.")
    }
    # integer track indices 1, 2, ..., nTracks in order of appearance
    trackID <- match(trackID, unique(trackID))
    if (is.unsorted(trackID)) {
      stop("Tracks must be contiguous, i.e. the observations of each track must be consecutive. Order your data by trackID.")
    }
    nTracks <- max(trackID)
    start <- which(!duplicated(trackID))
  } else {
    nTracks <- 1L
    start <- 1L
  }
  end <- c(start[-1] - 1L, nObs)

  ## handle Gamma
  Gamma <- AD(as.array(Gamma))
  dG <- dim(Gamma)

  # allow matrix input -> convert to 3D array
  if (length(dG) == 2) {
    Gamma <- AD(array(Gamma, dim = c(dG, 1)))
    dG <- dim(Gamma)
  }

  # basic square check
  if (!all(dG[1:2] == nStates)) {
    stop("Gamma must be square with dimensions (N, N).")
  }

  # gidx[t]: slice of Gamma for the transition into t (unused at the start of a track)
  type <- "homogeneous"
  if (!is.null(gidx)) { # tod, for periodic models
    if (length(gidx) != nObs) {
      stop("tod must have the same length as the number of rows in allprobs.")
    }
    if (anyNA(gidx) || any(gidx < 1 | gidx > dG[3] | gidx != round(gidx))) {
      stop("tod must contain integers in 1, ..., L, indexing the slices of Gamma (array of dimension c(N,N,L)).")
    }
    type <- "periodic"
  } else if (dG[3] == 1) {
    gidx <- rep(1L, nObs)
  } else if (!is.null(trackID) && dG[3] == nTracks) { # track-specific
    gidx <- trackID
    type <- "track-specific"
  } else if (dG[3] == nObs) { # time-varying
    gidx <- seq_len(nObs)
    type <- "inhomogeneous"
  } else if (is.null(trackID) && dG[3] == nObs - 1) { # time-varying, without first slice
    gidx <- c(1L, seq_len(nObs - 1))
    type <- "inhomogeneous"
  } else {
    stop("Gamma must have third dimension 1, nTracks (track-specific), nObs (time-varying),",
         " or nObs - 1 (time-varying, only without trackID).")
  }

  ## handle delta
  if (is.character(delta)) {
    if (!identical(delta, "stationary")) stop("delta must be numeric or \"stationary\".")
    delta <- stationary_delta(Gamma, gidx[start], type == "periodic")
  }
  delta <- AD(as.matrix(delta))
  d <- dim(delta)

  # allow N x 1 and transpose
  if (d[1] == nStates && d[2] == 1) {
    delta <- t(delta)
    d <- dim(delta)
  }

  # check dimensions depending on trackID
  if (is.null(trackID)) {
    # single track: must be 1 x N
    if (!(d[1] == 1 && d[2] == nStates)) {
      stop("delta must be a vector of length N or a matrix of size (1, N) or (N, 1) when no trackID is provided.")
    }
  } else {
    # multiple tracks: allow vector, 1xN, or kxN
    if (d[1] == 1 && d[2] == nStates) {
      delta <- delta[rep(1, nTracks), , drop = FALSE]
    } else if (!(d[1] == nTracks && d[2] == nStates)) {
      stop("delta must have dimensions (k, N), (1, N), or be a vector of length N when trackID is provided.")
    }
  }

  # if allprobs on log-scale, divide each row by its maximum:
  # exact, as each row scales the likelihood by a constant, and numerically stable
  shift <- 0
  if (logspace) {
    shift <- rowmax(allprobs)
    allprobs <- exp(allprobs - shift)
  }

  list(delta = delta, Gamma = Gamma, allprobs = allprobs, shift = shift,
       gidx = as.integer(gidx), start = start, end = end,
       trackID = trackID, type = type, bw = bw)
}

# delta = "stationary": initial distribution of each track (matrix c(nTracks, nStates)),
# given the slices g of Gamma at the track starts
# - periodic: Gamma[,,s] leads into s, so the distribution at s is row s + 1 of stationary_p()
#   (where Gamma[,,s] leads out of s)
# - else: stationary distribution of slice g (for time-varying Gamma, at the first covariate values)
stationary_delta <- function(Gamma, g, periodic = FALSE){
  # check only possible without AD (e.g. when the likelihood is run once before MakeADFun())
  if (!ad_context() && any(abs(colSums(aperm(Gamma, c(2, 1, 3))) - 1) > 1e-8)) {
    stop("delta = \"stationary\" needs transition probability matrices in Gamma (rows summing to 1), ",
         "e.g. not sub-stochastic matrices of MMPPs or state-space models.")
  }
  L <- dim(Gamma)[3]
  if (periodic && L > 1) {
    return(stationary_p(Gamma)[g %% L + 1, , drop = FALSE])
  }
  u <- unique(g) # each slice only once
  stationary(Gamma[, , u, drop = FALSE])[match(g, u), , drop = FALSE]
}

# REPORT quantities of the forward algorithm (with these names)
report_forward <- function(delta, Gamma, allprobs, trackID, type){
  REPORT(delta)
  REPORT(Gamma)
  REPORT(allprobs)
  if(!is.null(trackID)) REPORT(trackID)
  REPORT(type)
}

# Weights of the log-likelihood contributions log p(x_t | x_1, ..., x_t-1), only used by scores():
# - counting mode (.lama_env$count not NULL): counts the observations of all forward() calls
# - weights (.lama_env$weights): consecutive forward() calls take consecutive blocks
get_weights <- function(nObs){
  if(!is.null(.lama_env$count)) .lama_env$count <- .lama_env$count + nObs
  w <- .lama_env$weights
  if(is.null(w)) return(NULL)
  pos <- .lama_env$pos
  if(pos + nObs > length(w)) stop("More observations in forward() than weights (internal error, please report).")
  .lama_env$pos <- pos + nObs
  w[pos + seq_len(nObs)]
}

# Runs the forward algorithm on the output of prepare_forward()
# - if ad == FALSE, use C++ version for speed
# - if ad == TRUE, use R version which is converted into diff'able C++ code via RTMB
run_forward <- function(inp, ad = NULL){

  # if ad is not explicitly provided, check via RTMB:::ad_context()
  if(is.null(ad)){
    ad <- ad_context()
  }

  nStates <- ncol(inp$allprobs)

  # weights of the log-likelihood contributions, only set by scores() while taping (AD only)
  w <- get_weights(nrow(inp$allprobs))
  if(!is.null(w) && !ad) stop("Weights are only used with AD (internal error, please report).")
  wsum <- function(x) if(is.null(w)) sum(x) else sum(w * x)

  # value check, only possible without AD
  # (no row sum checks: MMPPs and state-space models use sub-stochastic Gamma)
  if(!ad && anyNA(inp$allprobs)) {
    stop("allprobs contains NA values. For missing observations, keep the rows of allprobs at 1, ",
         "e.g. by only filling allprobs[ind, ] with ind = which(!is.na(x)).")
  }
  if(!ad && any(inp$allprobs < 0)) {
    stop("allprobs contains negative values. For log-densities, use logspace = TRUE.")
  }

  # only one state: no hidden process
  if(nStates == 1) {
    return(as.numeric(wsum(log(inp$allprobs[, 1])) + wsum(inp$shift)))
  }

  if(!ad) { # non-AD version in C++, always exact (bw is ignored)

    l <- forward_cpp(inp$allprobs, inp$delta, inp$Gamma, inp$gidx - 1L, inp$start - 1L)

  } else { # AD version -> R code

    # printing suggestion to change TapeConfig (once per session)
    cfg <- RTMB::TapeConfig()
    if(!isTRUE(.lama_env$tape_msg_shown)) {
      if(cfg$matmul != "plain" & nStates <= 5) {
        message("Performance tip: Consider running TapeConfig(matmul = 'plain') before MakeADFun() to speed up the forward algorithm.\n")
        .lama_env$tape_msg_shown <- TRUE
      }
      if(cfg$matmul == "plain" & nStates > 50) {
        msg <- paste0("Your model has a lot of states (", nStates, "). ",
                      "Run TapeConfig(matmul = 'atomic') or TapeConfig(matmul = 'compact') before MakeADFun() to speed up the forward algorithm.\n")
        message(msg)
        .lama_env$tape_msg_shown <- TRUE
      }
    }
    if (!is.null(inp$bw) && inp$bw <= 5) {
      warning("bw = ", inp$bw, " is very small; the banded approximation may be ",
              "poor. Consider a larger bandwidth.", call. = FALSE)
    }

    l <- forward_ad(inp, w)
  }

  as.numeric(l + wsum(inp$shift))
}

# forward algorithm in R, taped by RTMB for automatic differentiation
# w: optional weights of the log-likelihood contributions (see run_forward())
forward_ad <- function(inp, w = NULL){

  Gamma <- inp$Gamma
  bw <- inp$bw
  nStates <- ncol(inp$allprobs)

  # log-likelihood contribution of observation t, weighted if w is given
  lc <- function(sumfoo, t) if(is.null(w)) log(sumfoo) else w[t] * log(sumfoo)

  # Initialising log-likelihood at 0
  l <- 0

  for(i in seq_along(inp$start)) {
    ind <- inp$start[i]:inp$end[i] # indices of track i

    allprobs_i <- inp$allprobs[ind, , drop = FALSE] # matrix
    delta_i <- inp$delta[i, , drop = FALSE] # matrix
    gidx_i <- inp$gidx[ind] # Gamma slices of track i

    # Initialize forward algorithm
    foo <- delta_i * allprobs_i[1, , drop = FALSE]
    sumfoo <- sum(foo)
    phi <- foo / sumfoo
    l <- l + lc(sumfoo, ind[1])

    if(is.null(bw)) {
      k <- length(ind)
    } else {
      k <- bw
    }

    # Regular forward algorithm
    # If banded, only runs first block
    for(t in seq_len(min(k, length(ind)))[-1]) {
      foo <- (phi %*% Gamma[, , gidx_i[t]]) * allprobs_i[t, , drop = FALSE]
      sumfoo <- sum(foo)
      phi <- foo / sumfoo
      l <- l + lc(sumfoo, ind[t])
    }

    # If banded, the approximation is as follows:
    # - separate the log likelihood into blocks of size bw
    # - for each block, starting at time t
    #   - initialise with uniform distribution at t-bw
    #   - use observations from t-bw to t-1 to update (scaled) forward variable to get close approx to true one
    #   - compute likelihood of block from t to t + bw using this approximation
    #   - this way, blocks are independent after a time lag of 2 * bw (at most)

    # If banded, run this banded version below for remaining blocks
    if(!is.null(bw) & (k < length(ind))) {
      stateDist <- rep(1, nStates) / nStates

      startInd <- k + 1
      endInd <- 2 * k
      nBlocks <- ceiling(length(ind) / k) - 1 # number of remaining blocks

      for(b in seq_len(nBlocks)) {
        # updating forward variable to get good approximation
        foo <- stateDist * allprobs_i[startInd-k, , drop = FALSE]
        foo <- foo / sum(foo)
        for(t in (startInd-k+1):(startInd - 1)) {
          foo <- (foo %*% Gamma[, , gidx_i[t]]) * allprobs_i[t, , drop = FALSE]
          foo <- foo / sum(foo)
        }
        # computing likelihood of this block based on approx forward variable
        end <- min(endInd, length(ind))
        for(t in startInd:end) {
          foo <- (foo %*% Gamma[, , gidx_i[t]]) * allprobs_i[t, , drop = FALSE]
          sumfoo <- sum(foo)
          foo <- foo / sumfoo
          l <- l + lc(sumfoo, ind[t])
        }
        startInd <- startInd + k
        endInd <- endInd + k
      }
    }
  }
  l
}


#' Forward algorithm with time-varying transition probability matrix
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model with time-varying transition probabilities using the \strong{forward algorithm} (Zucchini, MacDonald & Langrock, 2016).
#' 
#' \strong{Deprecated:} please use \code{\link{forward}} for all model types instead.
#'
#' @references 
#' Zucchini, W., MacDonald, I.L., & Langrock, R. (2016). 
#' \emph{Hidden Markov Models for Time Series: An Introduction Using R} 
#' (2nd ed.). Chapman & Hall/CRC.
#'
#' Mews, Koslik, & Langrock (2025). 
#' How to build your latent Markov model: The role of time and space. 
#' Statistical Modelling, 25(6), 481--507.
#'
#' @family forward algorithms
#' 
#' @param delta 
#' initial distribution; either
#'   \itemize{
#'     \item a vector of length \code{nStates}, or
#'     \item a matrix of dimension \code{c(nTracks, nStates)}, if \code{trackID} is provided. 
#'   }
#'   
#' @param Gamma
#' array of transition probability matrices of dimension \code{c(nStates, nStates, nObs)}, where
#' the first slice of each track is ignored as there is no transition into the start of a track.
#' 
#' For a single track, an array of dimension \code{c(nStates, nStates, nObs-1)} is also accepted.
#'
#' This function also supports continuous-time HMMs, where each slice is a Markov semigroup
#' \eqn{\Gamma(\Delta t) = \exp(Q \Delta t)} for generator \eqn{Q}.
#' 
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#' 
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks (see \sQuote{Details}).
#' The observations of each track need to be consecutive.
#' 
#' @param logspace 
#' logical; if \code{TRUE}, \code{allprobs} is assumed to be on the log-scale, improving numerical stability for small probabilities.
#' 
#' @param bw 
#' optional positive integer specifying the bandwidth for a banded approximation of the forward algorithm, inducing a banded Hessian w.r.t. the observations. 
#' Defaults to \code{NULL} (exact algorithm). Approximation error decays geometrically in \code{bw}.
#' Only used with automatic differentiation; otherwise the exact likelihood is computed.
#' 
#' @param report 
#' logical; if \code{TRUE} (default), \code{delta}, \code{Gamma}, \code{allprobs}, and \code{trackID} are reported from the fitted model. Requires \code{ad = TRUE}.
#' 
#' @param ad 
#' logical; whether to use automatic differentiation. Determined automatically — for debugging only.
#' 
#' @details
#' If \code{trackID} is provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{Gamma} must be an array of dimension \code{c(nStates, nStates, nObs)}, matching the number of rows of allprobs. For each track, the transition matrix at the beginning of the track will be ignored (as there is no transition between tracks).
#' Additionally, \code{delta} can be a vector (same initial distribution for each track) or a matrix of dimension \code{c(nTracks, nStates)} (different initial distribution for each track).
#' 
#' \strong{Note:} When there are multiple tracks, for compatibility with downstream functions like \code{\link{viterbi_g}}, \code{\link{stateprobs_g}} or \code{\link{pseudo_res}}, 
#' \code{forward_g} should only be called \strong{once} with a \code{trackID} argument.
#'
#' @return HMM log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' ## Simple usage
#' Gamma = array(c(0.9, 0.2, 0.1, 0.8), dim = c(2,2,10))
#' delta = c(0.5, 0.5)
#' allprobs = matrix(0.5, 10, 2)
#' forward_g(delta, Gamma, allprobs)
#' \donttest{
#' ## Full model fitting example
#' ## negative log likelihood function
#' nll = function(par, step, Z) {
#'  # parameter transformations for unconstrained optimisation
#'  beta = matrix(par[1:6], nrow = 2)
#'  Gamma = tpm_g(Z, beta) # multinomial logit link for each time point
#'  delta = stationary(Gamma[,,1]) # stationary HMM
#'  mu = exp(par[7:8])
#'  sigma = exp(par[9:10])
#'  # calculate all state-dependent probabilities
#'  allprobs = matrix(1, length(step), 2)
#'  ind = which(!is.na(step))
#'  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'  # simple forward algorithm to calculate log-likelihood
#'  -forward_g(delta, Gamma, allprobs)
#' }
#' 
#' ## fitting an HMM to the trex data
#' par = c(-1.5,-1.5,        # initial tpm intercepts (logit-scale)
#'         rep(0, 4),        # initial tpm slopes
#'         log(c(0.3, 2.5)), # initial means for step length (log-transformed)
#'         log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
#' mod = nlm(nll, par, step = trex$step[1:500], Z = cosinor(trex$tod[1:500]))
#' }
forward_g <- function(delta, 
                      Gamma, 
                      allprobs, 
                      trackID = NULL, 
                      logspace = FALSE,
                      bw = NULL,
                      report = TRUE,
                      ad = NULL
                      ){
  
  deprecated(
    "forward_g() is deprecated, please use forward() for all model types instead."
  )
  
  # forward() handles all types of Gamma
  forward(delta, Gamma, allprobs, trackID = trackID, logspace = logspace, 
          bw = bw, report = report, ad = ad)
}

#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with for periodically varying transition probability matrices
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model with periodically varying transition probabilities using the \strong{forward algorithm}.
#' 
#' \strong{Deprecated:} please use \code{\link{forward}} with argument \code{tod} instead.
#' 
#' @family forward algorithms
#' 
#' @details
#' When the transition probability matrix only varies periodically (e.g. as a function of time of day), there are only \eqn{L} unique matrices if \eqn{L} is the period length (e.g. \eqn{L=24} for hourly data and time-of-day variation).
#' Thus, it is much more efficient to only calculate these \eqn{L} matrices and index them by a time variable (e.g. time of day or day of year) instead of calculating such a matrix for each index in the data set (which would be redundant).
#' This function allows for that by only expecting a transition probability matrix for each time point in a period and an integer valued (\eqn{1, \dots, L}) time variable that maps the data index to the according time.
#'
#' @param delta initial or stationary distribution of length \code{nStates}, or matrix of dimension \code{c(nTracks, nStates)} for \code{nTracks} independent tracks, if \code{trackID} is provided
#' @param Gamma array of transition probability matrices of dimension \code{c(nStates, nStates, L)}.
#' 
#' Here we use the definition \eqn{\Pr(S_t=j \mid S_{t-1}=i) = \gamma_{ij}^{(t)}} such that the transition probabilities between time point \eqn{t-1} and \eqn{t} are an element of \eqn{\Gamma^{(t)}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)}
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length \code{nObs})
#' 
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions of dimension \code{c(nTracks, nObs)}, can be provided, such that each track starts with it's own initial distribution.
#' @param ad optional logical, indicating whether automatic differentiation with \code{RTMB} should be used. By default, the function determines this itself.
#' @param report logical, indicating whether \code{delta}, \code{Gamma}, \code{allprobs}, and potentially \code{trackID} should be reported from the fitted model. 
#' Defaults to \code{TRUE}, but only works if \code{ad = TRUE}, as it uses the \code{RTMB} package. 
#' 
#' \strong{Caution:} When there are multiple tracks, for compatibility with downstream functions like \code{\link{viterbi_p}}, \code{\link{stateprobs_p}} or \code{\link{pseudo_res}}, 
#' \code{forward_p} should only be called \strong{once} with a \code{trackID} argument.
#' @param logspace logical, indicating whether the probabilities/ densities in the \code{allprobs} matrix are on log-scale. If so, internal computations are also done on log-scale which is numerically more robust when the entries are very small.
#'
#' @return HMM log-likelihood for given data and parameters
#' @export
#' @importFrom RTMB REPORT
#'
#' @examples
#' ## negative log likelihood function
#' nll = function(par, step, tod) {
#'  # parameter transformations for unconstrained optimisation
#'  beta = matrix(par[1:6], nrow = 2)
#'  Gamma = tpm_p(1:24, beta = beta) # multinomial logit link for each time point
#'  delta = stationary_p(Gamma, tod[1] %% 24 + 1) # stationary HMM
#'  mu = exp(par[7:8])
#'  sigma = exp(par[9:10])
#'  # calculate all state-dependent probabilities
#'  allprobs = matrix(1, length(step), 2)
#'  ind = which(!is.na(step))
#'  for(j in 1:2) allprobs[ind,j] = dgamma2(step[ind], mu[j], sigma[j])
#'  # simple forward algorithm to calculate log-likelihood
#'  -forward_p(delta, Gamma, allprobs, tod)
#' }
#' 
#' ## fitting an HMM to the nessi data
#' par = c(-2,-2,            # initial tpm intercepts (logit-scale)
#'         rep(0, 4),        # initial tpm slopes
#'         log(c(0.3, 2.5)), # initial means for step length (log-transformed)
#'         log(c(0.2, 1.5))) # initial sds for step length (log-transformed)
#' mod = nlm(nll, par, step = trex$step[1:500], tod = trex$tod[1:500])
forward_p <- function(delta, 
                      Gamma, 
                      allprobs, 
                      tod, 
                      trackID = NULL, 
                      ad = NULL, 
                      report = TRUE,
                      logspace = FALSE){
  
  deprecated(
    "forward_p() is deprecated, please use forward(..., tod = tod) for periodic models instead."
  )
  
  # forward() handles periodic Gamma, indexed by tod
  forward(delta, Gamma, allprobs, trackID = trackID, logspace = logspace, 
          report = report, ad = ad, tod = tod)
}
# HSMMs -------------------------------------------------------------------


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for homogeneous hidden semi-Markov models
#' 
#' Calculates the (approximate) log-likelihood of a sequence of observations under a homogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#'
#' @details
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @references 
#' Langrock, R., & Zucchini, W. (2011). Hidden Markov models with arbitrary state dwell-time distributions. Computational Statistics & Data Analysis, 55(1), 715-724.
#'
#' Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
#'
#' @param dm list of length \code{nStates} containing vectors of dwell-time probability mass functions (PMFs) for each state. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' @param omega matrix of dimension \code{c(nStates, nStates)} of conditional transition probabilites, also called embedded transition probability matrix. 
#' 
#' Contains the transition probabilities given that the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Can be constructed using \code{\link{tpm_emb}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)} which will automatically be converted to the appropriate dimension.
#' @param trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{dm} can be a nested list, where the top layer contains k \code{dm} lists as described above. \code{omega} can then also be an array of dimension \code{c(nStates, nStates, nTracks)} with one conditional transition probability matrix for each track.
#' Furthermore, instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension \code{c(nTracks, nStates)}, can be provided, such that each track starts with it's own initial distribution.
#' @param delta optional vector of initial state probabilities of length \code{nStates}
#' 
#' By default, the stationary distribution is computed (which is typically recommended).
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return HSMM log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_hsmm <- function(dm, omega, allprobs,
                         trackID = NULL, delta = NULL, eps = 1e-10, report = TRUE){
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  agsizes = vapply(dm, length, 1L)
  
  N = ncol(allprobs) # number of HSMM states
  M = sum(agsizes) # total number of states of the approximating HMM
  
  stateInds <- rep(1:N, times = agsizes)
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    ## compute approximating tpm
    Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
    # Gamma = Matrix::Matrix(Gamma_sparse, sparse = FALSE) # dense for reporting
    
    ## if stationary, compute initial stationary distribution
    if(stationary){
      delta = stationary_sparse(Gamma_sparse)
      delta_sparse = methods::as(t(delta), "sparseMatrix")
    } else{ # if delta is provided, "stuff out" with zeros
      cols_to_fill = c(1, cumsum(agsizes[-N])+1)
      delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                          x = delta, dims = c(1, M))
    }
    
    ## report quantities for state decoding
    if(report){
      RTMB::REPORT(delta_sparse)
      RTMB::REPORT(Gamma_sparse)
      RTMB::REPORT(allprobs)
    }
    
    ## making AD work fine
    delta_sparse <- AD(delta_sparse)
    Gamma_sparse <- AD(Gamma_sparse)
    allprobs <- AD(allprobs)
    
    # forward algorithm
    # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
    foo = delta_sparse * allprobs[1, stateInds, drop = FALSE]
    # foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
    # sumfoo = sum(foo)
    sumfoo = rowSums(foo)
    phi = foo / sumfoo
    l = log(sumfoo)
    
    for(t in 2:nrow(allprobs)) {
      # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
      # foo = phi %*% Gamma_sparse %*% diag(rep(allprobs[t,], times = agsizes))
      foo = (phi %*% Gamma_sparse) * allprobs[t, stateInds, drop = FALSE]
      # sumfoo = sum(foo)
      sumfoo = rowSums(foo)
      phi = foo / sumfoo
      l = l + log(sumfoo)
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    if(length(dm) == k){ # dm for each track
      if(is.matrix(omega)){ # different dms but same omegas
        omega = array(omega, dim = c(dim(omega), k))
      } else if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != k) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number of tracks.")
      }
      
      Gamma_sparse = lapply(1:k, function(k) tpm_hsmm(omega[,,k], dm[[k]], eps = eps)) # build k Gammas based on omega and dm
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[k]]))) # build k deltas based on Gammas
        Delta_sparse = methods::as(Delta, "sparseMatrix")
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
    } else if(length(dm) == length(agsizes)){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, k)
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != k) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,k), matching the number of tracks.")
        
        Gamma_sparse = lapply(1:k, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[k]]))) # build k deltas based on Gammas
          Delta_sparse = methods::as(Delta, "sparseMatrix")
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
    }
    
    ## report quantities for state decoding
    if(report){
      RTMB::REPORT(Delta_sparse)
      RTMB::REPORT(Gamma_sparse)
      RTMB::REPORT(allprobs)
    }
    
    ## making AD work fine
    allprobs <- AD(allprobs)
    
    
    ## forward algorithm
    l = 0 # initialize log-likelihood
    for(i in 1:k) {
      ind = which(trackID == uID[i]) # indices of track i
      
      delta_i = AD(Delta_sparse[i, , drop = FALSE])
      Gamma_i = AD(Gamma_sparse[[i]])
      
      # foo = delta_i %*% Matrix::Diagonal(x = rep(allprobs[ind[1],], times = agsizes))
      foo = delta_i * allprobs[ind[1], stateInds, drop = FALSE]
      # sumfoo = sum(foo)
      sumfoo = rowSums(foo)
      phi = foo / sumfoo
      l_i = log(sumfoo)
      
      for(t in 2:length(ind)) {
        # foo = phi %*% Gamma_i %*% Matrix::Diagonal(x = rep(allprobs[ind[t],], times = agsizes))
        foo = (phi %*% Gamma_i) * allprobs[ind[t], stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l_i = l_i + log(sumfoo)
      }
      
      l = l + l_i
    }
  }
  
  as.numeric(l)
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with inhomogeneous state durations and/ or conditional transition probabilities
#'
#' Calculates the (approximate) log-likelihood of a sequence of observations under an inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @details
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @references Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
#'
#' @param dm list of length \code{nStates} containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#' 
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have n rows, where n is the number of observations. The number of columns again correponds to the size of the approximating state aggregates.
#' 
#' In the latter case, the first \code{max(sapply(dm, ncol)) - 1} observations will not be used because the first approximating transition probability matrix needs to be computed based on the first \code{max(sapply(dm, ncol))} covariate values (represented by \code{dm}).
#' @param omega matrix of dimension \code{c(nStates, nStates)} or array of dimension \code{c(nStates, nStates, nObs)} of conditional transition probabilites, also called embedded transition probability matrix.
#' 
#' It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)}
#' @param trackID trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension \code{c(nTracks, nStates)}, can be provided, such that each track starts with it's own initial distribution.
#' @param delta optional vector of initial state probabilities of length N
#' 
#' By default, instead of this, the stationary distribution is computed corresponding to the first approximating transition probability matrix of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
#' @param startInd optional integer index at which the forward algorithm starts. 
#' 
#' When approximating inhomogeneous HSMMs by inhomogeneous HMMs, the first transition probability matrix that can be constructed is at time \code{max(sapply(dm, ncol))} (as it depends on the previous covariate values).
#' Hence, when not provided, \code{startInd} is chosen to be \code{max(sapply(dm, ncol))}. Fixing \code{startInd} at a value \strong{larger} than max(aggregate sizes) is useful when models with different aggregate sizes are fitted to the same data and are supposed to be compared. In that case it is important that all models use the same number of observations.
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return HSMM log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_ihsmm <- function(dm, omega, allprobs,
                          trackID = NULL, delta = NULL, startInd = NULL,
                          eps = 1e-10, report = TRUE){
  
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  N = ncol(allprobs) # number of HSMM states
  
  # obtain aggregate sizes based on dm
  if(is.matrix(dm[[1]])){
    agsizes = sapply(dm, ncol)
  } else if(is.vector(dm[[1]])){
    agsizes = sapply(dm, length)
  }
  
  stateInds <- rep(1:N, times = agsizes)
  
  M = sum(agsizes) # total number of states of the approximating HMM
  n = nrow(allprobs) # number of observations
  maxag = max(agsizes) # maximum number of states in the approximating HMM
  # this is important, because the forward algo can only start at maxag because the first 1:(maxag-1) data points are needed for the first tpm
  
  # check if startInd is provided
  # if not, take the maximum aggregate size as startInd
  if(is.null(startInd)){
    startInd = maxag
  } else {
    if(startInd < maxag) stop("startInd must be at least the maximum aggregate size")
  }
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_ihsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_sparse(Gamma_sparse[[startInd - maxag + 1]])
        delta_sparse = methods::as(t(delta), "sparseMatrix")
      } else{ # if delta is provided, "stuff out" with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                            x = delta, dims = c(1, M))
      }
      
      ## report quantities for state decoding
      if(report) {
        RTMB::REPORT(delta_sparse)
        RTMB::REPORT(Gamma_sparse)
        RTMB::REPORT(allprobs)
      }
      
      ## making AD work fine
      delta_sparse <- AD(delta_sparse)
      # Gamma_sparse <- AD(Gamma_sparse)
      allprobs <- AD(allprobs)
      
      # forward algorithm
      # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
      # foo = delta_sparse %*% diag(rep(allprobs[startInd,], times = agsizes))
      foo = delta_sparse * allprobs[startInd, stateInds, drop = FALSE]
      # sumfoo = sum(foo)
      sumfoo = rowSums(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in (startInd + 1):nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        # foo = phi %*% Gamma_sparse[[t - maxag]] %*% diag(rep(allprobs[t,], times = agsizes))
        foo = (phi %*% Gamma_sparse[[t - maxag]]) * allprobs[t, stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(is.vector(dm[[1]])){
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE)
          delta_sparse = Matrix::Matrix(delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = delta, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, n)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
        
        Gamma_sparse = lapply(1:n, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          delta = stationary_sparse(Gamma_sparse[[1]]) # build k deltas based on Gammas
          delta_sparse = methods::as(t(delta), "sparseMatrix")
        } else{
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                              x = delta, dims = c(1, M))
        }
      }
      
      ## forward algorithm
      # foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      foo = delta_sparse * allprobs[1, stateInds, drop = FALSE]
      # sumfoo = sum(foo)
      sumfoo = rowSums(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse[[t-1]] %*% diag(rep(allprobs[t,], times = agsizes))
        foo = (phi %*% Gamma_sparse[[t-1]]) * allprobs[t, stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){ # dm varies
      if(nrow(dm[[1]]) != n) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      if(is.matrix(omega)){ # fixed omega
        omega = array(omega, dim = c(dim(omega), n))
      } else if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
      }
      
      # nested list of sparse approximating tpms
      Gamma_sparse = lapply(1:k, function(i){
        ind = which(trackID == uID[i])
        this_dm = lapply(1:N, function(j) dm[[j]][ind,])
        tpm_ihsmm(omega[,,ind], this_dm, eps = eps)
      })
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        trackInd = calc_trackInd(trackID)
        Delta = t(sapply(1:k, function(i) stationary_sparse(Gamma_sparse[[i]][[1]]))) # build k deltas based on Gammas
        Delta_sparse = methods::as(Delta, "sparseMatrix")
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
      
      ## making AD work fine
      allprobs <- AD(allprobs)
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = AD(Delta_sparse[i, , drop = FALSE])
        Gamma_i = AD(Gamma_sparse[[i]])
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        # foo = delta_i %*% diag(rep(allprobs[ind[startInd],], times = agsizes))
        foo = delta_i * allprobs[ind[startInd], stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in (startInd + 1):length(ind)) {
          foo = phi %*% Gamma_i[[t-maxag]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          foo = (phi %*% Gamma_i[[t-maxag]]) * allprobs[ind[t], stateInds, drop = FALSE]
          # sumfoo = sum(foo)
          sumfoo = rowSums(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
      
      ## second case: dm does not vary, only omega varies -> easier
    } else if(is.vector(dm[[1]])){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, n)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
        
        Gamma_sparse = lapply(1:n, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          trackInd = calc_trackInd(trackID)
          Delta = t(sapply(1:k, function(k) stationary_sparse(Gamma_sparse[[trackInd[k]]]))) # build k deltas based on Gammas
          Delta_sparse = methods::as(Delta, "sparseMatrix")
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
      
      ## making AD work fine
      allprobs <- AD(allprobs)
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = AD(Delta_sparse[i, , drop = FALSE])
        Gamma_i = Gamma_sparse[ind]
        
        # foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        foo = delta_i * allprobs[ind[1], stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          # foo = phi %*% Gamma_i[[t-1]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          foo = (phi %*% Gamma_i[[t-1]]) * allprobs[ind[t], stateInds, drop = FALSE]
          # sumfoo = sum(foo)
          sumfoo = rowSums(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
    }
  }
  
  as.numeric(l)
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with periodically inhomogeneous state durations and/ or conditional transition probabilities
#'
#' Calculates the (approximate) log-likelihood of a sequence of observations under a periodically inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelhood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary periodically.
#' In the special case of periodic variation (as compared to arbitrary covariate influence), this version is to be preferred over \code{\link{forward_ihsmm}} because it computes the \strong{correct periodically stationary distribution} and no observations are lost for the approximation.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @references Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
#'
#' @param dm list of length \code{nStates} containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#'
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probablity mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have L rows, where L is the cycle length. The number of columns again correpond to the size of the approximating state aggregates.
#' @param omega matrix of dimension \code{c(nStates, nStates)} or array of dimension \code{c(nStates, nStates, L)} of conditional transition probabilites, also called embedded transition probability matrix
#' 
#' It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)}
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length \code{nObs}).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension \code{c(nTracks, nStates)}, can be provided, such that each track starts with it's own initial distribution.
#' @param delta Optional vector of initial state probabilities of length \code{nStates}. By default, instead of this, the stationary distribution is computed corresponding to the first approximating t.p.m. of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, approximating transition probability matrix and \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#'
#' @return HSMM log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' # currently no examples
forward_phsmm <- function(dm, omega, allprobs, tod,
                          trackID = NULL, delta = NULL, eps = 1e-10, report = TRUE){
  
  # overloading assignment operators, currently necessary
  "[<-" <- ADoverload("[<-")
  "c" <- ADoverload("c")
  "diag<-" <- ADoverload("diag<-")
  ################################
  
  N = ncol(allprobs) # number of HSMM states
  
  # obtain aggregate sizes based on dm
  if(is.matrix(dm[[1]])){
    agsizes = sapply(dm, ncol)
  } else if(is.vector(dm[[1]])){
    agsizes = sapply(dm, length)
  }
  
  stateInds <- rep(1:N, times = agsizes)
  
  M = sum(agsizes) # total number of states of the approximating HMM
  n = nrow(allprobs) # number of observations
  L = length(unique(tod)) # cycle length -> number of unique tpms
  # this is important, because the forward algo can only start at maxag because the first 1:(maxag-1) data points are needed for the first tpm
  
  stationary = is.null(delta) # if delta is not provided, stationary distribution needs to be computed
  
  if(is.null(trackID)) {
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){
      if(nrow(dm[[1]]) != L) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      
      ## compute approximating tpm
      Gamma_sparse = tpm_phsmm(omega, dm, eps = eps)
      
      ## if stationary, compute initial stationary distribution
      if(stationary){
        delta = stationary_p_sparse(Gamma_sparse, t = tod[1])
        delta_sparse = methods::as(t(delta), "sparseMatrix")
      } else{ # if delta is provided, "stuff out" with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                            x = delta, dims = c(1, M))
      }
      
      ## report quantities for state decoding
      if(report) {
        RTMB::REPORT(delta_sparse)
        RTMB::REPORT(Gamma_sparse)
        RTMB::REPORT(allprobs)
      }
      
      ## making AD work fine
      delta_sparse <- AD(delta_sparse)
      # Gamma_sparse <- AD(Gamma_sparse)
      allprobs <- AD(allprobs)
      
      # forward algorithm
      # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
      # foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      foo = delta_sparse * allprobs[1, stateInds, drop = FALSE]
      # sumfoo = sum(foo)
      sumfoo = rowSums(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse %*% Matrix::Diagonal(x = rep(allprobs[t,], times = agsizes))
        # foo = phi %*% Gamma_sparse[[tod[t-1]]] %*% diag(rep(allprobs[t,], times = agsizes))
        foo = (phi %*% AD(Gamma_sparse[[tod[t-1]]])) * allprobs[t, stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
      
    } else if(is.vector(dm[[1]])){
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          delta = matrix(delta, nrow = 1, ncol = length(delta), byrow = TRUE)
          delta_sparse = Matrix::Matrix(delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = delta, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, L)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
        
        Gamma_sparse = lapply(1:L, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          delta = stationary_p_sparse(Gamma_sparse, t = tod[1])
          delta_sparse = methods::as(t(delta), "sparseMatrix")
        } else{
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          delta_sparse = Matrix::sparseMatrix(i = rep(1, N), j = cols_to_fill, 
                                              x = delta, dims = c(1, M))
        }
      }
      
      ## making AD work fine
      delta_sparse <- AD(delta_sparse)
      # Gamma_sparse <- AD(Gamma_sparse)
      allprobs <- AD(allprobs)
      
      ## forward algorithm
      # foo = delta_sparse %*% diag(rep(allprobs[1,], times = agsizes))
      foo = delta_sparse * allprobs[1, stateInds, drop = FALSE]
      # sumfoo = sum(foo)
      sumfoo = rowSums(foo)
      phi = foo / sumfoo
      l = log(sumfoo)
      
      for(t in 2:nrow(allprobs)) {
        # foo = phi %*% Gamma_sparse[[tod[t-1]]] %*% diag(rep(allprobs[t,], times = agsizes))
        foo = (phi %*% AD(Gamma_sparse[[tod[t-1]]])) * allprobs[t, stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l = l + log(sumfoo)
      }
    }
    
  } else if(!is.null(trackID)) {
    
    RTMB::REPORT(trackID) # report trackID for viterbi etc.
    
    uID = unique(trackID) # unique track IDs
    k = length(uID) # number of tracks
    
    ## first case: dm varies -> truely inhomogeneous HSMM
    if(is.matrix(dm[[1]])){ # dm varies
      if(nrow(dm[[1]]) != L) stop("dm needs to be a list of length N, either containing vectors or matrices of dwell-time PMFs.")
      if(length(dim(omega)) == 3){ # different dms and different omegas
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
      }
      
      Gamma_sparse = tpm_phsmm(omega, dm, eps = eps) # still only L unique tpms
      
      if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
        trackInd = calc_trackInd(trackID)
        Delta = stationary_p_sparse(Gamma_sparse) # calculate all periodically stationaries -> not very expensive because of recursion
        Delta = Delta[tod[trackInd], ] # corresponding initial periodically stationary for each track
        Delta_sparse = methods::as(Delta, "sparseMatrix") 
      } else{
        
        # if delta is provided
        delta = as.matrix(delta) # reshape to matrix for easier handling
        
        if(ncol(delta) != N) delta = t(delta) # transpose if necessary
        
        if(nrow(delta) == 1) {
          Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
        } else {
          Delta = delta
        }
        
        # "stuff out" delta with zeros
        cols_to_fill = c(1, cumsum(agsizes[-N])+1)
        
        row_indices <- rep(1:k, N)
        col_indices <- rep(cols_to_fill, each = k)
        
        values = as.vector(Delta)
        
        Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                             j = col_indices,   # Column indices
                                             x = values,        # Non-zero values (from dense matrix)
                                             dims = c(k, M))    # Final sparse matrix dimensions
      }
      
      ## making AD work fine
      allprobs <- AD(allprobs)
      Delta_sparse <- AD(Delta_sparse)
      # Gamma_sparse <- AD(Gamma_sparse)
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = AD(Delta_sparse[i, , drop = FALSE])
        thistod = tod[ind]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        # foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        foo = delta_i * allprobs[ind[1], stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          # foo = phi %*% Gamma_sparse[[thistod[t-1]]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          foo = (phi %*% AD(Gamma_sparse[[thistod[t-1]]])) * allprobs[ind[t], stateInds, drop = FALSE]
          # sumfoo = sum(foo)
          sumfoo = rowSums(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
      
      ## second case: dm does not vary, only omega varies -> easier
    } else if(is.vector(dm[[1]])){ # only one dm
      if(is.matrix(omega)){ # if omega is matrix, just build one Gamma
        Gamma_sparse = tpm_hsmm(omega, dm, eps = eps)
        
        ## if stationary, compute initial stationary distribution
        if(stationary){
          delta = stationary_sparse(Gamma_sparse)
          Delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE) # repeat rows for each track
          Delta_sparse = Matrix::Matrix(Delta, sparse = TRUE) # make sparse
        } else{ # if delta is provided, "stuff out" with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N]) + 1)
          
          if(is.matrix(delta)){
            values = as.vector(delta)
          } else{
            values = rep(delta, each = k)
          }
          Delta_sparse = Matrix::sparseMatrix(i = rep(1:k, N), j = rep(cols_to_fill, each=k), 
                                              x = values, dims = c(k, M))
        }
        Gamma_sparse = list(Gamma_sparse)
        Gamma_sparse = rep(Gamma_sparse, L)
        
      } else if(length(dim(omega)) == 3){ # different omegas but same dm
        if(dim(omega)[3] != L) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,L), matching the cycle length.")
        
        Gamma_sparse = lapply(1:L, function(k) tpm_hsmm(omega[,,k], dm, eps = eps)) # build k Gammas with fixed dm but varying omega
        
        if(stationary == TRUE){ # delta not provided but computed: each track starts in stationary distribution
          trackInd = calc_trackInd(trackID)
          Delta = stationary_p_sparse(Gamma_sparse) # calculate all periodically stationaries -> not very expensive because of recursion
          Delta = Delta[tod[trackInd], ] # corresponding initial periodically stationary for each track
          Delta_sparse = methods::as(Delta, "sparseMatrix") 
        } else{
          # if delta is provided
          delta = as.matrix(delta) # reshape to matrix for easier handling
          
          if(ncol(delta) != N) delta = t(delta) # transpose if necessary
          
          if(nrow(delta) == 1) {
            Delta = matrix(delta, nrow = k, ncol = N, byrow = TRUE)
          } else {
            Delta = delta
          }
          
          # "stuff out" delta with zeros
          cols_to_fill = c(1, cumsum(agsizes[-N])+1)
          
          row_indices <- rep(1:k, N)
          col_indices <- rep(cols_to_fill, each = k)
          
          values = as.vector(Delta)
          
          Delta_sparse <- Matrix::sparseMatrix(i = row_indices,   # Row indices
                                               j = col_indices,   # Column indices
                                               x = values,        # Non-zero values (from dense matrix)
                                               dims = c(k, M))    # Final sparse matrix dimensions
        }
      }
      
      ## making AD work fine
      allprobs <- AD(allprobs)
      Delta_sparse <- AD(Delta_sparse)
      # Gamma_sparse <- AD(Gamma_sparse)
      
      ## forward algorithm
      l = 0 # initialize log-likelihood
      for(i in 1:k) {
        ind = which(trackID == uID[i]) # indices of track i
        
        delta_i = AD(Delta_sparse[i, , drop = FALSE])
        thistod = tod[ind]
        
        # foo = delta_sparse %*% Matrix::Diagonal(x = rep(allprobs[1,], times = agsizes))
        # foo = delta_i %*% diag(rep(allprobs[ind[1],], times = agsizes))
        foo = delta_i * allprobs[ind[1], stateInds, drop = FALSE]
        # sumfoo = sum(foo)
        sumfoo = rowSums(foo)
        phi = foo / sumfoo
        l_i = log(sumfoo)
        
        for(t in 2:length(ind)) {
          # foo = phi %*% Gamma_sparse[[thistod[t-1]]] %*% diag(rep(allprobs[ind[t],], times = agsizes))
          foo = (phi %*% AD(Gamma_sparse[[thistod[t-1]]])) * allprobs[ind[t], stateInds, drop = FALSE]
          # sumfoo = sum(foo)
          sumfoo = rowSums(foo)
          phi = foo / sumfoo
          l_i = l_i + log(sumfoo)
        }
        
        l = l + l_i
      }
    }
  }
  
  as.numeric(l)
}
