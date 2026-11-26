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
#' The likelihood of observations \eqn{x_1, \dots, x_T} is
#' \deqn{\mathcal{L}(\boldsymbol{\theta}) = \boldsymbol{\delta}^{(1)} \mathbf{P}(x_1) \boldsymbol{\Gamma}^{(2)} \mathbf{P}(x_2) \boldsymbol{\Gamma}^{(3)} \cdots \boldsymbol{\Gamma}^{(T)} \mathbf{P}(x_T) \mathbf{1}^\top,}{L(theta) = delta^(1) P(x_1) Gamma^(2) P(x_2) Gamma^(3) ... Gamma^(T) P(x_T) 1',}
#' where \eqn{\mathbf{P}(x_t)}{P(x_t)} is a diagonal matrix with the \eqn{t}-th row of \code{allprobs} on its diagonal,
#' \eqn{\boldsymbol{\Gamma}^{(t)}}{Gamma^(t)} is the transition probability matrix into time \eqn{t},
#' and \eqn{\mathbf{1}}{1} is a row vector of ones.
#' The forward algorithm evaluates this product recursively and on the log-scale.
#'
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
  tracks <- prepare_tracks(trackID, nObs)
  trackID <- tracks$trackID
  nTracks <- tracks$nTracks
  start <- tracks$start

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
  slices <- slice_index(dG[3], nObs, tracks, gidx, "slices of Gamma (array of dimension c(N,N,L))")
  if (is.null(slices)) {
    stop("Gamma must have third dimension 1, nTracks (track-specific), nObs (time-varying),",
         " or nObs - 1 (time-varying, only without trackID).")
  }
  gidx <- slices$index
  type <- slices$type

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

  scaled <- scale_allprobs(allprobs, logspace)

  list(delta = delta, Gamma = Gamma, allprobs = scaled$allprobs, shift = scaled$shift,
       gidx = gidx, start = start, end = tracks$end,
       trackID = trackID, type = type, bw = bw)
}

# Checks trackID, returns
# - trackID: integer track indices 1, 2, ..., nTracks in order of appearance (NULL if not given)
# - nTracks, start, end: number of tracks, first and last index of each track
prepare_tracks <- function(trackID, nObs){
  if (is.null(trackID)) {
    return(list(trackID = NULL, nTracks = 1L, start = 1L, end = nObs))
  }
  if (length(trackID) != nObs) {
    stop("Length of trackID must equal number of rows in allprobs.")
  }
  trackID <- match(trackID, unique(trackID))
  if (is.unsorted(trackID)) {
    stop("Tracks must be contiguous, i.e. the observations of each track must be consecutive. Order your data by trackID.")
  }
  start <- which(!duplicated(trackID))
  end <- c(start[-1] - 1L, nObs)
  list(trackID = trackID, nTracks = max(trackID), start = start, end = end)
}

# Which slice (of nSlices) belongs to each observation? Returns the index and the type of model
# - tracks: output of prepare_tracks()
# - tod: integer time of day, for periodic models
# - what: name of the slices in error messages
# NULL if the number of slices fits no model
slice_index <- function(nSlices, nObs, tracks, tod = NULL, what = "slices"){
  hasID <- !is.null(tracks$trackID)
  type <- "homogeneous"
  if (!is.null(tod)) { # periodic
    if (length(tod) != nObs) {
      stop("tod must have the same length as the number of rows in allprobs.")
    }
    if (anyNA(tod) || any(tod < 1 | tod > nSlices | tod != round(tod))) {
      stop("tod must contain integers in 1, ..., L, indexing the ", what, ".")
    }
    index <- tod
    type <- "periodic"
  } else if (nSlices == 1) {
    index <- rep(1L, nObs)
  } else if (hasID && nSlices == tracks$nTracks) { # track-specific
    index <- tracks$trackID
    type <- "track-specific"
  } else if (nSlices == nObs) { # time-varying
    index <- seq_len(nObs)
    type <- "inhomogeneous"
  } else if (!hasID && nSlices == nObs - 1) { # time-varying, without first slice
    index <- c(1L, seq_len(nObs - 1))
    type <- "inhomogeneous"
  } else {
    return(NULL)
  }
  list(index = as.integer(index), type = type)
}

# if allprobs on log-scale, divide each row by its maximum:
# exact, as each row scales the likelihood by a constant, and numerically stable
scale_allprobs <- function(allprobs, logspace = FALSE){
  if (!logspace) return(list(allprobs = allprobs, shift = 0))
  shift <- rowmax(allprobs)
  list(allprobs = exp(allprobs - shift), shift = shift)
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
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions of dimension \code{c(nTracks, nStates)}, can be provided, such that each track starts with its own initial distribution.
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


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models
#' 
#' Calculates the (approximate) log-likelihood of a sequence of observations under a hidden semi-Markov model using a modified \strong{forward algorithm}.
#' Dwell-time distributions and conditional transition probabilities can be fixed, track-specific, time-varying or periodic.
#' 
#' @family forward algorithms
#'
#' @details
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers.
#' For direct numerical maximum likelihood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' A stay that begins at time point \code{t} is governed by \code{omega[,,t]} (which state is entered) and row \code{t} of \code{dm} (how long the stay lasts).
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' @references 
#' Langrock, R., & Zucchini, W. (2011). Hidden Markov models with arbitrary state dwell-time distributions. Computational Statistics & Data Analysis, 55(1), 715-724.
#'
#' Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
#'
#' @param dm list of length \code{nStates} containing the dwell-time probability mass functions (PMFs) of the states, where each element is either
#'   \itemize{
#'     \item a vector, if the dwell-time distribution is fixed, or
#'     \item a matrix with \code{nTracks} rows (track-specific), \code{nObs} rows (time-varying), or \code{L} rows (periodic, if \code{tod} is provided), where each row is a PMF.
#'   }
#' The vector lengths or numbers of columns correspond to the approximating state aggregate sizes, hence there should be little probability mass not covered by these.
#' @param omega matrix of dimension \code{c(nStates, nStates)} of conditional transition probabilities, also called embedded transition probability matrix, 
#' or array with \code{nTracks}, \code{nObs} or \code{L} of such matrices as slices.
#' 
#' Contains the transition probabilities given that the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Can be constructed using \code{\link{tpm_emb}} and \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)} which will automatically be converted to the appropriate dimension.
#' @param trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' In this case, \code{dm} and \code{omega} can be track-specific and \code{delta} can be a matrix of initial distributions, of dimension \code{c(nTracks, nStates)}, such that each track starts with its own initial distribution.
#' @param delta optional vector of initial state probabilities of length \code{nStates}, in which case a new stay begins at the first observation used in each track.
#' 
#' By default, the stationary distribution is computed (which is typically recommended).
#' @param eps small value to avoid numerical issues in the approximating transition matrix construction. Usually, this should not be changed.
#' @param report logical, indicating whether initial distribution, \code{dm}, \code{omega} and the \code{allprobs} matrix should be reported from the fitted model. Defaults to \code{TRUE}.
#' @param tod optional integer vector of length \code{nObs} containing the time of day (or any other cycle), with values in \code{1, ..., L}, for periodic models.
#' 
#' With periodic \code{dm}, it needs to advance by one in each row within a track.
#' @param use_all_obs logical, only relevant for time-varying \code{dm}. 
#' 
#' By default, the first \code{max(sapply(dm, ncol)) - 1} observations of each track are not used, as they may belong to stays that began before the first observation, with unknown dwell-time distributions. 
#' If \code{TRUE}, all observations are used, which needs \code{delta}.
#' @param first_obs optional integer, the first observation that is used in each track. 
#' 
#' By default \code{1}, or \code{max(sapply(dm, ncol))} for time-varying \code{dm}, in which case it cannot be smaller. 
#' Useful to fit models with different aggregate sizes or types of \code{dm} to the same observations, such that their likelihoods can be compared.
#' @param logspace logical, indicating whether \code{allprobs} contains log-probabilities. Defaults to \code{FALSE}.
#'
#' @return HSMM log-likelihood for given data and parameters
#' @export
#' @import RTMB
#'
#' @examples
#' ## negative log-likelihood of a 2-state HSMM with shifted Poisson dwell times
#' nll = function(par) {
#'   getAll(par)
#'   sizes = c(20, 20) # state aggregate sizes
#'   dm = lapply(1:2, function(i) dpois(1:sizes[i] - 1, exp(loglambda[i])))
#'   omega = matrix(c(0, 1, 1, 0), 2) # 2 states: no parameters
#'   allprobs = matrix(1, length(step), 2)
#'   ind = which(!is.na(step))
#'   for(j in 1:2) allprobs[ind, j] = dgamma2(step[ind], exp(logmu[j]), exp(logsigma[j]))
#'   -forward_hsmm(dm, omega, allprobs)
#' }
#' step = trex$step[1:200]
#' par = list(loglambda = log(c(5, 5)), logmu = log(c(0.3, 2.5)), logsigma = log(c(0.2, 1.5)))
#' obj = MakeADFun(nll, par, silent = TRUE)
#' opt = nlminb(obj$par, obj$fn, obj$gr)
forward_hsmm <- function(dm, omega, allprobs,
                         trackID = NULL, delta = NULL, eps = 1e-10, report = TRUE,
                         tod = NULL, use_all_obs = FALSE, first_obs = NULL, logspace = FALSE){
  
  inp <- prepare_hsmm(dm, omega, allprobs, trackID, delta, tod, use_all_obs, first_obs,
                      logspace = logspace, eps = eps)
  if(report) {
    report_hsmm(inp, dm, omega, if(logspace) exp(allprobs) else allprobs, trackID, tod)
    if(inp$type %in% c("homogeneous", "track-specific")) report_sparse(inp, sparse_Gammas(inp))
  }
  
  run_forward_hsmm(inp)
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with inhomogeneous state durations and/ or conditional transition probabilities
#'
#' Calculates the (approximate) log-likelihood of a sequence of observations under an inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @details
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelihood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' \strong{Deprecated:} please use \code{\link{forward_hsmm}} for all types of HSMMs instead. There, a stay that begins at \code{t} is governed by \code{omega[,,t]}, here by \code{omega[,,t-1]} (\code{omega[,,t]} is for the step from \code{t} into \code{t + 1}).
#'
#' @references Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
#'
#' @param dm list of length \code{nStates} containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#' 
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probability mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have \code{nObs} rows. The number of columns again corresponds to the size of the approximating state aggregates.
#' 
#' In the latter case, the first \code{max(sapply(dm, ncol)) - 1} observations will not be used because the first approximating transition probability matrix needs to be computed based on the first \code{max(sapply(dm, ncol))} covariate values (represented by \code{dm}).
#' @param omega matrix of dimension \code{c(nStates, nStates)} or array of dimension \code{c(nStates, nStates, nObs)} of conditional transition probabilities, also called embedded transition probability matrix.
#' 
#' It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)}
#' @param trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension \code{c(nTracks, nStates)}, can be provided, such that each track starts with its own initial distribution.
#' @param delta optional vector of initial state probabilities of length \code{nStates}
#' 
#' By default, instead of this, the stationary distribution corresponding to the first approximating transition probability matrix of each track is computed. Contrary to the homogeneous case, this is not theoretically motivated but just for convenience.
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
  
  deprecated(
    "forward_ihsmm() is deprecated, please use forward_hsmm() instead, where a stay that begins at t is governed by omega[,,t] (here: omega[,,t-1])."
  )
  
  # here, omega[,,t] is for the step from t into t + 1, in forward_hsmm() for the step into t
  omega_into <- omega
  if(length(dim(omega)) == 3) {
    n <- nrow(allprobs)
    if(dim(omega)[3] != n) stop("omega needs to be either a matrix of dimension c(N,N) or an array of dimension c(N,N,n), matching the number of observations")
    omega_into <- omega[, , c(1, seq_len(n - 1)), drop = FALSE]
  }
  
  # startInd is only used for dm varying over time
  varying <- is.matrix(dm[[1]])
  Ni <- sapply(dm, function(d) if(is.matrix(d)) ncol(d) else length(d))
  if(!is.null(startInd) && startInd < max(Ni)) stop("startInd must be at least the maximum aggregate size")
  if(!varying) startInd <- NULL
  
  inp <- prepare_hsmm(dm, omega_into, allprobs, trackID, delta, first_obs = startInd, eps = eps)
  if(report) {
    report_hsmm(inp, dm, omega, allprobs, trackID)
    omega_out <- TRUE # omega[, , t] for the step out of t, needed for decoding
    REPORT(omega_out)
    if(varying && is.null(trackID)) report_sparse(inp, tpm_ihsmm(omega, dm, eps = eps))
  }
  
  run_forward_hsmm(inp)
}


#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} for hidden semi-Markov models with periodically inhomogeneous state durations and/ or conditional transition probabilities
#'
#' Calculates the (approximate) log-likelihood of a sequence of observations under a periodically inhomogeneous hidden semi-Markov model using a modified \strong{forward algorithm}.
#' 
#' @family forward algorithms
#' 
#' @description
#' Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs, where the state duration distribution is explicitly modelled by a distribution on the positive integers. This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary with covariates.
#' For direct numerical maximum likelihood estimation, HSMMs can be represented as HMMs on an enlarged state space (of size \eqn{M}) and with structured transition probabilities.
#' 
#' This function can be used to fit HSMMs where the state-duration distribution and/ or the conditional transition probabilities vary periodically.
#' In the special case of periodic variation (as compared to arbitrary covariate influence), this version is to be preferred over \code{\link{forward_ihsmm}} because it computes the \strong{correct periodically stationary distribution} and no observations are lost for the approximation.
#' 
#' This function is designed to be used with automatic differentiation based on the \code{R} package \code{RTMB}. It will be very slow without it!
#'
#' \strong{Deprecated:} please use \code{\link{forward_hsmm}} with \code{tod} instead. There, a stay that begins at time of day \code{s} is governed by \code{omega[,,s]}, here by \code{omega[,,s-1]} (\code{omega[,,s]} is for the step out of time of day \code{s}).
#'
#' @references Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state dwell-time distributions. Computational Statistics & Data Analysis, 209, 108171.
#'
#' @param dm list of length \code{nStates} containing matrices (or vectors) of dwell-time probability mass functions (PMFs) for each state.
#'
#' If the dwell-time PMFs are constant, the vectors are the PMF of the dwell-time distribution fixed in time. The vector lengths correspond to the approximating state aggregate sizes, hence there should be little probability mass not covered by these.
#' 
#' If the dwell-time PMFs are inhomogeneous, the matrices need to have L rows, where L is the cycle length. The number of columns again correspond to the size of the approximating state aggregates.
#' @param omega matrix of dimension \code{c(nStates, nStates)} or array of dimension \code{c(nStates, nStates, L)} of conditional transition probabilities, also called embedded transition probability matrix
#' 
#' It contains the transition probabilities given the current state is left. Hence, the diagonal elements need to be zero and the rows need to sum to one. Such a matrix can be constructed using \code{\link{tpm_emb}} and an array using \code{\link{tpm_emb_g}}.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension \code{c(nObs, nStates)}
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length \code{nObs}).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of length \code{nObs} containing IDs
#' 
#' If provided, the total log-likelihood will be the sum of each track's likelihood contribution.
#' Instead of a single vector \code{delta} corresponding to the initial distribution, a \code{delta} matrix of initial distributions, of dimension \code{c(nTracks, nStates)}, can be provided, such that each track starts with its own initial distribution.
#' @param delta optional vector of initial state probabilities of length \code{nStates}. By default, the periodically stationary distribution at the first time point is computed.
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
  
  deprecated(
    "forward_phsmm() is deprecated, please use forward_hsmm(..., tod = tod) instead, where a stay that begins at time of day s is governed by omega[,,s] (here: omega[,,s-1])."
  )
  
  # here, omega[,,s] is for the step out of time of day s, in forward_hsmm() for the step into s
  periodic_omega <- length(dim(omega)) == 3
  omega_into <- omega
  if(periodic_omega) {
    L <- dim(omega)[3]
    omega_into <- omega[, , c(L, seq_len(L - 1)), drop = FALSE]
  }
  
  inp <- prepare_hsmm(dm, omega_into, allprobs, trackID, delta, tod, eps = eps)
  # slice of the previous time of day, also if tod does not advance by one (possible with fixed dm)
  if(periodic_omega) inp$oidx <- c(1, tod[-length(tod)] %% L + 1)
  
  if(report) {
    report_hsmm(inp, dm, omega, allprobs, trackID, tod)
    omega_out <- TRUE # omega[, , s] for the step out of time of day s, needed for decoding
    REPORT(omega_out)
    if(is.matrix(dm[[1]]) && is.null(trackID)) report_sparse(inp, tpm_phsmm(omega, dm, eps = eps))
  }
  
  run_forward_hsmm(inp)
}
