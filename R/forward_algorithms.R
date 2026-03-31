# HMMs, ct-HMMs and MMPPs -------------------------------------------------

#' Forward algorithm to calculate the HMM log-likelihood
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model
#' using the \strong{forward algorithm} (Zucchini, MacDonald & Langrock, 2016).
#'
#' @references Zucchini, W., MacDonald, I.L., & Langrock, R. (2016). \emph{Hidden Markov Models for
#'   Time Series: An Introduction Using R} (2nd ed.). Chapman & Hall/CRC.
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
#' transition probability matrix; either
#'   \itemize{
#'     \item a matrix of dimension \code{c(nStates, nStates)}, or
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} if \code{trackID} is provided, or
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, in which case \code{\link{forward_g}} is called internally.
#'   }
#'   
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#' 
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks (see \sQuote{Details}).
#' 
#' @param logspace 
#' logical; if \code{TRUE}, \code{allprobs} is assumed to be on the log-scale, improving numerical stability for small probabilities. Only supported with \code{RTMB}.
#' 
#' @param bw 
#' optional positive integer specifying the bandwidth for a banded approximation of the forward algorithm, inducing a banded Hessian w.r.t. the observations. 
#' Defaults to \code{NULL} (exact algorithm). Approximation error decays geometrically in \code{bw}.
#' 
#' @param report 
#' logical; if \code{TRUE} (default), \code{delta}, \code{Gamma}, \code{allprobs}, and \code{trackID} are reported from the fitted model. Requires \code{ad = TRUE}.
#' 
#' @param ad 
#' logical; whether to use automatic differentiation. Determined automatically — for debugging only.
#'   
#' @details
#' If \code{trackID} is provided, the total log-likelihood is the sum of each track's likelihood
#' contribution. In this case, \code{Gamma} can be
#'   \itemize{
#'     \item a matrix (same transition probabilities for each track),
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} (track-specific transition probability matrix), or
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, in which case \code{\link{forward_g}} is called internally.
#'   }
#' Additionally, \code{delta} can be a vector (same initial distribution for each track) or a
#' matrix of dimension \code{c(nTracks, nStates)} (track-specific initial distributions).
#' 
#' \strong{Note:} When there are multiple tracks, for compatibility with downstream functions like \code{\link{viterbi}}, \code{\link{stateprobs}} or \code{\link{pseudo_res}}, 
#' \code{forward} should only be called \strong{once} with a \code{trackID} argument.
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
                    ad = NULL
                    ){
  
  # Check allprobs 
  if(!is.matrix(allprobs)) {
    stop("allprobs needs to be a matrix of dimension c(n, N).")
  } 
  nObs <- nrow(allprobs)
  nStates <- ncol(allprobs)
  
  # Handle time-varying Gamma by dispatching to forward_g
  if (is.array(Gamma) && length(dim(Gamma)) == 3 && (dim(Gamma)[3] == nObs | dim(Gamma)[3] == nObs-1)) {
    return(forward_g(delta, Gamma, allprobs, trackID = trackID,
                     logspace = logspace, bw = bw, report = report, ad = ad))
  }
  
  # exit if only one state
  if (nStates == 1) return(sum(log(allprobs)))
  
  # check trackID
  if (!is.null(trackID)) {
    if (length(trackID) != nObs) {
      stop("Length of trackID must equal number of rows in allprobs.")
    }
    # coerce to integer track indices
    trackID <- as.integer(as.factor(trackID))
    uID <- unique(trackID)
    nTracks <- length(uID)
    trackID <- match(trackID, uID) # ensure trackID is 1, 2, ..., nTracks
  } else {
    uID <- 1L
    nTracks <- 1L
  }
  
  # if allprobs on log-scale, report exp(allprobs)
  if(logspace){
    logallprobs <- allprobs
    allprobs <- exp(logallprobs)
  }
  
  # report quantities for easy use later
  if(report) {
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
    if(!is.null(trackID)) REPORT(trackID)
    type <- "homogeneous"
    REPORT(type)
    
    # if (length(.obs_info$calls) > 0) {
    #   obs_info <- .obs_info$calls
    #   REPORT(obs_info)
    #   rm(list = ls(.obs_info), envir = .obs_info)
    # }
  }
  
  # if ad is not explicitly provided, check via RTMB:::ad_context()
  if(is.null(ad)){
    ad <- ad_context()
  }
  
  
  ##### Dimension matching for input arguments #####
  
  ## handle delta
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

  # handle Gamma
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
  
  # expand across tracks if necessary
  if (is.null(trackID)) {
    if (dG[3] != 1) {
      stop("Gamma must have third dimension 1 when no trackID is provided.")
    }
  } else {
    if (dG[3] == 1 && nTracks > 1) {
      Gamma <- AD(array(Gamma, dim = c(nStates, nStates, nTracks)))
    } else if (dG[3] != nTracks) {
      stop("Gamma must have third dimension equal to number of tracks.")
    }
  }
  
  ### Inputs are cleaned up now: 
  # - now delta is matrix of dimension c(nTracks, nStates)
  # - and Gamma is array of dimension c(nStates, nStates, nTracks)

  ### Now run forward algorithm: 
  # - if ad == FALSE, use C++ version for speed
  # - if ad == TRUE, use R version which is converted into diff'able C++ code via RTMB
  
  if(!ad) { # non-AD version in C++
    
    # if(logspace) {
    #   stop("Logspace computations are not supported in the non-AD version.")
    # }
    
    if(is.null(trackID)) {
      l <- forward_cpp_h(allprobs, delta[1, ], Gamma[,, 1])
    } else {
      l <- forward_cpp_h_tracks(allprobs, delta, Gamma, calc_trackInd(trackID))
    }
    
  } else if(ad) { # AD version -> R code
    
    # printing suggestion to change TapeConfig
    cfg <- RTMB::TapeConfig()
    if(cfg$matmul != "plain" & nStates <= 5) {
      message("Performance tip: Consider running TapeConfig(matmul = 'plain') before MakeADFun() to speed up the forward algorithm.\n")
    }
    if(cfg$matmul == "plain" & nStates > 50) {
      msg <- paste0("Your model has a lot of states (", nStates, "). ",
                    "Run TapeConfig(matmul = 'atomic') or TapeConfig(matmul = 'compact') before MakeADFun() to speed up the forward algorithm.\n")
      message(msg)
    }
    
    # AD overloading to avoid trouble 
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    # Initialising log-likelihood at 0
    l <- 0 
    
    if(is.null(trackID)) {
      trackID <- rep(1, nObs) # if no trackID, only one track
    }
    # outer loop only has one iteration in this case
    
    ### forward algorithm code:
    # - if logspace == FALSE, regular computations
    # - if logspace == TRUE, computations on log-scale for numerical precision
    if(!logspace) { # regular, non-logspace computations
      for(i in seq_len(nTracks)) {
        ind <- which(trackID == uID[i]) # indices of track i
        
        allprobs_i <- allprobs[ind, , drop = FALSE] # matrix
        delta_i <- delta[i, , drop = FALSE] # matrix 
        Gamma_i <- Gamma[, , i] # matrix
        
        # Initialize forward algorithm
        foo <- delta_i * allprobs_i[1, , drop = FALSE]
        sumfoo <- sum(foo)
        phi <- foo / sumfoo
        l <- l + log(sumfoo)
        
        if(is.null(bw)) {
          k <- length(ind)
        } else {
          k <- bw
        }
        
        # Regular forward algorithm
        # If banded, only runs first block
        for(t in 2:min(k, length(ind))) {
          foo <- (phi %*% Gamma_i) * allprobs_i[t, , drop = FALSE]
          sumfoo <- sum(foo)
          phi <- foo / sumfoo
          l <- l + log(sumfoo)
        }
        
        # If banded, the approximation is as follows:
        # - compute state distribution of the chain at each time point, independent of the observations
        # - separate the log likelihood into blocks of size bw
        # - for each block, starting at time t
        #   - initialise with statedist_{t-bw}
        #   - use observations from t-bw to t-1 to update (scaled) forward variable to get close approx to true one
        #   - compute likelihood of block from t to t + bw using this approximation
        #   - this way, blocks are independent after a time lag of 2 * bw (at most)
        
        # If banded, run this banded version below for remaining blocks
        if(!is.null(bw) & (k < length(ind))) {
          # State distribution needs to be computed to initialise blocks
          stateDist <- AD(matrix(NaN, nrow = length(ind), ncol = nStates))
          stateDist[1, ] <- as.numeric(delta_i)
          for(l in 2:length(ind)) {
            stateDist[l, ] <- stateDist[l-1, ] %*% Gamma_i
          }
          
          startInd <- k + 1
          endInd <- 2 * k
          nBlocks <- ceiling(length(ind) / k) - 1 # number of remaining blocks
        
          for(b in seq_len(nBlocks)) {
            # updating forward variable to get good approximation
            foo <- stateDist[startInd-k, , drop = FALSE] * allprobs_i[startInd-k, , drop = FALSE]
            foo <- foo / sum(foo)
            for(t in (startInd-k+1):(startInd - 1)) { 
              foo <- (foo %*% Gamma_i) * allprobs_i[t, , drop = FALSE]
              foo <- foo / sum(foo)
            }
            # computing likelihood of this block based on approx forward variable
            end <- min(endInd, length(ind))
            for(t in startInd:end) { 
              foo <- (foo %*% Gamma_i) * allprobs_i[t, , drop = FALSE]
              sumfoo <- sum(foo)
              foo <- foo / sumfoo
              l <- l + log(sumfoo)
            }
            startInd <- startInd + k
            endInd <- endInd + k
          }
        }
      }
      
    } else { # logspace computations
      # when state-dependent probabilities are small, this is advantageous for numerical precision
      
      # to avoid log(0)
      delta <- delta + 1e-7
      delta <- delta / rowSums(delta)
      logdelta <- log(delta)
      
      for(i in seq_len(nTracks)) {
        ind <- which(trackID == uID[i]) # indices of track i
        
        logallprobs_i <- logallprobs[ind, , drop = FALSE] # matrix
        logdelta_i <- logdelta[i, , drop = FALSE] # matrix 
        Gamma_i <- Gamma[, , i] # matrix
        
        # Initialize forward algorithm
        logfoo <- logdelta_i + logallprobs_i[1, , drop = FALSE]
        logsumfoo <- logspace_add(logfoo)
        l <- l + logsumfoo
        logfoo <- logfoo - logsumfoo

        if(is.null(bw)) {
          k <- length(ind)
        } else {
          k <- bw
        }
        
        # Regular forward algorithm
        # If banded, only runs first block
        for(t in 2:min(k, length(ind))) {
          logfoo <- log((exp(logfoo) %*% Gamma_i)) + logallprobs_i[t, , drop = FALSE]
          logsumfoo <- logspace_add(logfoo)
          l <- l + logsumfoo
          logfoo <- logfoo - logsumfoo
        }
        
        # For explanation of the banded approximation, see above (non-logspace version)
        
        # If banded, run this baned version below for remaining blocks
        if(!is.null(bw) & (k < length(ind))) {
          # state distribution needs to be computed to initialise blocks
          stateDist <- AD(matrix(NaN, nrow = length(ind), ncol = nStates))
          stateDist[1, ] <- as.numeric(exp(logdelta_i))
          for(i in 2:length(ind)) {
            stateDist[i, ] <- stateDist[i-1, ] %*% Gamma_i
          }
          
          startInd <- k + 1
          endInd <- 2 * k
          nBlocks <- ceiling(length(ind) / k) - 1 # number of remaining blocks
          # Compute likelihood contribution of each block, independent from all other blocks
          for(b in seq_len(nBlocks)) {
            # updating forward variable to get good approximation
            logfoo <- log(stateDist[startInd-k, , drop = FALSE]) + logallprobs_i[startInd-k, , drop = FALSE]
            logsumfoo <- logspace_add(logfoo)
            logfoo <- logfoo - logsumfoo
            for(t in (startInd-k+1):(startInd - 1)) { 
              logfoo <- log(exp(logfoo) %*% Gamma_i) + logallprobs_i[t, , drop = FALSE]
              logsumfoo <- logspace_add(logfoo)
              logfoo <- logfoo - logsumfoo
            }
            # computing likelihood of this block based on approx forward variable
            end <- min(endInd, length(ind))
            for(t in startInd:end) { 
              logfoo <- log(exp(logfoo) %*% Gamma_i) + logallprobs_i[t, , drop = FALSE]
              logsumfoo <- logspace_add(logfoo)
              l <- l + logsumfoo
              logfoo <- logfoo - logsumfoo
            }
            startInd <- startInd + k
            endInd <- endInd + k
          }
        }
      }
    }
  }
  as.numeric(l)
}


#' Forward algorithm with time-varying transition probability matrix
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model with time-varying transition probabilities using the \strong{forward algorithm} (Zucchini, MacDonald & Langrock, 2016).
#'
#' @references Zucchini, W., MacDonald, I.L., & Langrock, R. (2016). \emph{Hidden Markov Models for
#'   Time Series: An Introduction Using R} (2nd ed.). Chapman & Hall/CRC.
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
#' 
#' @param logspace 
#' logical; if \code{TRUE}, \code{allprobs} is assumed to be on the log-scale, improving numerical stability for small probabilities. Only supported with \code{RTMB}.
#' 
#' @param bw 
#' optional positive integer specifying the bandwidth for a banded approximation of the forward algorithm, inducing a banded Hessian w.r.t. the observations. 
#' Defaults to \code{NULL} (exact algorithm). Approximation error decays geometrically in \code{bw}.
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
  
  # Check allprobs 
  if(!is.matrix(allprobs)) {
    stop("allprobs needs to be a matrix of dimension c(n, N).")
  } 
  nObs <- nrow(allprobs)
  nStates <- ncol(allprobs)
  
  # exit if only one state
  if (nStates == 1) return(sum(log(allprobs)))
  
  # check trackID
  if (!is.null(trackID)) {
    if (length(trackID) != nObs) {
      stop("Length of trackID must equal number of rows in allprobs.")
    }
    # coerce to integer track indices
    trackID <- as.integer(as.factor(trackID))
    uID <- unique(trackID)
    nTracks <- length(uID)
    trackID <- match(trackID, uID) # ensure trackID is 1, 2, ..., nTracks
  } else {
    uID <- 1L
    nTracks <- 1L
  }
  
  # if allprobs on log-scale, report exp(allprobs)
  if(logspace){
    logallprobs <- allprobs
    allprobs <- exp(logallprobs)
  }
  
  # report quantities for easy use later
  if(report) {
    REPORT(delta)
    REPORT(Gamma)
    REPORT(allprobs)
    if(!is.null(trackID)) REPORT(trackID)
    # if(!is.null(attributes(Gamma)$time) && attributes(Gamma)$time == "continuous") {
    #   type <- "continuous_time"
    # } else {
    #   type <- "inhomogeneous"
    # }
    type <- "inhomogeneous"
    REPORT(type)
    
    # if (length(.obs_info$calls) > 0) {
    #   obs_info <- .obs_info$calls
    #   REPORT(obs_info)
    #   rm(list = ls(.obs_info), envir = .obs_info)
    # }
  }
  
  # if ad is not explicitly provided, check via RTMB:::ad_context()
  if(is.null(ad)){
    ad <- ad_context()
  }
  
  
  ##### Dimension matching for input arguments #####
  
  ## handle delta
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
  
  # handle Gamma
  
  Gamma <- AD(as.array(Gamma))
  dG <- dim(Gamma)
  
  # if matrix input, break and suggest using forward()
  if (length(dG) == 2) {
    stop("For time-constant transition probability matrix, please use forward().")
  }
  
  # basic square check
  if (!all(dG[1:2] == nStates)) {
    stop("Each Gamma slice must be square with dimensions (N, N).")
  }
  
  nGamma <- dG[3]
  
  if (is.null(trackID)) {
    # single track case: allow nObs - 1
    if (nGamma == nObs - 1) {
      # prepend an unused dummy slice (e.g. identity transition)
      Gamma_pad <- AD(array(NaN, dim = c(nStates, nStates, nObs)))
      Gamma_pad[,,2:nObs] <- Gamma
      Gamma <- Gamma_pad
    } else if (nGamma != nObs) {
      stop("For single-track data, Gamma must have third dimension nObs or nObs - 1.")
    }
  } else {
    # multi-track case: must match total number of observations
    if (nGamma != nObs) {
      stop("For multi-track data, Gamma must have third dimension equal to nObs.")
    }
  }
  
  ### Inputs are cleaned up now: 
  # - now delta is matrix of dimension c(nTracks, nStates)
  # - and Gamma is array of dimension c(nStates, nStates, nObs)
  
  ### Now run forward algorithm: 
  # - if ad == FALSE, use C++ version for speed
  # - if ad == TRUE, use R version which is converted into diff'able C++ code via RTMB
  
  if(!ad) { # non-AD version in C++
    
    # if(logspace) {
    #   stop("Logspace computations are not supported in the non-AD version.")
    # }
    
    if(is.null(trackID)) {
      l <- forward_cpp_g(allprobs, delta[1, ], Gamma[,,-1])
    } else {
      l <- forward_cpp_g_tracks(allprobs, delta, Gamma, calc_trackInd(trackID))
    }
    
  } else if(ad) { # AD version -> R code
    
    # printing suggestion to change TapeConfig
    cfg <- RTMB::TapeConfig()
    if(cfg$matmul != "plain" & nStates <= 5) {
      message("Performance tip: Consider running TapeConfig(matmul = 'plain') before MakeADFun() to speed up the forward algorithm.\n")
    }
    if(cfg$matmul == "plain" & nStates > 50) {
      msg <- paste0("Your model has a lot of states (", nStates, "). ",
                    "Run TapeConfig(matmul = 'atomic') or TapeConfig(matmul = 'compact') before MakeADFun() to speed up the forward algorithm.\n")
      message(msg)
    }
    
    # AD overloading to avoid trouble 
    "[<-" <- ADoverload("[<-")
    "c" <- ADoverload("c")
    "diag<-" <- ADoverload("diag<-")
    
    # Initialising log-likelihood at 0
    l <- 0 
    
    if(is.null(trackID)) {
      trackID <- rep(1, nObs) # if no trackID, only one track
    }
    # outer loop only has one iteration in this case
    
    ### forward algorithm code:
    # - if logspace == FALSE, regular computations
    # - if logspace == TRUE, computations on log-scale for numerical precision
    if(!logspace) { # regular, non-logspace computations
      for(i in seq_len(nTracks)) {
        ind <- which(trackID == uID[i]) # indices of track i
        
        allprobs_i <- allprobs[ind, , drop = FALSE] # matrix
        delta_i <- delta[i, , drop = FALSE] # matrix 
        Gamma_i <- Gamma[, , ind, drop = FALSE] # select relevant slices
        
        # Initialize forward algorithm
        foo <- delta_i * allprobs_i[1, , drop = FALSE]
        sumfoo <- sum(foo)
        phi <- foo / sumfoo
        l <- l + log(sumfoo)
        
        if(is.null(bw)) {
          k <- length(ind)
        } else {
          k <- bw
        }
        
        # Regular forward algorithm
        # If banded, only runs first block
        for(t in 2:min(k, length(ind))) { # starts at t, so first slice of Gamma is never used
          foo <- (phi %*% Gamma_i[,,t]) * allprobs_i[t, , drop = FALSE]
          sumfoo <- sum(foo)
          phi <- foo / sumfoo
          l <- l + log(sumfoo)
        }
        
        # If banded, the approximation is as follows:
        # - compute state distribution of the chain at each time point, independent of the observations
        # - separate the log likelihood into blocks of size bw
        # - for each block, starting at time t
        #   - initialise with statedist_{t-bw}
        #   - use observations from t-bw to t-1 to update (scaled) forward variable to get close approx to true one
        #   - compute likelihood of block from t to t + bw using this approximation
        #   - this way, blocks are independent after a time lag of 2 * bw (at most)
        
        # If banded, run this banded version below for remaining blocks
        if(!is.null(bw) & (k < length(ind))) {
          # State distribution needs to be computed to initialise blocks
          # stateDist <- stateDist_banded(delta_i, Gamma_i, bw = bw)
          stateDist <- rep(1, nStates) / nStates
          
          startInd <- k + 1
          endInd <- 2 * k
          nBlocks <- ceiling(length(ind) / k) - 1 # number of remaining blocks
          
          for(b in seq_len(nBlocks)) {
            # updating forward variable to get good approximation
            # foo <- stateDist[startInd-k, , drop = FALSE] * allprobs_i[startInd-k, , drop = FALSE]
            foo <- stateDist * allprobs_i[startInd-k, , drop = FALSE]
            foo <- foo / sum(foo)
            for(t in (startInd-k+1):(startInd - 1)) { 
              foo <- (foo %*% Gamma_i[, , t]) * allprobs_i[t, , drop = FALSE]
              foo <- foo / sum(foo)
            }
            # computing likelihood of this block based on approx forward variable
            end <- min(endInd, length(ind))
            for(t in startInd:end) { 
              foo <- (foo %*% Gamma_i[, , t]) * allprobs_i[t, , drop = FALSE]
              sumfoo <- sum(foo)
              foo <- foo / sumfoo
              l <- l + log(sumfoo)
            }
            startInd <- startInd + k
            endInd <- endInd + k
          }
        }
      }
      
    } else { # logspace computations
      # when state-dependent probabilities are small, this is advantageous for numerical precision
      
      # to avoid log(0)
      delta <- delta + 1e-7
      delta <- delta / rowSums(delta)
      logdelta <- log(delta)
      
      for(i in seq_len(nTracks)) {
        ind <- which(trackID == uID[i]) # indices of track i
        
        logallprobs_i <- logallprobs[ind, , drop = FALSE] # matrix
        logdelta_i <- logdelta[i, , drop = FALSE] # matrix 
        Gamma_i <- Gamma[, , ind, drop = FALSE] # select relevant slices
        
        # Initialize forward algorithm
        logfoo <- logdelta_i + logallprobs_i[1, , drop = FALSE]
        logsumfoo <- logspace_add(logfoo)
        l <- l + logsumfoo
        logfoo <- logfoo - logsumfoo
        
        if(is.null(bw)) {
          k <- length(ind)
        } else {
          k <- bw
        }
        
        # Regular forward algorithm
        # If banded, only runs first block
        for(t in 2:min(k, length(ind))) {
          logfoo <- log((exp(logfoo) %*% Gamma_i[, , t])) + logallprobs_i[t, , drop = FALSE]
          logsumfoo <- logspace_add(logfoo)
          l <- l + logsumfoo
          logfoo <- logfoo - logsumfoo
        }
        
        # For explanation of the banded approximation, see above (non-logspace version)
        
        # If banded, run this baned version below for remaining blocks
        if(!is.null(bw) & (k < length(ind))) {
          # state distribution needs to be computed to initialise blocks
          # stateDist <- stateDist_banded(exp(logdelta_i), Gamma_i, bw = bw)
          stateDist <- rep(1, nStates) / nStates
          
          startInd <- k + 1
          endInd <- 2 * k
          nBlocks <- ceiling(length(ind) / k) - 1 # number of remaining blocks
          # Compute likelihood contribution of each block, independent from all other blocks
          for(b in seq_len(nBlocks)) {
            # updating forward variable to get good approximation
            # logfoo <- log(stateDist[startInd-k, , drop = FALSE]) + logallprobs_i[startInd-k, , drop = FALSE]
            logfoo <- log(stateDist) + logallprobs_i[startInd-k, , drop = FALSE]
            logsumfoo <- logspace_add(logfoo)
            logfoo <- logfoo - logsumfoo
            for(t in (startInd-k+1):(startInd - 1)) { 
              logfoo <- log(exp(logfoo) %*% Gamma_i[, , t]) + logallprobs_i[t, , drop = FALSE]
              logsumfoo <- logspace_add(logfoo)
              logfoo <- logfoo - logsumfoo
            }
            # computing likelihood of this block based on approx forward variable
            end <- min(endInd, length(ind))
            for(t in startInd:end) { 
              logfoo <- log(exp(logfoo) %*% Gamma_i[, , t]) + logallprobs_i[t, , drop = FALSE]
              logsumfoo <- logspace_add(logfoo)
              l <- l + logsumfoo
              logfoo <- logfoo - logsumfoo
            }
            startInd <- startInd + k
            endInd <- endInd + k
          }
        }
      }
    }
  }
  as.numeric(l)
}

# allows for computing state distribution of an inhomogeneous Markov chain in a banded mannor
# when Gamma has a parameter per slice, these are independent after lag 2 * bw
stateDist_banded <- function(delta, Gamma, bw){
  nObs <- dim(Gamma)[3] # number of observations
  nStates <- dim(Gamma)[1] # number of states
  nBlocks <- ceiling(nObs / bw) - 2 # number of remaining tracks
  rho <- rep(1, nStates) / nStates # fixed intermediate state distribution
  
  # Initialise Delta matrix to store state distributions
  Delta <- matrix(NaN, nrow = nObs, ncol = nStates)
  
  # Set initial state distribution
  Delta[1, ] <- as.numeric(delta)
  # Compute the state distribution at time t
  for(t in 2:(2*bw)){ 
    Delta[t, ] <- Delta[t - 1, ] %*% Gamma[, , t]
  }
  
  # remaining blocks, always initialise with rho to keep independent
  startInd <- bw + 1
  endInd <- bw + 2 * bw
  for(b in seq_len(nBlocks)){
    foo <- rho
    end <- min(endInd, nObs)
    for(t in startInd:end){
      foo <- foo %*% Gamma[, , t]
      if(t > startInd + bw - 1){
        Delta[t, ] <- foo
      }
    }
    startInd <- startInd + bw
    endInd <- endInd + bw
  }
  Delta
}



#' \href{https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock}{Forward algorithm} with for periodically varying transition probability matrices
#'
#' Calculates the log-likelihood of a sequence of observations under a hidden Markov model with periodically varying transition probabilities using the \strong{forward algorithm}.
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
#'  delta = stationary_p(Gamma, tod[1]) # stationary HMM
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
  
  utod = unique(tod) # unique time of days
  L = length(utod) # cycle length
  
  if(dim(Gamma)[3] != L){
    stop("Gamma needs to be an array of dimension c(N,N,L), matching the number of unique time points in tod.") 
  } 
  
  Gammanew = Gamma[,,tod]
  
  if(report){
    REPORT(tod)
    type <- "periodic"
    REPORT(type)
  }
  
  forward_g(delta = delta, 
            Gamma = Gammanew, 
            allprobs = allprobs, 
            trackID = trackID, 
            ad = ad, 
            report = report,
            logspace = logspace)
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
