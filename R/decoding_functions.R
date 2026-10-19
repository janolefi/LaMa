

# Takes delta, Gamma, allprobs, trackID and tod from mod (if provided) and
# brings them into the form of the forward algorithm, see prepare_forward()
prepare_decoding <- function(delta, Gamma, allprobs, trackID = NULL, tod = NULL, mod = NULL){
  if(!is.null(mod)){
    for(x in c("delta", "Gamma", "allprobs")) {
      if(is.null(mod[[x]])) stop("Model object contains no ", x, ".")
    }
    delta <- mod$delta
    Gamma <- mod$Gamma
    allprobs <- mod$allprobs
    trackID <- mod$trackID
    # periodic models (LaMa <= 2.1.1 reported tod with an expanded Gamma, type inhomogeneous)
    if(!identical(mod$type, "inhomogeneous")) tod <- mod$tod
  }
  prepare_forward(delta, Gamma, allprobs, trackID, gidx = tod)
}

# Viterbi algorithm in C++ on the output of prepare_decoding()
run_viterbi <- function(inp){
  as.integer(viterbi_cpp(inp$allprobs, inp$delta, inp$Gamma, inp$gidx - 1L, inp$start - 1L))
}


# Global decoding ---------------------------------------------------------


#' Viterbi algorithm for state decoding in HMMs
#'
#' The Viterbi algorithm decodes the most probable state sequence of an HMM.
#'
#' @family decoding functions
#'
#' @param delta 
#' initial distribution; either
#'   \itemize{
#'     \item a vector of length \code{nStates}, or
#'     \item a matrix of dimension \code{c(nTracks, nStates)} if \code{trackID} is provided
#'   }
#' @param Gamma 
#' transition probability matrix; either
#'   \itemize{
#'     \item a matrix of dimension \code{c(nStates, nStates)},
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} if \code{trackID} is provided, or
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, in which case \code{\link{viterbi_g}} is called internally
#'   }
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#'
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{RTMB::report()} or \code{\link{qreml}} can be passed directly.
#'
#' @return vector of decoded states of length \code{nObs}
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' states = viterbi(delta, Gamma, allprobs)
viterbi = function(delta, Gamma, allprobs, trackID = NULL,
                   mod = NULL){
  inp <- prepare_decoding(delta, Gamma, allprobs, trackID, mod = mod)
  run_viterbi(inp)
}


#' Viterbi algorithm for state decoding in inhomogeneous HMMs
#'
#' The Viterbi algorithm decodes the most probable state sequence of an HMM.
#'
#' @family decoding functions
#'
#' @param delta 
#' initial distribution; either
#'   \itemize{
#'     \item a vector of length \code{nStates}, or
#'     \item a matrix of dimension \code{c(nTracks, nStates)} if \code{trackID} is provided
#'   }
#'
#' @param Gamma 
#' array of transition probability matrices of dimension \code{c(nStates, nStates, nObs)},
#'   where the first slice of each track is ignored as there is no transition into the start of a track.
#'   For a single track, an array of dimension \code{c(nStates, nStates, nObs-1)} is also accepted.
#'
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#'
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward_g}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{RTMB::report()} or \code{\link{qreml}} can be passed directly.
#'
#' @return vector of decoded states of length \code{nObs}
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' Gamma = tpm_g(runif(10), matrix(c(-2,-2,1,-1), nrow = 2))
#' allprobs = matrix(runif(20), nrow = 10, ncol = 2)
#' states = viterbi_g(delta, Gamma[,,-1], allprobs)
viterbi_g = function(delta, Gamma, allprobs, trackID = NULL,
                     mod = NULL){
  inp <- prepare_decoding(delta, Gamma, allprobs, trackID, mod = mod)
  run_viterbi(inp)
}


#' Viterbi algorithm for state decoding in periodically inhomogeneous HMMs
#' 
#' The Viterbi algorithm allows one to decode the most probable state sequence of an HMM.
#'
#' @family decoding functions
#'
#' @param delta initial distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#' 
#' This could e.g. be the periodically stationary distribution (for each track).
#' @param Gamma array of transition probability matrices for each time point in the cycle of dimension c(N,N,L), where L is the length of the cycle
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n)
#' 
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#'
#' @return vector of decoded states of length n
#' @export
#'
#' @examples
#' delta = c(0.5, 0.5)
#' beta = matrix(c(-2, 1, -1,
#'                 -2, -1, 1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:24, 24, beta)
#' 
#' tod = rep(1:24, 5)
#' n = length(tod)
#' 
#' allprobs = matrix(runif(2*n), nrow = n, ncol = 2)
#' states = viterbi_p(delta, Gamma, allprobs, tod)
viterbi_p = function(delta, Gamma, allprobs, tod, trackID = NULL,
                     mod = NULL){
  inp <- prepare_decoding(delta, Gamma, allprobs, trackID, tod, mod)
  run_viterbi(inp)
}



# Local decoding ----------------------------------------------------------


#' Calculate conditional local state probabilities in HMMs
#' 
#' Computes
#' \deqn{\Pr(\text{State}_t = j \mid X_1, ..., X_T)}
#' for a given HMM.
#'
#' @family decoding functions
#' 
#' @param delta 
#' initial distribution; either
#'   \itemize{
#'     \item a vector of length \code{nStates}, or
#'     \item a matrix of dimension \code{c(nTracks, nStates)} if \code{trackID} is provided
#'   }
#' @param Gamma 
#' transition probability matrix; either
#'   \itemize{
#'     \item a matrix of dimension \code{c(nStates, nStates)},
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} if \code{trackID} is provided, or
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, in which case \code{\link{stateprobs_g}} is called internally
#'   }
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#'
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{RTMB::report()} or \code{\link{qreml}} can be passed directly.
#'
#' @param forecast 
#' logical, indicating if forecast probabilities \eqn{\Pr(\text{State}_t = j \mid X_1, ..., X_{t-1})} should be calculated instead.
#' 
#' @return matrix of conditional state probabilities of dimension \code{c(nObs, nStates)}
#' @export
#'
#' @examples
#' Gamma = tpm(c(-1,-2))
#' delta = stationary(Gamma)
#' allprobs = matrix(runif(10), nrow = 10, ncol = 2)
#' 
#' probs = stateprobs(delta, Gamma, allprobs)
stateprobs = function(delta, Gamma, allprobs, trackID = NULL,
                      mod = NULL, forecast = FALSE) {
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }

    # dispatch on the reported model type
    if(identical(mod$type, "periodic")) return(stateprobs_p(mod = mod, forecast = forecast))
    if(identical(mod$type, "inhomogeneous")) return(stateprobs_g(mod = mod, forecast = forecast))

    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    allprobs = mod$allprobs
    trackID = mod$trackID
  }

  n = nrow(allprobs)
  N = ncol(allprobs)

  # Handle time-varying Gamma by dispatching to stateprobs_g (if model type is unknown)
  if (is.null(mod$type) && is.array(Gamma) && length(dim(Gamma)) == 3 && (dim(Gamma)[3] == n | dim(Gamma)[3] == n-1)) {
    return(stateprobs_g(delta, Gamma, allprobs, trackID = trackID, forecast = forecast))
  }
  
  # inflating Gamma to use stateprobs_g
  if(is.null(trackID)){
    Gammanew = array(Gamma, dim = c(N, N, n-1))
  } else{
    uID = unique(trackID)
    k = length(uID) # number of tracks
    
    if(is.matrix(Gamma)){
      Gammanew = array(Gamma, dim = c(N, N, n))
    } else if(is.array(Gamma)){
      if(dim(Gamma)[3] != k) stop("Number of distinct transition matrices does not match the number of tracks.")
      
      ## construct integer trackID
      integerID = match(trackID, uID)
      
      Gammanew = Gamma[,,integerID]
      Gammanew = Gammanew[,, -1] # ignoring first slice
    }
  }
  
  stateprobs_g(delta, Gammanew, allprobs, trackID, forecast = forecast)
}


#' Calculate conditional local state probabilities for inhomogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(\text{State}_t = j \mid X_1, ..., X_T)}
#' for inhomogeneous HMMs
#' 
#' @family decoding functions
#' 
#' @param delta 
#' initial distribution; either
#'   \itemize{
#'     \item a vector of length \code{nStates}, or
#'     \item a matrix of dimension \code{c(nTracks, nStates)} if \code{trackID} is provided
#'   }
#'
#' @param Gamma 
#' array of transition probability matrices of dimension \code{c(nStates, nStates, nObs)},
#'   where the first slice of each track is ignored as there is no transition into the start of a track.
#'   For a single track, an array of dimension \code{c(nStates, nStates, nObs-1)} is also accepted.
#'
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#'
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward_g}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{RTMB::report()} or \code{\link{qreml}} can be passed directly.
#'
#' @param forecast 
#' logical, indicating if forecast probabilities \eqn{\Pr(\text{State}_t = j \mid X_1, ..., X_{t-1})} should be calculated instead.
#'
#' @return matrix of conditional state probabilities of dimension \code{c(nObs, nStates)}
#' @export
#'
#' @examples
#' Gamma = tpm_g(runif(10), matrix(c(-1,-1,1,-2), nrow = 2, byrow = TRUE))
#' delta = c(0.5, 0.5)
#' allprobs = matrix(runif(20), nrow = 10, ncol = 2)
#' 
#' probs = stateprobs_g(delta, Gamma[,,-1], allprobs)
stateprobs_g = function(delta, Gamma, allprobs, trackID = NULL,
                        mod = NULL, forecast = FALSE) {
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    # Gamma = mod$Gamma[,, -1] # ignoring first slice
    allprobs = mod$allprobs
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  # initialising state vector
  stateprobs = matrix(NA, nrow = n, ncol = N)
  
  # If ID is provided, several tracks need to be decoded separately
  if(!is.null(trackID)){
    uID = unique(trackID) # unique IDs
    k = length(uID) # number of tracks
    
    if(is.vector(delta)){ # single initial distribution
      delta = matrix(delta, nrow = k, ncol = length(delta), byrow = TRUE)
    } else if(is.matrix(delta)){
      if(dim(delta)[1] != k){
        if(dim(delta)[1] == 1){
          delta = matrix(c(delta), nrow = k, ncol = length(delta), byrow = TRUE)
        } else{
          stop("Delta needs to be either a vector of length N or a matrix of dimension c(k,N), matching the number tracks.")
        }
      }
    }
    
    if(dim(delta)[2] != N) stop("Initial distribution(s) do not match the number of states implied by allprobs.")
    
    # loop over individual tracks
    for(i in seq_len(length(uID))){
      id_i = which(trackID == uID[i])
      allprobs_i = allprobs[id_i, ]
      
      if(length(id_i) == 1) stop("All tracks must be at least of length 2.")
      
      if(nrow(allprobs_i) == 2){
        # viterbi algorithm for track of length 2 only
        Gamma_i = as.matrix(Gamma[,,id_i[-1]])
        
        lalpha = matrix(NA, 2, N)
        lbeta = matrix(NA, 2, N)
        
        # computing forward variables on log scale
        foo = delta[i,] * allprobs_i[1,]
        l = log(sum(foo))
        foo = foo / sum(foo)
        lalpha[1,] = log(foo) + l
        
        foo = (foo %*% Gamma_i) * allprobs_i[2,]
        l = l + log(sum(foo))
        foo = foo / sum(foo)
        lalpha[2,] = log(foo) + l
        
        if(!forecast) {
          # computing backward variables on log scale
          foo = rep(1, N)
          lbeta[2,] = rep(0, N)
          foo = Gamma_i %*% diag(allprobs_i[2, ]) %*% foo
          lbeta[1,] = log(foo)
          
          c = max(lalpha[2, ])
          llk = c + log(sum(exp(lalpha[2, ] - c)))
          stateprobs[id_i, ] = exp(lalpha + lbeta - llk)
        } else {
          phi1 <- exp(lalpha[1,])
          phi2 <- as.vector(phi1 %*% Gamma_i)
          stateprobs[id_i, ] <- rbind(delta[i,], phi2) # forecast at the start is delta
        }
   
      } else{
        # regurlar local decoding
        lalpha = logalpha_cpp(allprobs_i, delta[i,], Gamma[,,id_i[-1]])
        
        if(!forecast) {
          lbeta = logbeta_cpp(allprobs_i, Gamma[,,id_i[-1]])
          c = max(lalpha[nrow(lalpha), ])
          llk = c + log(sum(exp(lalpha[nrow(lalpha), ] - c)))
          stateprobs[id_i, ] = exp(lalpha + lbeta - llk)
        } else {
          rowmax <- apply(lalpha, 1, max)
          Alpha <- exp(lalpha - rowmax)
          Phi <- Alpha / rowSums(Alpha)
          stateprobs[id_i[1], ] <- delta[i, ] # forecast at the start is delta
          for(t in seq_along(id_i)[-1]) {
            stateprobs[id_i[t], ] <- Phi[t-1, ] %*% Gamma[,,id_i[t]]
          }
        }
      }
    }
    
  } else{
    if(!is.vector(delta)){
      stop("If trackID is not provided, delta needs to be a vector of length N.")
    }
    
    if(dim(Gamma)[3] == n){
      # message("Igoring the first slice of Gamma, as there are only n-1 transitions in a time series of length n.")
      # not using the first slice of Gamma, if n slices are provided
      Gamma = Gamma[,,-1]
    }
    
    lalpha = logalpha_cpp(allprobs, delta, Gamma)
    
    if(!forecast) {
      lbeta = logbeta_cpp(allprobs, Gamma)
      
      c = max(lalpha[nrow(lalpha),])
      llk = c + log(sum(exp(lalpha[nrow(lalpha),] - c)))
      
      stateprobs = exp(lalpha + lbeta - llk)
    } else{
      rowmax <- apply(lalpha, 1, max)
      Alpha <- exp(lalpha - rowmax)
      Phi <- Alpha / rowSums(Alpha)
      stateprobs[1, ] <- delta # forecast at the start is delta
      for(t in 2:n) {
        stateprobs[t, ] <- Phi[t-1, ] %*% Gamma[,,t-1]
      }
    }
  }
  
  # rowSums should already be one, but just to be safe
  stateprobs <- stateprobs / rowSums(stateprobs)
  
  colnames(stateprobs) <- paste0("S", 1:N)
  stateprobs
}


#' Calculate conditional local state probabilities for periodically inhomogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(S_t = j \mid X_1, ..., X_T)}
#' for periodically inhomogeneous HMMs
#'  
#' @family decoding functions
#' 
#' @param delta initial or stationary distribution of length N, or matrix of dimension c(k,N) for k independent tracks, if \code{trackID} is provided
#'
#' This could e.g. be the periodically stationary distribution (for each track) as computed by \code{\link{stationary_p}}.
#' @param Gamma array of transition probability matrices for each time point in the cycle of dimension c(N,N,L), where L is the length of the cycle.
#' @param allprobs matrix of state-dependent probabilities/ density values of dimension c(n, N)
#' @param tod (Integer valued) variable for cycle indexing in 1, ..., L, mapping the data index to a generalised time of day (length n).
#' For half-hourly data L = 48. It could, however, also be day of year for daily data and L = 365.
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from running \code{RTMB::report()} or from \code{\link{qreml}} directly to this function.
#' @param forecast logical, indicating if forecast probabilities \eqn{\Pr(S_t = j \mid X_1, ..., X_{t-1})} should be calculated instead.
#'
#' @return matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:L, L, beta, degree = 1)
#' delta = stationary_p(Gamma, 1)
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' tod = rep(1:24, 5)[1:100]
#' 
#' probs = stateprobs_p(delta, Gamma, allprobs, tod)

stateprobs_p = function(delta, Gamma, allprobs, tod, trackID = NULL,
                        mod = NULL, forecast = FALSE) {
  
  # check if a model with delta, Gamma and allprobs is provided
  if(!is.null(mod)){
    if(is.null(mod$delta)){
      stop("Model object contains no initial distribution.")
    }
    if(is.null(mod$Gamma)){
      stop("Model object contains no transition matrix.")
    }
    if(is.null(mod$allprobs)){
      stop("Model object contains no state-dependent probabilities.")
    }
    if(is.null(mod$tod)){
      stop("Model object contains no cyclic indexing variable.")
    }
    
    # if suitable model object is provided, overwrite inputs with model object
    delta = mod$delta
    Gamma = mod$Gamma
    allprobs = mod$allprobs
    tod = mod$tod
    trackID = mod$trackID
  }
  
  n = nrow(allprobs)
  N = ncol(allprobs)
  
  Gammanew = Gamma[,,tod] # select the transition matrix for the current time of day
  
  if(is.null(trackID)){
    Gammanew = Gammanew[,,-1]
  }
  
  stateprobs_g(delta, Gammanew, allprobs, trackID, forecast = forecast)
}

