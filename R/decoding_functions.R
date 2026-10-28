

# Takes delta, Gamma, allprobs, trackID and tod from mod (if provided) and
# brings them into the form of the forward algorithm, see prepare_forward()
prepare_decoding <- function(delta, Gamma, allprobs, trackID = NULL, tod = NULL, mod = NULL){
  # model object as first argument, e.g. viterbi(mod)
  if(is.null(mod) && is.list(delta)) mod <- delta
  if(!is.null(mod)){
    for(x in c("delta", "Gamma", "allprobs")) {
      if(is.null(mod[[x]])) stop("Model object contains no ", x, ".")
    }
    delta <- mod$delta
    Gamma <- mod$Gamma
    allprobs <- mod$allprobs
    trackID <- mod$trackID
    # periodic models (LaMa <= 2.1.1 reported tod with an expanded Gamma, type inhomogeneous)
    tod <- if(!identical(mod$type, "inhomogeneous")) mod$tod
  }
  prepare_forward(delta, Gamma, allprobs, trackID, gidx = tod)
}

# Viterbi algorithm in C++ on the output of prepare_decoding()
run_viterbi <- function(inp){
  as.integer(viterbi_cpp(inp$allprobs, inp$delta, inp$Gamma, inp$gidx - 1L, inp$start - 1L))
}

# local state probabilities in C++ on the output of prepare_decoding()
run_stateprobs <- function(inp, forecast = FALSE){
  probs <- stateprobs_cpp(inp$allprobs, inp$delta, inp$Gamma, inp$gidx - 1L, inp$start - 1L, forecast)
  probs <- probs / rowSums(probs) # normalise rows
  colnames(probs) <- paste0("S", seq_len(ncol(probs)))
  probs
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
#'     \item a vector of length \code{nStates},
#'     \item a matrix of dimension \code{c(nTracks, nStates)} if \code{trackID} is provided, or
#'     \item a model object as for \code{mod}, so that \code{mod} can be passed as the first argument
#'   }
#' @param Gamma 
#' transition probability matrix; either
#'   \itemize{
#'     \item a matrix of dimension \code{c(nStates, nStates)},
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} if \code{trackID} is provided,
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, or
#'     \item an array of dimension \code{c(nStates, nStates, L)} for periodic transition probabilities, if \code{tod} is provided
#'   }
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#'
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks.
#' The observations of each track need to be consecutive.
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} can be passed directly.
#'
#' @param tod 
#' optional integer vector of length \code{nObs} with values in 1, ..., L, indexing the slices of \code{Gamma} for periodic models
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
                   mod = NULL, tod = NULL){
  inp <- prepare_decoding(delta, Gamma, allprobs, trackID, tod, mod)
  run_viterbi(inp)
}


#' Viterbi algorithm for state decoding in inhomogeneous HMMs
#'
#' The Viterbi algorithm decodes the most probable state sequence of an HMM.
#' 
#' \strong{Deprecated:} please use \code{\link{viterbi}} for all model types instead.
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
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks.
#' The observations of each track need to be consecutive.
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward_g}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} can be passed directly.
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
  
  deprecated(
    "viterbi_g() is deprecated, please use viterbi() for all model types instead."
  )
  
  viterbi(delta, Gamma, allprobs, trackID = trackID, mod = mod)
}


#' Viterbi algorithm for state decoding in periodically inhomogeneous HMMs
#' 
#' The Viterbi algorithm allows one to decode the most probable state sequence of an HMM.
#' 
#' \strong{Deprecated:} please use \code{\link{viterbi}} with argument \code{tod} instead.
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
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately. The observations of each track need to be consecutive.
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} directly to this function.
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
  
  deprecated(
    "viterbi_p() is deprecated, please use viterbi(..., tod = tod) for periodic models instead."
  )
  
  viterbi(delta, Gamma, allprobs, trackID = trackID, mod = mod, tod = tod)
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
#'     \item a vector of length \code{nStates},
#'     \item a matrix of dimension \code{c(nTracks, nStates)} if \code{trackID} is provided, or
#'     \item a model object as for \code{mod}, so that \code{mod} can be passed as the first argument
#'   }
#' @param Gamma 
#' transition probability matrix; either
#'   \itemize{
#'     \item a matrix of dimension \code{c(nStates, nStates)},
#'     \item an array of dimension \code{c(nStates, nStates, nTracks)} if \code{trackID} is provided,
#'     \item an array of dimension \code{c(nStates, nStates, nObs)} for time-varying transition probabilities, or
#'     \item an array of dimension \code{c(nStates, nStates, L)} for periodic transition probabilities, if \code{tod} is provided
#'   }
#' @param allprobs 
#' matrix of state-dependent probabilities or density values of dimension \code{c(nObs, nStates)}
#'
#' @param trackID 
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks.
#' The observations of each track need to be consecutive.
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} can be passed directly.
#'
#' @param forecast 
#' logical, indicating if forecast probabilities \eqn{\Pr(\text{State}_t = j \mid X_1, ..., X_{t-1})} should be calculated instead.
#' 
#' @param tod 
#' optional integer vector of length \code{nObs} with values in 1, ..., L, indexing the slices of \code{Gamma} for periodic models
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
                      mod = NULL, forecast = FALSE, tod = NULL) {
  inp <- prepare_decoding(delta, Gamma, allprobs, trackID, tod, mod)
  run_stateprobs(inp, forecast)
}


#' Calculate conditional local state probabilities for inhomogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(\text{State}_t = j \mid X_1, ..., X_T)}
#' for inhomogeneous HMMs
#' 
#' \strong{Deprecated:} please use \code{\link{stateprobs}} for all model types instead.
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
#' optional vector of length \code{nObs} containing \code{nTracks} unique IDs that separate tracks.
#' The observations of each track need to be consecutive.
#'
#' @param mod 
#' optional model object containing \code{delta}, \code{Gamma}, \code{allprobs}, and
#'   optionally \code{trackID}. When using \code{RTMB::MakeADFun} or \code{\link{qreml}} with
#'   \code{\link{forward_g}} in the likelihood, these are reported automatically after model fitting
#'   and the object returned by \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} can be passed directly.
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
  
  deprecated(
    "stateprobs_g() is deprecated, please use stateprobs() for all model types instead."
  )
  
  stateprobs(delta, Gamma, allprobs, trackID = trackID, mod = mod, forecast = forecast)
}


#' Calculate conditional local state probabilities for periodically inhomogeneous HMMs
#' 
#' Computes
#' \deqn{\Pr(S_t = j \mid X_1, ..., X_T)}
#' for periodically inhomogeneous HMMs
#' 
#' \strong{Deprecated:} please use \code{\link{stateprobs}} with argument \code{tod} instead.
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
#' @param trackID optional vector of k track IDs, if multiple tracks need to be decoded separately. The observations of each track need to be consecutive.
#' @param mod optional model object containing initial distribution \code{delta}, transition probability matrix \code{Gamma}, matrix of state-dependent probabilities \code{allprobs}, and potentially a \code{trackID} variable
#' 
#' If you are using automatic differentiation either with \code{RTMB::MakeADFun} or \code{\link{qreml}} and include \code{\link{forward_p}} in your likelihood function, the objects needed for state decoding are automatically reported after model fitting.
#' Hence, you can pass the model object obtained from \code{obj$report()}, \code{\link{lama}} or \code{\link{qreml}} directly to this function.
#' @param forecast logical, indicating if forecast probabilities \eqn{\Pr(S_t = j \mid X_1, ..., X_{t-1})} should be calculated instead.
#'
#' @return matrix of conditional state probabilities of dimension c(n,N)
#' @export
#'
#' @examples
#' L = 24
#' beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
#' Gamma = tpm_p(1:L, L, beta, degree = 1)
#' delta = stationary_p(Gamma, 2) # for tod[1] = 1
#' allprobs = matrix(runif(200), nrow = 100, ncol = 2)
#' tod = rep(1:24, 5)[1:100]
#' 
#' probs = stateprobs_p(delta, Gamma, allprobs, tod)

stateprobs_p = function(delta, Gamma, allprobs, tod, trackID = NULL,
                        mod = NULL, forecast = FALSE) {
  
  deprecated(
    "stateprobs_p() is deprecated, please use stateprobs(..., tod = tod) for periodic models instead."
  )
  
  stateprobs(delta, Gamma, allprobs, trackID = trackID, mod = mod, forecast = forecast, tod = tod)
}

