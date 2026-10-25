# LaMa (development version)

## New features

* statedens() builds the matrix of state-dependent densities of one data stream for allprobs, 
  with parameters given per state, for all states or per observation and state; missing 
  observations are handled for each stream separately, and log = TRUE gives log-densities 
  for forward(..., logspace = TRUE); pseudo_res(mod) works with it
* forward() without AD gives an error for negative allprobs, suggesting logspace = TRUE
* plot(mod, type = "statedens") shows the state-dependent densities of each data stream built with 
  statedens() on a histogram, weighted by the local state probabilities (or decoded states, the 
  stationary distribution or given weights), which also works with covariates, several tracks and 
  inhomogeneous models; the weighted densities are returned invisibly
* plot(mod, stream = ...) selects data streams, for residuals and state-dependent densities
* simulate(mod) simulates states and observations of models built with statedens() on the observed design: 
  same tracks and lengths, with delta, Gamma and the parameters as reported (conditional on the covariates, 
  for continuous-time models on the observed time points); missing observations stay missing, and rdens 
  gives random number generators for densities without a matching r-function; a data stream that is NA 
  throughout is a placeholder, so a likelihood function can be used to simulate data and then be fitted to them
* simulate(mod, random = TRUE) and obj$simulate() of RTMB work with forward() and statedens(): for data marked 
  by x <- OBS(x), random effects (declared by u %~% dnorm(...)) are drawn first, then the states and data 
  given them; forward() reports the simulated states as states
* print() and summary() for models from report() and qreml(): print() gives a short overview, 
  summary() the estimated state process (transition probability or generator matrix, stationary 
  distribution and mean dwell times, averaged over tracks if track-specific) and the reported 
  quantities, with those of length nStates in one table; summary.qremlModel() adds the smooth 
  terms and returns an object now
* qreml() models have class c("qremlModel", "LaMaModel"), so all LaMaModel methods apply to them
* forward() reports type "track-specific" for one Gamma per track (before "homogeneous")

## Changes

* report(): the log-likelihood is called llk (before ll), as for qreml()
* report() and qreml(): nobs does not count missing observations (rows of allprobs that are all 1), 
  which changes BIC for data with missing observations
* report() stores the maximum gradient component mgc to check convergence

* pseudo_res(mod) computes pseudo-residuals from the model object alone: the likelihood function 
  is evaluated again with its densities replaced by their CDFs, after checking that this is valid; 
  several data streams give one set of residuals each; also for qreml() models and densities in 
  helper functions

* scores(mod) computes the score contribution of each observation for models using forward(), 
  e.g. for robust (sandwich) standard errors; it takes the model from report() or the RTMB 
  object, respects map, and the likelihood function does not need to change
* report() also stores obj (the same object, not a copy), so MCreport() can be called on 
  the model object from report() or qreml()

## Bug fixes

* report(): the log-likelihood was evaluated at the last parameters the optimiser tried, not the best
* summary.qremlModel() failed for reported lists or character vectors, and if Gamma was not reported
  but delta was

* pseudo_res(): missing observations give NA in all variants (before 0 with normal = FALSE, 
  -Inf for non-randomised discrete residuals)
* pseudo_res(): parameters given per observation (matrix or vector of length nObs) work with 
  missing observations

* forward(): with track-specific delta or Gamma, the AD version matched them to the 
  wrong tracks if the track IDs were not sorted (e.g. "ID1", ..., "ID10")
* forward(): tracks of length one failed with AD
* logspace = TRUE is now exact (delta is no longer modified) and also stable without AD
* one-state models are correct with logspace = TRUE and report their quantities
* forward_p() works if tod does not contain all values 1, ..., L
* forward_g() works without AD for two observations
* stateprobs_g(forecast = TRUE) with trackID used the transition matrix of the previous 
  time point, which also affected pseudo_res() for inhomogeneous models with tracks
* stateprobs() passes forecast on for time-varying Gamma
* viterbi() and stateprobs() work for tracks of length one
* stateprobs(forecast = TRUE): the first row of each track is delta, not the filtered 
  probabilities given x_1, so pseudo_res() is correct for the first observation

## Minor changes

* forward(), forward_g() and forward_p() share one implementation (one C++ and one 
  R version for AD); forward_p() no longer expands Gamma for the computation
* tracks need to be contiguous, otherwise forward() gives an informative error
* without AD, forward() gives an error for NA in allprobs
* the TapeConfig performance tip is shown only once per session
* viterbi(), stateprobs() and their _g and _p versions use the inputs of forward() and one 
  C++ implementation each, so Gamma is no longer expanded; as in forward(), tracks need 
  to be contiguous
* forward_p() reports Gamma with its L slices and type "periodic"
* forward(), viterbi() and stateprobs() have an argument tod for periodic models, so they 
  handle all model types; forward_g(), forward_p(), viterbi_g(), viterbi_p(), stateprobs_g() 
  and stateprobs_p() are deprecated and give a message
* viterbi() and stateprobs() also take the model object as first argument, 
  e.g. viterbi(mod) instead of viterbi(mod = mod)
* tpm_g() and generator_g() are deprecated in favour of tpm(beta, Z) and generator(beta, Z), 
  which also accept a single covariate vector Z; the sparse argument of tpm_g() is gone

# LaMa 2.1.1 (2026-05-13)

## Major changes

* vignettes 1-3 have been updated to work with automatic differentiation

* tpm(), forward(), stateprobs(), viterbi() can now be used with inhomogeneous models:
  call their inhomogeneous counterparts internally
  
* function MCreport() has been added for sampling 
  - from the distribution of the MLE (fixed effects only models)
  - from the joint post. distribution of parameters and random effects (random effects models)
  
## Minor changes
 
* some functions have been renamed, e.g. stationary_cont() -> stationary_ct() 
  (old versions still exist)
  
* some internal cleanup  