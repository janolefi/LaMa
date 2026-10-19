# LaMa (development version)

## Bug fixes

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