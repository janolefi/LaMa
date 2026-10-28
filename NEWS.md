# LaMa 3.0.0

This is a major release: models are now built around statedens() and lama(), which give 
print(), summary(), plot(), simulate(), pseudo_res() and scores() from the model object alone. 
Several functions are renamed (old names still work but are deprecated), a few are removed, 
and some outputs change (e.g. llk instead of ll, sp instead of lambda), so existing code may 
need small adjustments.

## New features

* statedens() builds the matrix of state-dependent densities of one data stream for allprobs, 
  with parameters given per state, for all states or per observation and state; missing 
  observations are handled for each stream separately, and log = TRUE gives log-densities 
  for forward(..., logspace = TRUE); pseudo_res(mod) works with it
* forward() without AD gives an error for negative allprobs, suggesting logspace = TRUE
* forward(delta = "stationary") computes the initial distribution from Gamma: the stationary distribution 
  (per track if track-specific, of the first slice of each track if time-varying) or, with tod, the 
  periodically stationary distribution at the first time point of each track
* plot(mod, type = "statedens") shows the state-dependent densities of each data stream built with 
  statedens() on a histogram, weighted by the local state probabilities (or decoded states, the 
  stationary distribution or given weights), which also works with covariates, several tracks and 
  inhomogeneous models; the weighted densities are returned invisibly
* plot(mod, stream = ...) selects data streams, for residuals and state-dependent densities
* plot(mod, type = "statedens", xlim = ...) evaluates the density curves on the plotted range, not only the range of the data
* simulate(mod) simulates states and observations of models built with statedens() on the observed design: 
  same tracks and lengths, with delta, Gamma and the parameters as reported (conditional on the covariates, 
  for continuous-time models on the observed time points); missing observations stay missing unless fill_na = TRUE, and rdens 
  gives random number generators for densities without a matching r-function; a data stream that is NA 
  throughout is a placeholder, so a likelihood function can be used to simulate data and then be fitted to them
* simulate(mod, random = TRUE) and obj$simulate() of RTMB work with forward() and statedens(): statedens() 
  registers its data like OBS(x), so random effects (with their density evaluated first, e.g. u %~% dnorm(...)) 
  are drawn first, then the states and data given them; forward() reports the simulated states as states
* print() and summary() for models from lama() and qreml(): print() gives a short overview, 
  summary() the estimated state process (transition probability or generator matrix, stationary 
  distribution and mean dwell times, averaged over tracks if track-specific) and the reported 
  quantities, with those of length nStates in one table; summary.qremlModel() adds the smooth 
  terms and returns an object now
* qreml() models have class c("qremlModel", "LaMaModel"), so all LaMaModel methods apply to them
* forward() reports type "track-specific" for one Gamma per track (before "homogeneous")
* multistart(obj) optimises an RTMB object from several starting values (obj$par and perturbations 
  of it, or of the first fit scaled by the unit information with scale = "unit_info"), stops by a 
  Bayesian stopping rule or once the best optimum is reached often enough, and returns the best fit 
  with a table of all runs; lama(obj) then builds the model from the best fit

* modmat() returns sp0, initial penalty strengths from mgcv::initial.sp() (times 100, 
  capped at 1e6), one per penalty matrix and named as S; rep_sp() repeats them as lambda for 
  smooths used more than once, e.g. one per state
* The penalty strengths are now called sp: qreml() takes them from dat$sp if present, otherwise 
  from dat$lambda as before, so lambda is free for other things, e.g. a Poisson mean; spname 
  still overrides this. penalty() and penalty2() call their argument sp, lambda still works
* print() for objects from modmat() shows the formulas and the initial penalty strengths
* smoothdens() builds the objects for smooth (P-spline) densities from a data vector, or one object per 
  stream from a data frame or list; density(beta) in the returned object gives the allprobs matrix 
  and predict() the estimated densities, e.g. on the grid xseq

## Changes

* make_matrices() is renamed to modmat(); make_matrices() still works
* report() is renamed to lama(), as it does much more than calling obj$report(); report() still works but is deprecated
* old functions still work but are deprecated, with a message pointing to the new implementation: 
  forward_s() and forward_sp() (use forward_hsmm() and forward_phsmm()), tpm_hsmm2() and tpm_phsmm2() 
  (use tpm_hsmm() and tpm_phsmm()), tpm_g2() (use tpm()), qreml_old() (use qreml()), 
  stationary_cont() and tpm_cont() (use stationary_ct() and tpm_ct())
* deprecation messages are shown once per session, not at every call
* lama(): the log-likelihood is called llk (before ll), as for qreml()
* lama() and qreml(): nobs does not count missing observations (rows of allprobs that are all 1), 
  which changes BIC for data with missing observations
* lama() stores the maximum gradient component mgc to check convergence

* pseudo_res(mod) computes pseudo-residuals from the model object alone: the likelihood function 
  is evaluated again with its densities replaced by their CDFs, after checking that this is valid; 
  several data streams give one set of residuals each; also for qreml() models and densities in 
  helper functions

* scores(mod) computes the score contribution of each observation for models using forward(), 
  e.g. for robust (sandwich) standard errors; it takes the model from lama() or the RTMB 
  object, respects map, and the likelihood function does not need to change
* lama() also stores obj (the same object, not a copy), so MCreport() can be called on 
  the model object from lama() or qreml()

* pred_matrix() is deprecated, use predict() on the object from modmat() instead
* make_matrices_old() is removed
* smooth_dens_construct() is deprecated, use smoothdens() instead
* make_matrices_dens() and smooth_dens_construct() no longer print messages; basis$knots 
  includes the boundary knots

## Bug fixes

* Examples and vignettes: the periodically stationary initial distribution matches forward()'s tod, 
  where Gamma[,,tod[t]] is the transition into t, i.e. stationary_p(Gamma, tod[1] %% L + 1)
* lama(): the log-likelihood was evaluated at the last parameters the optimiser tried, not the best
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

* modmat(): smooths with several penalties but not one per margin (e.g. bs = "fs", "ad" 
  or t2(full = TRUE)) kept only the first, so later smooths got the wrong penalty matrices
* modmat(): a nested formula list failed if one data stream had no smooth terms
* modmat(): tensor product margins over several covariates are named e.g. "x,y" (two 
  such margins failed); unpenalised smooths (fx = TRUE) give an informative error
* predict() on objects from modmat() ignored exclude

* smooth_dens_construct(): coefficients of one state are a matrix (before a vector, so cbind(coef, 0) 
  was wrong); circular data with custom knots failed; angles outside [-pi, pi) are wrapped; 
  par accepts kappa next to concentration
* make_matrices_dens(): the basis is normalised on a midpoint grid (before a fixed 1e-3 offset, 
  which gave negative weights for data on a small scale), which changes the weights very slightly

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
* modmat() sets up each basis once and no longer keeps gam$X (returned as Z), 
  so it is faster and needs less memory

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
