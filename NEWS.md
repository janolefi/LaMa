------------------------------------------------------------------------
LaMa 2.1.1 (2026-05-13)
------------------------------------------------------------------------

o vignettes 1-3 have been updated to work with automatic differentiation

o tpm(), forward(), stateprobs(), viterbi() can now be used with inhomogeneous models:
  call their inhomogeneous counterparts internally
  
o some functions have been renamed, e.g. stationary_cont() -> stationary_ct() 
  (old versions still exist)
  
o function MCreport() has been added for sampling 
  - from the distribution of the MLE (fixed effects only models)
  - from the joint post. distribution of parameters and random effects (random effects moels)
  
o some internal cleanup  