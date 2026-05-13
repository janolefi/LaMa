# Changelog

## LaMa 2.1.1 (2026-05-13)

CRAN release: 2026-05-13

### Major changes

- vignettes 1-3 have been updated to work with automatic differentiation

- tpm(), forward(), stateprobs(), viterbi() can now be used with
  inhomogeneous models: call their inhomogeneous counterparts internally

- function MCreport() has been added for sampling

  - from the distribution of the MLE (fixed effects only models)
  - from the joint post. distribution of parameters and random effects
    (random effects models)

### Minor changes

- some functions have been renamed, e.g. stationary_cont() -\>
  stationary_ct() (old versions still exist)

- some internal cleanup
