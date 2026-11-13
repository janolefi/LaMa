## Submission

Minor release 3.1.0. I apologise for submitting again so soon after 3.0.0. I need this 
version on CRAN to submit a paper whose code relies on it.

The changes concern mostly `qreml()`, the penalised likelihood estimation of models with smooth 
effects:

* new `joint_uncertainty()` for the joint uncertainty of the parameters and the penalty 
  strengths, which `MCreport()` can now account for (this uses the new `par` argument of 
  `RTMBdist::mcreport()`, hence `RTMBdist (>= 1.2.1)`)
* a more accurate outer Hessian w.r.t. the log penalty strengths
* more robust convergence (lower bound on the log penalty strengths, fewer unnecessary refits)
* new `coef()`, `vcov()` and `residuals()` methods for fitted models

Further, `tpm_ct()` is much faster for continuous-time models.

## Test environments

* local macOS (aarch64), R release
* GitHub Actions: macOS (release), Windows (release), Ubuntu (devel, release, oldrel-1)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

RTMBdist (Suggests) was checked against this version: no problems.
