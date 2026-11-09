# Changelog

## LaMa 3.1.0

- [`tpm_ct()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md) is
  much faster: with AD and 2 states it uses an exact closed form instead
  of
  [`Matrix::expm()`](https://rdrr.io/pkg/Matrix/man/expm-methods.html)
  for each interval (about 100x faster gradients), and for a constant
  `Q` it computes the matrix exponential only once per unique time
  difference
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) no
  longer evaluates the Hessian a second time for the returned model,
  which saves about `length(par)` gradient evaluations per fit, and the
  final message now reads `"Final sp: ..."` as no fit follows it unless
  the last iterate was not the best one
- [`MCreport()`](https://janolefi.github.io/LaMa/reference/MCreport.md)
  for [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
  models samples around the final estimate, and with `joint = TRUE` (or
  `Q` from `joint_uncertainty()`) accounts for the uncertainty in the
  penalty strengths; needs RTMBdist \>= 1.2.1
- `residuals(mod)` gives the pseudo-residuals, calling
  `pseudo_res(mod = mod)`
- `coef(mod)` and `vcov(mod)` give the estimated fixed parameters and
  their covariance from the inverse Hessian (of the marginal likelihood
  with random effects, conditional on the smoothing parameters for
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
  models), so `confint(mod)` works too
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) no
  longer tries doubled steps after improving ones: each trial costs a
  full inner fit, and across the case studies tried it almost never
  reduced the number of inner fits
- `qreml(lsp_min = -15)` bounds `log(lambda)` from below, so a penalty
  strength on a smooth the data cannot inform stops falling and the
  iteration can converge; a penalty strength ending there gives a
  warning, and [`summary()`](https://rdrr.io/r/base/summary.html) notes
  it: a larger basis dimension is suggested if no other penalty acts on
  the same coefficients, otherwise the penalty is reported as redundant
  (e.g. the wiggliness penalty of a tensor product with a random effect
  margin), which a larger basis does not help
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md): a
  step that changes the restricted likelihood by less than
  `1e-6 * (1 + |criterion|)` counts as unchanged, instead of setting off
  up to `max_halve` refits at essentially the same penalty strengths at
  the optimum
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) also
  converges once the inner optimiser no longer moves from its warm
  start, so a penalty strength creeping towards zero along a flat
  restricted likelihood no longer runs on until `lsp_min` or `maxiter`
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
  returns `n_fits`, the number of inner fits including rejected steps,
  shown by [`summary()`](https://rdrr.io/r/base/summary.html) next to
  the number of outer iterations
- [`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md) and
  [`generator()`](https://janolefi.github.io/LaMa/reference/generator.md)
  with a design matrix give a one-time message while taping if
  `TapeConfig(matmul = 'plain')` would record more than 1e8
  multiply-adds for the linear predictor (counting only nonzero entries
  of `Z`, which RTMB skips otherwise), suggesting
  `TapeConfig(matmul = 'atomic')` as the tape may not fit into memory
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
  returns `edf_hist`, the effective degrees of freedom of each smooth at
  the start of each outer iteration (column 1 is the initial fit)
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md):
  `outer_hessian()` now differentiates the approximate gradient of the
  restricted log-likelihood w.r.t. `log(lambda)` by central differences
  (2 inner fits per free penalty strength, step `h = 0.1`), instead of
  one-sided differences with an absolute step of 1e-4 in `lambda`, which
  was inaccurate and could give an indefinite Hessian. The fits are warm
  started at the final estimate with a tighter inner tolerance. The
  result is now on the log scale
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md):
  `outer_hessian()` no longer leaves the objective at the last perturbed
  penalty strength, so samples or Hessians computed from `mod$obj`
  afterwards use the fitted penalty strengths again
- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
  returns `joint_uncertainty()`, the joint precision matrix of the
  parameters and the log penalty strengths (the first-order correction
  of Wood, Pya and Saefken, 2016, for the uncertainty in the penalty
  strengths, in the form of TMB’s `sdreport(getJointPrecision = TRUE)`),
  which stays sparse with `spHess = TRUE`; with `precision = FALSE` the
  corrected covariance of the parameters
- [`sdreport_outer()`](https://janolefi.github.io/LaMa/reference/sdreport_outer.md)
  works with the outer Hessian w.r.t. `log(lambda)` and obtains the
  standard deviations of `lambda` and `1/lambda` from it by the delta
  method

## LaMa 3.0.0

CRAN release: 2026-10-28

This is a major release: models are now built around
[`statedens()`](https://janolefi.github.io/LaMa/reference/statedens.md)
and [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md), which
give [`print()`](https://rdrr.io/r/base/print.html),
[`summary()`](https://rdrr.io/r/base/summary.html),
[`plot()`](https://rdrr.io/r/graphics/plot.default.html),
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`pseudo_res()`](https://janolefi.github.io/LaMa/reference/pseudo_res.md)
and [`scores()`](https://janolefi.github.io/LaMa/reference/scores.md)
from the model object alone. Several functions are renamed (old names
still work but are deprecated), a few are removed, and some outputs
change (e.g. `llk` instead of `ll`, `sp` instead of `lambda`), so
existing code may need small adjustments.

### New features

- [`statedens()`](https://janolefi.github.io/LaMa/reference/statedens.md)
  builds the matrix of state-dependent densities of one data stream for
  `allprobs`, with parameters given per state, for all states or per
  observation and state; missing observations are handled for each
  stream separately, and `log = TRUE` gives log-densities for
  `forward(..., logspace = TRUE)`; `pseudo_res(mod)` works with it

- [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  without AD gives an error for negative `allprobs`, suggesting
  `logspace = TRUE`

- `forward(delta = "stationary")` computes the initial distribution from
  `Gamma`: the stationary distribution (per track if track-specific, of
  the first slice of each track if time-varying) or, with `tod`, the
  periodically stationary distribution at the first time point of each
  track

- `plot(mod, type = "statedens")` shows the state-dependent densities of
  each data stream built with
  [`statedens()`](https://janolefi.github.io/LaMa/reference/statedens.md)
  on a histogram, weighted by the local state probabilities (or decoded
  states, the stationary distribution or given weights), which also
  works with covariates, several tracks and inhomogeneous models; the
  weighted densities are returned invisibly

- `plot(mod, stream = ...)` selects data streams, for residuals and
  state-dependent densities

- `plot(mod, type = "statedens", xlim = ...)` evaluates the density
  curves on the plotted range, not only the range of the data

- `simulate(mod)` simulates states and observations of models built with
  [`statedens()`](https://janolefi.github.io/LaMa/reference/statedens.md)
  on the observed design: same tracks and lengths, with `delta`, `Gamma`
  and the parameters as reported (conditional on the covariates, for
  continuous-time models on the observed time points); missing
  observations stay missing unless `fill_na = TRUE`, and `rdens` gives
  random number generators for densities without a matching r-function;
  a data stream that is `NA` throughout is a placeholder, so a
  likelihood function can be used to simulate data and then be fitted to
  them

- `simulate(mod, random = TRUE)` and `obj$simulate()` of RTMB work with
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  and
  [`statedens()`](https://janolefi.github.io/LaMa/reference/statedens.md):
  [`statedens()`](https://janolefi.github.io/LaMa/reference/statedens.md)
  registers its data like `OBS(x)`, so random effects (with their
  density evaluated first, e.g. `u %~% dnorm(...)`) are drawn first,
  then the states and data given them;
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  reports the simulated states as `states`

- [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rdrr.io/r/base/summary.html) for models from
  [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md) and
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md):
  [`print()`](https://rdrr.io/r/base/print.html) gives a short overview,
  [`summary()`](https://rdrr.io/r/base/summary.html) the estimated state
  process (transition probability or generator matrix, stationary
  distribution and mean dwell times, averaged over tracks if
  track-specific) and the reported quantities, with those of length
  `nStates` in one table;
  [`summary.qremlModel()`](https://janolefi.github.io/LaMa/reference/summary.qremlModel.md)
  adds the smooth terms and returns an object now

- [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) models
  have class `c("qremlModel", "LaMaModel")`, so all `LaMaModel` methods
  apply to them

- [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  reports type `"track-specific"` for one `Gamma` per track (before
  `"homogeneous"`)

- `multistart(obj)` optimises an RTMB object from several starting
  values (`obj$par` and perturbations of it, or of the first fit scaled
  by the unit information with `scale = "unit_info"`), stops by a
  Bayesian stopping rule or once the best optimum is reached often
  enough, and returns the best fit with a table of all runs; `lama(obj)`
  then builds the model from the best fit

- [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md)
  returns `sp0`, initial penalty strengths from
  [`mgcv::initial.sp()`](https://rdrr.io/pkg/mgcv/man/initial.sp.html)
  (times 100, capped at 1e6), one per penalty matrix and named as `S`;
  [`rep_sp()`](https://janolefi.github.io/LaMa/reference/rep_sp.md)
  repeats them as `lambda` for smooths used more than once, e.g. one per
  state

- The penalty strengths are now called `sp`:
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) takes
  them from `dat$sp` if present, otherwise from `dat$lambda` as before,
  so `lambda` is free for other things, e.g. a Poisson mean; `spname`
  still overrides this.
  [`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md)
  and
  [`penalty2()`](https://janolefi.github.io/LaMa/reference/penalty2.md)
  call their argument `sp`, `lambda` still works

- [`print()`](https://rdrr.io/r/base/print.html) for objects from
  [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md)
  shows the formulas and the initial penalty strengths

- [`smoothdens()`](https://janolefi.github.io/LaMa/reference/smoothdens.md)
  builds the objects for smooth (P-spline) densities from a data vector,
  or one object per stream from a data frame or list; `density(beta)` in
  the returned object gives the `allprobs` matrix and
  [`predict()`](https://rdrr.io/r/stats/predict.html) the estimated
  densities, e.g. on the grid `xseq`

### Changes

- [`make_matrices()`](https://janolefi.github.io/LaMa/reference/modmat.md)
  is renamed to
  [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md);
  [`make_matrices()`](https://janolefi.github.io/LaMa/reference/modmat.md)
  still works

- [`report()`](https://janolefi.github.io/LaMa/reference/lama.md) is
  renamed to
  [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md), as it
  does much more than calling `obj$report()`;
  [`report()`](https://janolefi.github.io/LaMa/reference/lama.md) still
  works but is deprecated

- old functions still work but are deprecated, with a message pointing
  to the new implementation:
  [`forward_s()`](https://janolefi.github.io/LaMa/reference/forward_s.md)
  and
  [`forward_sp()`](https://janolefi.github.io/LaMa/reference/forward_sp.md)
  (use
  [`forward_hsmm()`](https://janolefi.github.io/LaMa/reference/forward_hsmm.md)
  and
  [`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)),
  [`tpm_hsmm2()`](https://janolefi.github.io/LaMa/reference/tpm_hsmm2.md)
  and
  [`tpm_phsmm2()`](https://janolefi.github.io/LaMa/reference/tpm_phsmm2.md)
  (use
  [`tpm_hsmm()`](https://janolefi.github.io/LaMa/reference/tpm_hsmm.md)
  and
  [`tpm_phsmm()`](https://janolefi.github.io/LaMa/reference/tpm_phsmm.md)),
  [`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md) (use
  [`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md)),
  [`qreml_old()`](https://janolefi.github.io/LaMa/reference/qreml_old.md)
  (use [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)),
  [`stationary_cont()`](https://janolefi.github.io/LaMa/reference/stationary_ct.md)
  and
  [`tpm_cont()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md)
  (use
  [`stationary_ct()`](https://janolefi.github.io/LaMa/reference/stationary_ct.md)
  and [`tpm_ct()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md))

- deprecation messages are shown once per session, not at every call

- [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md): the
  log-likelihood is called `llk` (before `ll`), as for
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)

- [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md) and
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md):
  `nobs` does not count missing observations (rows of `allprobs` that
  are all 1), which changes BIC for data with missing observations

- [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md) stores
  the maximum gradient component `mgc` to check convergence

- `pseudo_res(mod)` computes pseudo-residuals from the model object
  alone: the likelihood function is evaluated again with its densities
  replaced by their CDFs, after checking that this is valid; several
  data streams give one set of residuals each; also for
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) models
  and densities in helper functions

- `scores(mod)` computes the score contribution of each observation for
  models using
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
  e.g. for robust (sandwich) standard errors; it takes the model from
  [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md) or the
  RTMB object, respects `map`, and the likelihood function does not need
  to change

- [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md) also
  stores `obj` (the same object, not a copy), so
  [`MCreport()`](https://janolefi.github.io/LaMa/reference/MCreport.md)
  can be called on the model object from
  [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md) or
  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)

- [`pred_matrix()`](https://janolefi.github.io/LaMa/reference/pred_matrix.md)
  is deprecated, use [`predict()`](https://rdrr.io/r/stats/predict.html)
  on the object from
  [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md)
  instead

- `make_matrices_old()` is removed

- [`smooth_dens_construct()`](https://janolefi.github.io/LaMa/reference/smooth_dens_construct.md)
  is deprecated, use
  [`smoothdens()`](https://janolefi.github.io/LaMa/reference/smoothdens.md)
  instead

- [`make_matrices_dens()`](https://janolefi.github.io/LaMa/reference/make_matrices_dens.md)
  and
  [`smooth_dens_construct()`](https://janolefi.github.io/LaMa/reference/smooth_dens_construct.md)
  no longer print messages; `basis$knots` includes the boundary knots

### Bug fixes

- Examples and vignettes: the periodically stationary initial
  distribution matches
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)’s
  `tod`, where `Gamma[,,tod[t]]` is the transition into `t`,
  i.e. `stationary_p(Gamma, tod[1] %% L + 1)`

- [`lama()`](https://janolefi.github.io/LaMa/reference/lama.md): the
  log-likelihood was evaluated at the last parameters the optimiser
  tried, not the best

- [`summary.qremlModel()`](https://janolefi.github.io/LaMa/reference/summary.qremlModel.md)
  failed for reported lists or character vectors, and if `Gamma` was not
  reported but `delta` was

- [`pseudo_res()`](https://janolefi.github.io/LaMa/reference/pseudo_res.md):
  missing observations give `NA` in all variants (before `0` with
  `normal = FALSE`, `-Inf` for non-randomised discrete residuals)

- [`pseudo_res()`](https://janolefi.github.io/LaMa/reference/pseudo_res.md):
  parameters given per observation (matrix or vector of length `nObs`)
  work with missing observations

- [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md):
  with track-specific `delta` or `Gamma`, the AD version matched them to
  the wrong tracks if the track IDs were not sorted
  (e.g. `"ID1", ..., "ID10"`)

- [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md):
  tracks of length one failed with AD

- `logspace = TRUE` is now exact (`delta` is no longer modified) and
  also stable without AD

- one-state models are correct with `logspace = TRUE` and report their
  quantities

- [`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md)
  works if `tod` does not contain all values `1, ..., L`

- [`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md)
  works without AD for two observations

- `stateprobs_g(forecast = TRUE)` with `trackID` used the transition
  matrix of the previous time point, which also affected
  [`pseudo_res()`](https://janolefi.github.io/LaMa/reference/pseudo_res.md)
  for inhomogeneous models with tracks

- [`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
  passes `forecast` on for time-varying `Gamma`

- [`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md)
  and
  [`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
  work for tracks of length one

- `stateprobs(forecast = TRUE)`: the first row of each track is `delta`,
  not the filtered probabilities given `x_1`, so
  [`pseudo_res()`](https://janolefi.github.io/LaMa/reference/pseudo_res.md)
  is correct for the first observation

- [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md):
  smooths with several penalties but not one per margin
  (e.g. `bs = "fs"`, `"ad"` or `t2(full = TRUE)`) kept only the first,
  so later smooths got the wrong penalty matrices

- [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md): a
  nested formula list failed if one data stream had no smooth terms

- [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md):
  tensor product margins over several covariates are named e.g. `"x,y"`
  (two such margins failed); unpenalised smooths (`fx = TRUE`) give an
  informative error

- [`predict()`](https://rdrr.io/r/stats/predict.html) on objects from
  [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md)
  ignored `exclude`

- [`smooth_dens_construct()`](https://janolefi.github.io/LaMa/reference/smooth_dens_construct.md):
  coefficients of one state are a matrix (before a vector, so
  `cbind(coef, 0)` was wrong); circular data with custom knots failed;
  angles outside \[-pi, pi) are wrapped; `par` accepts `kappa` next to
  `concentration`

- [`make_matrices_dens()`](https://janolefi.github.io/LaMa/reference/make_matrices_dens.md):
  the basis is normalised on a midpoint grid (before a fixed 1e-3
  offset, which gave negative weights for data on a small scale), which
  changes the weights very slightly

### Minor changes

- [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
  [`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md)
  and
  [`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md)
  share one implementation (one C++ and one R version for AD);
  [`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md)
  no longer expands `Gamma` for the computation
- tracks need to be contiguous, otherwise
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  gives an informative error
- without AD,
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  gives an error for `NA` in `allprobs`
- the `TapeConfig` performance tip is shown only once per session
- [`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md),
  [`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
  and their `_g` and `_p` versions use the inputs of
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md)
  and one C++ implementation each, so `Gamma` is no longer expanded; as
  in
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
  tracks need to be contiguous
- [`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md)
  reports `Gamma` with its `L` slices and type `"periodic"`
- [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
  [`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md)
  and
  [`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
  have an argument `tod` for periodic models, so they handle all model
  types;
  [`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md),
  [`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md),
  [`viterbi_g()`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
  [`viterbi_p()`](https://janolefi.github.io/LaMa/reference/viterbi_p.md),
  [`stateprobs_g()`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md)
  and
  [`stateprobs_p()`](https://janolefi.github.io/LaMa/reference/stateprobs_p.md)
  are deprecated and give a message
- [`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md)
  and
  [`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md)
  also take the model object as first argument, e.g. `viterbi(mod)`
  instead of `viterbi(mod = mod)`
- [`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md) and
  [`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md)
  are deprecated in favour of `tpm(beta, Z)` and `generator(beta, Z)`,
  which also accept a single covariate vector `Z`; the `sparse` argument
  of [`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md) is
  gone
- [`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md) sets
  up each basis once and no longer keeps `gam$X` (returned as `Z`), so
  it is faster and needs less memory

## LaMa 2.1.1 (2026-05-13)

CRAN release: 2026-05-13

### Major changes

- vignettes 1-3 have been updated to work with automatic differentiation

- [`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
  [`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
  [`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md),
  [`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md)
  can now be used with inhomogeneous models: call their inhomogeneous
  counterparts internally

- function
  [`MCreport()`](https://janolefi.github.io/LaMa/reference/MCreport.md)
  has been added for sampling

  - from the distribution of the MLE (fixed effects only models)
  - from the joint post. distribution of parameters and random effects
    (random effects models)

### Minor changes

- some functions have been renamed,
  e.g. [`stationary_cont()`](https://janolefi.github.io/LaMa/reference/stationary_ct.md)
  -\>
  [`stationary_ct()`](https://janolefi.github.io/LaMa/reference/stationary_ct.md)
  (old versions still exist)

- some internal cleanup
