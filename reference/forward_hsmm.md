# [Forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) for hidden semi-Markov models

Calculates the (approximate) log-likelihood of a sequence of
observations under a hidden semi-Markov model using a modified **forward
algorithm**. Dwell-time distributions and conditional transition
probabilities can be fixed, track-specific, time-varying or periodic.

## Usage

``` r
forward_hsmm(
  dm,
  omega,
  allprobs,
  trackID = NULL,
  delta = NULL,
  eps = 1e-10,
  report = TRUE,
  tod = NULL,
  use_all_obs = FALSE,
  first_obs = NULL,
  logspace = FALSE
)
```

## Arguments

- dm:

  list of length `nStates` containing the dwell-time probability mass
  functions (PMFs) of the states, where each element is either

  - a vector, if the dwell-time distribution is fixed, or

  - a matrix with `nTracks` rows (track-specific), `nObs` rows
    (time-varying), or `L` rows (periodic, if `tod` is provided), where
    each row is a PMF.

  The vector lengths or numbers of columns correspond to the
  approximating state aggregate sizes, hence there should be little
  probability mass not covered by these.

- omega:

  matrix of dimension `c(nStates, nStates)` of conditional transition
  probabilities, also called embedded transition probability matrix, or
  array with `nTracks`, `nObs` or `L` of such matrices as slices.

  Contains the transition probabilities given that the current state is
  left. Hence, the diagonal elements need to be zero and the rows need
  to sum to one. Can be constructed using
  [`tpm_switch`](https://janolefi.github.io/LaMa/reference/tpm_switch.md).

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  `c(nObs, nStates)` which will automatically be converted to the
  appropriate dimension.

- trackID:

  optional vector of length `nObs` containing IDs

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. In this case, `dm` and `omega` can be
  track-specific and `delta` can be a matrix of initial distributions,
  of dimension `c(nTracks, nStates)`, such that each track starts with
  its own initial distribution.

- delta:

  optional vector of initial state probabilities of length `nStates`, in
  which case a new stay begins at the first observation used in each
  track.

  By default, the stationary distribution is computed (which is
  typically recommended).

- eps:

  small value to avoid numerical issues in the approximating transition
  matrix construction. Usually, this should not be changed.

- report:

  logical, indicating whether initial distribution, `dm`, `omega` and
  the `allprobs` matrix should be reported from the fitted model.
  Defaults to `TRUE`.

- tod:

  optional integer vector of length `nObs` containing the time of day
  (or any other cycle), with values in `1, ..., L`, for periodic models.

  With periodic `dm`, it needs to advance by one in each row within a
  track.

- use_all_obs:

  logical, only relevant for time-varying `dm`.

  By default, the first `max(sapply(dm, ncol)) - 1` observations of each
  track are not used, as they may belong to stays that began before the
  first observation, with unknown dwell-time distributions. If `TRUE`,
  all observations are used, which needs `delta`.

- first_obs:

  optional integer, the first observation that is used in each track.

  By default `1`, or `max(sapply(dm, ncol))` for time-varying `dm`, in
  which case it cannot be smaller. Useful to fit models with different
  aggregate sizes or types of `dm` to the same observations, such that
  their likelihoods can be compared.

- logspace:

  logical, indicating whether `allprobs` contains log-probabilities.
  Defaults to `FALSE`.

## Value

HSMM log-likelihood for given data and parameters

## Details

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs,
where the state duration distribution is explicitly modelled by a
distribution on the positive integers. For direct numerical maximum
likelihood estimation, HSMMs can be represented as HMMs on an enlarged
state space (of size \\M\\) and with structured transition
probabilities.

A stay that begins at time point `t` is governed by `omega[,,t]` (which
state is entered) and row `t` of `dm` (how long the stay lasts).

This function is designed to be used with automatic differentiation
based on the `R` package `RTMB`. It will be very slow without it!

## References

Langrock, R., & Zucchini, W. (2011). Hidden Markov models with arbitrary
state dwell-time distributions. Computational Statistics & Data
Analysis, 55(1), 715-724.

Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state
dwell-time distributions. Computational Statistics & Data Analysis, 209,
108171.

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
[`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md),
[`forward_phsmm()`](https://janolefi.github.io/LaMa/reference/forward_phsmm.md)

## Examples

``` r
## negative log-likelihood of a 2-state HSMM with shifted Poisson dwell times
nll = function(par) {
  getAll(par)
  sizes = c(20, 20) # state aggregate sizes
  dm = lapply(1:2, function(i) dpois(1:sizes[i] - 1, exp(loglambda[i])))
  omega = matrix(c(0, 1, 1, 0), 2) # 2 states: no parameters
  allprobs = matrix(1, length(step), 2)
  ind = which(!is.na(step))
  for(j in 1:2) allprobs[ind, j] = dgamma2(step[ind], exp(logmu[j]), exp(logsigma[j]))
  -forward_hsmm(dm, omega, allprobs)
}
step = trex$step[1:200]
par = list(loglambda = log(c(5, 5)), logmu = log(c(0.3, 2.5)), logsigma = log(c(0.2, 1.5)))
obj = MakeADFun(nll, par, silent = TRUE)
opt = nlminb(obj$par, obj$fn, obj$gr)
```
