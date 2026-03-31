# [Forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) for hidden semi-Markov models with periodically inhomogeneous state durations and/ or conditional transition probabilities

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs,
where the state duration distribution is explicitly modelled by a
distribution on the positive integers. This function can be used to fit
HSMMs where the state-duration distribution and/ or the conditional
transition probabilities vary with covariates. For direct numerical
maximum likelhood estimation, HSMMs can be represented as HMMs on an
enlarged state space (of size \\M\\) and with structured transition
probabilities.

This function can be used to fit HSMMs where the state-duration
distribution and/ or the conditional transition probabilities vary
periodically. In the special case of periodic variation (as compared to
arbitrary covariate influence), this version is to be preferred over
[`forward_ihsmm`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md)
because it computes the **correct periodically stationary distribution**
and no observations are lost for the approximation.

This function is designed to be used with automatic differentiation
based on the `R` package `RTMB`. It will be very slow without it!

## Usage

``` r
forward_phsmm(
  dm,
  omega,
  allprobs,
  tod,
  trackID = NULL,
  delta = NULL,
  eps = 1e-10,
  report = TRUE
)
```

## Arguments

- dm:

  list of length `nStates` containing matrices (or vectors) of
  dwell-time probability mass functions (PMFs) for each state.

  If the dwell-time PMFs are constant, the vectors are the PMF of the
  dwell-time distribution fixed in time. The vector lengths correspond
  to the approximating state aggregate sizes, hence there should be
  little probablity mass not covered by these.

  If the dwell-time PMFs are inhomogeneous, the matrices need to have L
  rows, where L is the cycle length. The number of columns again
  correpond to the size of the approximating state aggregates.

- omega:

  matrix of dimension `c(nStates, nStates)` or array of dimension
  `c(nStates, nStates, L)` of conditional transition probabilites, also
  called embedded transition probability matrix

  It contains the transition probabilities given the current state is
  left. Hence, the diagonal elements need to be zero and the rows need
  to sum to one. Such a matrix can be constructed using
  [`tpm_emb`](https://janolefi.github.io/LaMa/reference/tpm_emb.md) and
  an array using
  [`tpm_emb_g`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md).

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  `c(nObs, nStates)`

- tod:

  (Integer valued) variable for cycle indexing in 1, ..., L, mapping the
  data index to a generalised time of day (length `nObs`). For
  half-hourly data L = 48. It could, however, also be day of year for
  daily data and L = 365.

- trackID:

  optional vector of length `nObs` containing IDs

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. Instead of a single vector `delta`
  corresponding to the initial distribution, a `delta` matrix of initial
  distributions, of dimension `c(nTracks, nStates)`, can be provided,
  such that each track starts with it's own initial distribution.

- delta:

  Optional vector of initial state probabilities of length `nStates`. By
  default, instead of this, the stationary distribution is computed
  corresponding to the first approximating t.p.m. of each track is
  computed. Contrary to the homogeneous case, this is not theoretically
  motivated but just for convenience.

- eps:

  small value to avoid numerical issues in the approximating transition
  matrix construction. Usually, this should not be changed.

- report:

  logical, indicating whether initial distribution, approximating
  transition probability matrix and `allprobs` matrix should be reported
  from the fitted model. Defaults to `TRUE`.

## Value

HSMM log-likelihood for given data and parameters

## Details

Calculates the (approximate) log-likelihood of a sequence of
observations under a periodically inhomogeneous hidden semi-Markov model
using a modified **forward algorithm**.

## References

Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state
dwell-time distributions. Computational Statistics & Data Analysis, 209,
108171.

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/LaMa/reference/forward.md),
[`forward_g()`](https://janolefi.github.io/LaMa/reference/forward_g.md),
[`forward_hsmm()`](https://janolefi.github.io/LaMa/reference/forward_hsmm.md),
[`forward_ihsmm()`](https://janolefi.github.io/LaMa/reference/forward_ihsmm.md),
[`forward_p()`](https://janolefi.github.io/LaMa/reference/forward_p.md)

## Examples

``` r
# currently no examples
```
