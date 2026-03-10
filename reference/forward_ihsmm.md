# [Forward algorithm](https://www.taylorfrancis.com/books/mono/10.1201/b20790/hidden-markov-models-time-series-walter-zucchini-iain-macdonald-roland-langrock) for hidden semi-Markov models with inhomogeneous state durations and/ or conditional transition probabilities

Calculates the (approximate) log-likelihood of a sequence of
observations under an inhomogeneous hidden semi-Markov model using a
modified **forward algorithm**.

## Usage

``` r
forward_ihsmm(
  dm,
  omega,
  allprobs,
  trackID = NULL,
  delta = NULL,
  startInd = NULL,
  eps = 1e-10,
  report = TRUE
)
```

## Arguments

- dm:

  list of length N containing matrices (or vectors) of dwell-time
  probability mass functions (PMFs) for each state.

  If the dwell-time PMFs are constant, the vectors are the PMF of the
  dwell-time distribution fixed in time. The vector lengths correspond
  to the approximating state aggregate sizes, hence there should be
  little probablity mass not covered by these.

  If the dwell-time PMFs are inhomogeneous, the matrices need to have n
  rows, where n is the number of observations. The number of columns
  again correponds to the size of the approximating state aggregates.

  In the latter case, the first `max(sapply(dm, ncol)) - 1` observations
  will not be used because the first approximating transition
  probability matrix needs to be computed based on the first
  `max(sapply(dm, ncol))` covariate values (represented by `dm`).

- omega:

  matrix of dimension c(N,N) or array of dimension c(N,N,n) of
  conditional transition probabilites, also called embedded transition
  probability matrix.

  It contains the transition probabilities given the current state is
  left. Hence, the diagonal elements need to be zero and the rows need
  to sum to one. Such a matrix can be constructed using
  [`tpm_emb`](https://janolefi.github.io/reference/tpm_emb.md) and an
  array using
  [`tpm_emb_g`](https://janolefi.github.io/reference/tpm_emb_g.md).

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N)

- trackID:

  trackID optional vector of length n containing IDs

  If provided, the total log-likelihood will be the sum of each track's
  likelihood contribution. Instead of a single vector `delta`
  corresponding to the initial distribution, a `delta` matrix of initial
  distributions, of dimension c(k,N), can be provided, such that each
  track starts with it's own initial distribution.

- delta:

  optional vector of initial state probabilities of length N

  By default, instead of this, the stationary distribution is computed
  corresponding to the first approximating transition probability matrix
  of each track is computed. Contrary to the homogeneous case, this is
  not theoretically motivated but just for convenience.

- startInd:

  optional integer index at which the forward algorithm starts.

  When approximating inhomogeneous HSMMs by inhomogeneous HMMs, the
  first transition probability matrix that can be constructed is at time
  `max(sapply(dm, ncol))` (as it depends on the previous covariate
  values). Hence, when not provided, `startInd` is chosen to be
  `max(sapply(dm, ncol))`. Fixing `startInd` at a value **larger** than
  max(aggregate sizes) is useful when models with different aggregate
  sizes are fitted to the same data and are supposed to be compared. In
  that case it is important that all models use the same number of
  observations.

- eps:

  small value to avoid numerical issues in the approximating transition
  matrix construction. Usually, this should not be changed.

- report:

  logical, indicating whether initial distribution, approximating
  transition probability matrix and `allprobs` matrix should be reported
  from the fitted model. Defaults to `TRUE`.

## Value

log-likelihood for given data and parameters

## Details

Hidden semi-Markov models (HSMMs) are a flexible extension of HMMs,
where the state duration distribution is explicitly modelled by a
distribution on the positive integers. This function can be used to fit
HSMMs where the state-duration distribution and/ or the conditional
transition probabilities vary with covariates. For direct numerical
maximum likelhood estimation, HSMMs can be represented as HMMs on an
enlarged state space (of size \\M\\) and with structured transition
probabilities.

This function is designed to be used with automatic differentiation
based on the `R` package `RTMB`. It will be very slow without it!

## References

Koslik, J. O. (2025). Hidden semi-Markov models with inhomogeneous state
dwell-time distributions. Computational Statistics & Data Analysis, 209,
108171.

## See also

Other forward algorithms:
[`forward()`](https://janolefi.github.io/reference/forward.md),
[`forward2()`](https://janolefi.github.io/reference/forward2.md),
[`forward_g()`](https://janolefi.github.io/reference/forward_g.md),
[`forward_g2()`](https://janolefi.github.io/reference/forward_g2.md),
[`forward_hsmm()`](https://janolefi.github.io/reference/forward_hsmm.md),
[`forward_p()`](https://janolefi.github.io/reference/forward_p.md),
[`forward_phsmm()`](https://janolefi.github.io/reference/forward_phsmm.md)

## Examples

``` r
# currently no examples
```
