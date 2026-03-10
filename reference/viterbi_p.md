# Viterbi algorithm for state decoding in periodically inhomogeneous HMMs

The Viterbi algorithm allows one to decode the most probable state
sequence of an HMM.

## Usage

``` r
viterbi_p(delta, Gamma, allprobs, tod, trackID = NULL, mod = NULL)
```

## Arguments

- delta:

  initial distribution of length N, or matrix of dimension c(k,N) for k
  independent tracks, if `trackID` is provided

  This could e.g. be the periodically stationary distribution (for each
  track).

- Gamma:

  array of transition probability matrices for each time point in the
  cycle of dimension c(N,N,L), where L is the length of the cycle

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N)

- tod:

  (Integer valued) variable for cycle indexing in 1, ..., L, mapping the
  data index to a generalised time of day (length n)

  For half-hourly data L = 48. It could, however, also be day of year
  for daily data and L = 365.

- trackID:

  optional vector of k track IDs, if multiple tracks need to be decoded
  separately

- mod:

  optional model object containing initial distribution `delta`,
  transition probability matrix `Gamma`, matrix of state-dependent
  probabilities `allprobs`, and potentially a `trackID` variable

  If you are using automatic differentiation either with
  `RTMB::MakeADFun` or
  [`qreml`](https://janolefi.github.io/reference/qreml.md) and include
  [`forward_p`](https://janolefi.github.io/reference/forward_p.md) in
  your likelihood function, the objects needed for state decoding are
  automatically reported after model fitting. Hence, you can pass the
  model object obtained from running `RTMB::report()` or from
  [`qreml`](https://janolefi.github.io/reference/qreml.md) directly to
  this function.

## Value

vector of decoded states of length n

## See also

Other decoding functions:
[`stateprobs()`](https://janolefi.github.io/reference/stateprobs.md),
[`stateprobs_g()`](https://janolefi.github.io/reference/stateprobs_g.md),
[`stateprobs_p()`](https://janolefi.github.io/reference/stateprobs_p.md),
[`viterbi()`](https://janolefi.github.io/reference/viterbi.md),
[`viterbi_g()`](https://janolefi.github.io/reference/viterbi_g.md)

## Examples

``` r
delta = c(0.5, 0.5)
beta = matrix(c(-2, 1, -1,
                -2, -1, 1), nrow = 2, byrow = TRUE)
Gamma = tpm_p(1:24, 24, beta)

tod = rep(1:24, 5)
n = length(tod)

allprobs = matrix(runif(2*n), nrow = n, ncol = 2)
states = viterbi_p(delta, Gamma, allprobs, tod)
```
