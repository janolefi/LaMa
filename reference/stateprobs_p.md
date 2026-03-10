# Calculate conditional local state probabilities for periodically inhomogeneous HMMs

Computes \$\$\Pr(S_t = j \mid X_1, ..., X_T)\$\$ for periodically
inhomogeneous HMMs

## Usage

``` r
stateprobs_p(
  delta,
  Gamma,
  allprobs,
  tod,
  trackID = NULL,
  mod = NULL,
  forecast = FALSE
)
```

## Arguments

- delta:

  initial or stationary distribution of length N, or matrix of dimension
  c(k,N) for k independent tracks, if `trackID` is provided

  This could e.g. be the periodically stationary distribution (for each
  track) as computed by
  [`stationary_p`](https://janolefi.github.io/LaMa/reference/stationary_p.md).

- Gamma:

  array of transition probability matrices for each time point in the
  cycle of dimension c(N,N,L), where L is the length of the cycle.

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N)

- tod:

  (Integer valued) variable for cycle indexing in 1, ..., L, mapping the
  data index to a generalised time of day (length n). For half-hourly
  data L = 48. It could, however, also be day of year for daily data and
  L = 365.

- trackID:

  optional vector of k track IDs, if multiple tracks need to be decoded
  separately

- mod:

  optional model object containing initial distribution `delta`,
  transition probability matrix `Gamma`, matrix of state-dependent
  probabilities `allprobs`, and potentially a `trackID` variable

  If you are using automatic differentiation either with
  `RTMB::MakeADFun` or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) and
  include
  [`forward_p`](https://janolefi.github.io/LaMa/reference/forward_p.md)
  in your likelihood function, the objects needed for state decoding are
  automatically reported after model fitting. Hence, you can pass the
  model object obtained from running `RTMB::report()` or from
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) directly
  to this function.

- forecast:

  logical, indicating if forecast probabilities \\\Pr(S_t = j \mid X_1,
  ..., X_t)\\ should be calculated instead.

## Value

matrix of conditional state probabilities of dimension c(n,N)

## See also

Other decoding functions:
[`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md),
[`stateprobs_g()`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md),
[`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md),
[`viterbi_g()`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
[`viterbi_p()`](https://janolefi.github.io/LaMa/reference/viterbi_p.md)

## Examples

``` r
L = 24
beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
Gamma = tpm_p(1:L, L, beta, degree = 1)
delta = stationary_p(Gamma, 1)
allprobs = matrix(runif(200), nrow = 100, ncol = 2)
tod = rep(1:24, 5)[1:100]

probs = stateprobs_p(delta, Gamma, allprobs, tod)
```
