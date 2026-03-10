# Calculate conditional local state probabilities for inhomogeneous HMMs

Computes \$\$\Pr(S_t = j \mid X_1, ..., X_T)\$\$ for inhomogeneous HMMs

## Usage

``` r
stateprobs_g(
  delta,
  Gamma,
  allprobs,
  trackID = NULL,
  mod = NULL,
  forecast = FALSE
)
```

## Arguments

- delta:

  initial or stationary distribution of length N, or matrix of dimension
  c(k,N) for k independent tracks, if `trackID` is provided

- Gamma:

  array of transition probability matrices of dimension c(N,N,n-1), as
  in a time series of length n, there are only n-1 transitions

  If an array of dimension c(N,N,n) for a single track is provided, the
  first slice will be ignored.

  If `trackID` is provided, `Gamma` needs to be an array of dimension
  c(N,N,n), where n is the number of rows in `allprobs`. Then for each
  track the first transition matrix will be ignored.

- allprobs:

  matrix of state-dependent probabilities/ density values of dimension
  c(n, N)

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
  [`forward_g`](https://janolefi.github.io/LaMa/reference/forward_g.md)
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
[`stateprobs_p()`](https://janolefi.github.io/LaMa/reference/stateprobs_p.md),
[`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md),
[`viterbi_g()`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
[`viterbi_p()`](https://janolefi.github.io/LaMa/reference/viterbi_p.md)

## Examples

``` r
Gamma = tpm_g(runif(10), matrix(c(-1,-1,1,-2), nrow = 2, byrow = TRUE))
delta = c(0.5, 0.5)
allprobs = matrix(runif(20), nrow = 10, ncol = 2)

probs = stateprobs_g(delta, Gamma[,,-1], allprobs)
```
