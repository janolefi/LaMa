# Viterbi algorithm for state decoding in inhomogeneous HMMs

The Viterbi algorithm allows one to decode the most probable state
sequence of an HMM.

## Usage

``` r
viterbi_g(delta, Gamma, allprobs, trackID = NULL, mod = NULL)
```

## Arguments

- delta:

  initial distribution of length N, or matrix of dimension c(k,N) for k
  independent tracks, if `trackID` is provided

- Gamma:

  array of transition probability matrices of dimension c(N,N,n-1), as
  in a time series of length n, there are only n-1 transitions

  If an array of dimension c(N,N,n) is provided for a single track, the
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

## Value

vector of decoded states of length n

## See also

Other decoding functions:
[`stateprobs()`](https://janolefi.github.io/LaMa/reference/stateprobs.md),
[`stateprobs_g()`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md),
[`stateprobs_p()`](https://janolefi.github.io/LaMa/reference/stateprobs_p.md),
[`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md),
[`viterbi_p()`](https://janolefi.github.io/LaMa/reference/viterbi_p.md)

## Examples

``` r
delta = c(0.5, 0.5)
Gamma = tpm_g(runif(10), matrix(c(-2,-2,1,-1), nrow = 2))
allprobs = matrix(runif(20), nrow = 10, ncol = 2)
states = viterbi_g(delta, Gamma[,,-1], allprobs)
```
