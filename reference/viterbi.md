# Viterbi algorithm for state decoding in homogeneous HMMs

The Viterbi algorithm allows one to decode the most probable state
sequence of an HMM.

## Usage

``` r
viterbi(delta, Gamma, allprobs, trackID = NULL, mod = NULL)
```

## Arguments

- delta:

  initial distribution of length N, or matrix of dimension c(k,N) for k
  independent tracks, if `trackID` is provided

- Gamma:

  transition probability matrix of dimension c(N,N) or array of
  transition probability matrices of dimension c(N,N,k) if `trackID` is
  provided

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
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md) in
  your likelihood function, the objects needed for state decoding are
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
[`viterbi_g()`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
[`viterbi_p()`](https://janolefi.github.io/LaMa/reference/viterbi_p.md)

## Examples

``` r
delta = c(0.5, 0.5)
Gamma = matrix(c(0.9, 0.1, 0.2, 0.8), nrow = 2, byrow = TRUE)
allprobs = matrix(runif(200), nrow = 100, ncol = 2)
states = viterbi(delta, Gamma, allprobs)
```
