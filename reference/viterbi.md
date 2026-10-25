# Viterbi algorithm for state decoding in HMMs

The Viterbi algorithm decodes the most probable state sequence of an
HMM.

## Usage

``` r
viterbi(delta, Gamma, allprobs, trackID = NULL, mod = NULL, tod = NULL)
```

## Arguments

- delta:

  initial distribution; either

  - a vector of length `nStates`,

  - a matrix of dimension `c(nTracks, nStates)` if `trackID` is
    provided, or

  - a model object as for `mod`, so that `mod` can be passed as the
    first argument

- Gamma:

  transition probability matrix; either

  - a matrix of dimension `c(nStates, nStates)`,

  - an array of dimension `c(nStates, nStates, nTracks)` if `trackID` is
    provided,

  - an array of dimension `c(nStates, nStates, nObs)` for time-varying
    transition probabilities, or

  - an array of dimension `c(nStates, nStates, L)` for periodic
    transition probabilities, if `tod` is provided

- allprobs:

  matrix of state-dependent probabilities or density values of dimension
  `c(nObs, nStates)`

- trackID:

  optional vector of length `nObs` containing `nTracks` unique IDs that
  separate tracks. The observations of each track need to be
  consecutive.

- mod:

  optional model object containing `delta`, `Gamma`, `allprobs`, and
  optionally `trackID`. When using `RTMB::MakeADFun` or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) with
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md) in
  the likelihood, these are reported automatically after model fitting
  and the object returned by `RTMB::lama()` or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) can be
  passed directly.

- tod:

  optional integer vector of length `nObs` with values in 1, ..., L,
  indexing the slices of `Gamma` for periodic models

## Value

vector of decoded states of length `nObs`

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
