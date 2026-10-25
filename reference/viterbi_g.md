# Viterbi algorithm for state decoding in inhomogeneous HMMs

The Viterbi algorithm decodes the most probable state sequence of an
HMM.

## Usage

``` r
viterbi_g(delta, Gamma, allprobs, trackID = NULL, mod = NULL)
```

## Arguments

- delta:

  initial distribution; either

  - a vector of length `nStates`, or

  - a matrix of dimension `c(nTracks, nStates)` if `trackID` is provided

- Gamma:

  array of transition probability matrices of dimension
  `c(nStates, nStates, nObs)`, where the first slice of each track is
  ignored as there is no transition into the start of a track. For a
  single track, an array of dimension `c(nStates, nStates, nObs-1)` is
  also accepted.

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
  [`forward_g`](https://janolefi.github.io/LaMa/reference/forward_g.md)
  in the likelihood, these are reported automatically after model
  fitting and the object returned by `obj$report()`,
  [`lama`](https://janolefi.github.io/LaMa/reference/lama.md) or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) can be
  passed directly.

## Value

vector of decoded states of length `nObs`

## Details

**Deprecated:** please use
[`viterbi`](https://janolefi.github.io/LaMa/reference/viterbi.md) for
all model types instead.

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
#> tpm_g() is deprecated, please use tpm(beta, Z) for all model types instead.
allprobs = matrix(runif(20), nrow = 10, ncol = 2)
states = viterbi_g(delta, Gamma[,,-1], allprobs)
#> viterbi_g() is deprecated, please use viterbi() for all model types instead.
```
