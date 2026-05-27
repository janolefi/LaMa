# Calculate conditional local state probabilities in HMMs

Computes \$\$\Pr(\text{State}\_t = j \mid X_1, ..., X_T)\$\$ for a given
HMM.

## Usage

``` r
stateprobs(
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

  initial distribution; either

  - a vector of length `nStates`, or

  - a matrix of dimension `c(nTracks, nStates)` if `trackID` is provided

- Gamma:

  transition probability matrix; either

  - a matrix of dimension `c(nStates, nStates)`,

  - an array of dimension `c(nStates, nStates, nTracks)` if `trackID` is
    provided, or

  - an array of dimension `c(nStates, nStates, nObs)` for time-varying
    transition probabilities, in which case
    [`stateprobs_g`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md)
    is called internally

- allprobs:

  matrix of state-dependent probabilities or density values of dimension
  `c(nObs, nStates)`

- trackID:

  optional vector of length `nObs` containing `nTracks` unique IDs that
  separate tracks

- mod:

  optional model object containing `delta`, `Gamma`, `allprobs`, and
  optionally `trackID`. When using `RTMB::MakeADFun` or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) with
  [`forward`](https://janolefi.github.io/LaMa/reference/forward.md) in
  the likelihood, these are reported automatically after model fitting
  and the object returned by `RTMB::report()` or
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) can be
  passed directly.

- forecast:

  logical, indicating if forecast probabilities \\\Pr(\text{State}\_t =
  j \mid X_1, ..., X\_{t-1})\\ should be calculated instead.

## Value

matrix of conditional state probabilities of dimension
`c(nObs, nStates)`

## See also

Other decoding functions:
[`stateprobs_g()`](https://janolefi.github.io/LaMa/reference/stateprobs_g.md),
[`stateprobs_p()`](https://janolefi.github.io/LaMa/reference/stateprobs_p.md),
[`viterbi()`](https://janolefi.github.io/LaMa/reference/viterbi.md),
[`viterbi_g()`](https://janolefi.github.io/LaMa/reference/viterbi_g.md),
[`viterbi_p()`](https://janolefi.github.io/LaMa/reference/viterbi_p.md)

## Examples

``` r
Gamma = tpm(c(-1,-2))
delta = stationary(Gamma)
allprobs = matrix(runif(10), nrow = 10, ncol = 2)

probs = stateprobs(delta, Gamma, allprobs)
```
