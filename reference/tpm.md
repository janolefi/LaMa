# Build the transition probability matrix from unconstrained parameter vector

Markov chains are parametrised in terms of a transition probability
matrix \\\Gamma\\, for which each row contains a conditional probability
distribution of the next state given the current state. Hence, each row
has entries between 0 and 1 that need to sum to one.

For numerical optimisation, we parameterise in terms of unconstrained
parameters, thus this function computes said matrix from an
unconstrained parameter vector via the inverse multinomial logistic link
(also known as softmax) applied to each row.

## Usage

``` r
tpm(
  beta,
  Z = NULL,
  Eta = NULL,
  byrow = FALSE,
  ref = NULL,
  ad = NULL,
  report = TRUE,
  param = NULL
)
```

## Arguments

- beta:

  parameters; either

  - a vector of length `nStates * (nStates-1)`, or

  - a matrix of dimension `c(nStates * (nStates-1), p+1)` if design
    matrix `Z` is also provided.

- Z:

  optional covariate design matrix with or without intercept column,
  i.e. of dimension `c(nObs, p)` or `c(nObs, p+1)`. If provided, `beta`
  needs to be a matrix of dimension `c(nStates * (nStates-1), p+1)`.

- Eta:

  optional pre-calculated matrix of linear predictors of dimension
  `c(nObs, nStates * (nStates-1))`. If provided, `Z` and `beta` will be
  ignored.

- byrow:

  logical indicating if each transition probability matrix should be
  filled by row. Defaults to `FALSE`, but should be set to `TRUE` if one
  wants to work with a matrix of beta parameters returned by popular HMM
  packages like `moveHMM`, `momentuHMM`, or `hmmTMB`.

- ref:

  optional integer vector of length `nStates` giving, for each row, the
  column index of the reference state (its predictor is fixed to 0).
  Defaults to the diagonal (`ref = 1:nStates`).

- ad:

  logical; whether to use automatic differentiation. Determined
  automatically — for debugging only.

- report:

  logical; if `TRUE` (default), `delta`, `Gamma`, `allprobs`, and
  `trackID` are reported from the fitted model. Requires `ad = TRUE`.

- param:

  depricated, please use argument `beta` instead.

## Value

Transition probability matrix of dimension `c(nStates, nStates)` or
array of such matrices of dimension `c(nStates, nStates, nObs)` if `Z`
or `Eta` is provided.

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md),
[`tpm_ct()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# 2 states: 2 free off-diagonal elements
par1 = rep(-1, 2)
Gamma1 = tpm(par1)

# 3 states: 6 free off-diagonal elements
par2 = rep(-2, 6)
Gamma2 = tpm(par2)
```
