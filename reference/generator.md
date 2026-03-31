# Build the generator matrix of a continuous-time Markov chain

This function builds the **infinitesimal generator matrix** for a
**continuous-time Markov chain** from an unconstrained parameter vector.

## Usage

``` r
generator(
  beta,
  Z = NULL,
  Eta = NULL,
  byrow = FALSE,
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

  logical indicating if the generator matrix should be filled by row

- report:

  logical, indicating whether the generator matrix `Q` should be
  reported from the fitted model. Defaults to `TRUE`, but only works if
  when automatic differentiation with `RTMB` is used.

- param:

  depricated, please use argument `beta` instead.

## Value

infinitesimal generator matrix of dimension `c(nStates, nStates)` or
array of such matrices of dimension `c(nStates, nStates, nObs)` if `Z`
or `Eta` is provided.

## Details

Off-diagonal entries are calculated as \\\exp(\beta_i)\\ to ensure
positivity. The diagonal entries are then set to the negative row sums,
which is required for generator matrices.

If a design matrix `Z` or a matrix of linear predictors `Eta` is
provided, the function will automatically call
[`generator_g`](https://janolefi.github.io/LaMa/reference/generator_g.md)
to build the generator matrix based on the design matrix and coefficient
matrix. In that case, the argument `beta` needs to be a matrix of
coefficients of dimension `c(nStates * (nStates-1), p+1)`, where the
first column contains the intercepts.

## See also

Other transition probability matrix functions:
[`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_ct()`](https://janolefi.github.io/LaMa/reference/tpm_ct.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# 2 states: 2 free off-diagonal elements
generator(rep(-1, 2))
#>            S1         S2
#> S1 -0.3678794  0.3678794
#> S2  0.3678794 -0.3678794
# 3 states: 6 free off-diagonal elements
generator(rep(-2, 6))
#>            S1         S2         S3
#> S1 -0.2706706  0.1353353  0.1353353
#> S2  0.1353353 -0.2706706  0.1353353
#> S3  0.1353353  0.1353353 -0.2706706
```
