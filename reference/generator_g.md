# Build generator matrices of a continuous-time Markov chain

This function builds **infinitesimal generator matrices** for a
**continuous-time Markov chain** based on a design matrix and
coefficient matrix.

## Usage

``` r
generator_g(Z, beta, Eta = NULL, byrow = FALSE, report = TRUE)
```

## Arguments

- Z:

  Covariate design matrix with or without intercept column, i.e. of
  dimension `c(nObs, p)` or `c(nObs, p+1)`. If not provided, intercept
  column is added automatically.

- beta:

  Matrix of coefficients for the off-diagonal elements of the generator
  matrix of dimension `c(nStates * (nStates-1), p+1)`. First columns
  contains the intercepts.

- Eta:

  optional pre-calculated matrix of linear predictors of dimension
  `c(nObs, nStates * (nStates-1))`. If provided, no `Z` and `beta` are
  necessary and will be ignored.

- byrow:

  logical indicating if the generator matrices should be filled by row

- report:

  logical, indicating whether the generator matrices `Q` should be
  reported from the fitted model. Defaults to `TRUE`, but only works if
  when automatic differentiation with `RTMB` is used.

## Value

array of infinitesimal generator matrices of dimension
`c(nStates, nStates, nObs)`

## Details

Off-diagonal entries are calculated as \\\exp(Z \beta_i)\\ to ensure
positivity. The diagonal entries are then set to the negative row sums,
which is required for generator matrices.

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
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
