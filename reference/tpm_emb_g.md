# Build all embedded transition probability matrices of an inhomogeneous HSMM

**Deprecated:** please use
[`tpm_switch`](https://janolefi.github.io/LaMa/reference/tpm_switch.md)
with argument `Z` instead, i.e. `tpm_switch(beta, Z)`.

## Usage

``` r
tpm_emb_g(Z, beta, report = TRUE)
```

## Arguments

- Z:

  covariate design matrix with or without intercept column, i.e. of
  dimension `c(nObs, p)` or `c(nObs, p+1)`

- beta:

  matrix of coefficients of dimension `c(nStates * (nStates-2), p+1)`,
  with the intercepts in the first column

- report:

  logical, indicating whether the coefficient matrix beta should be
  reported from the fitted model. Defaults to `TRUE`.

## Value

array of embedded/ conditional transition probability matrices of
dimension `c(nStates, nStates, nObs)`

## Examples

``` r
beta = matrix(c(rep(0, 3), -0.2, 0.2, 0.1), nrow = 3)
omega = tpm_emb_g(cbind(1, rnorm(100)), beta)
#> tpm_emb_g() is deprecated, please use tpm_switch(beta, Z) instead.
```
