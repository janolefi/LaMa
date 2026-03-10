# Build all transition probability matrices of an inhomogeneous HMM

In an HMM, we often model the influence of covariates on the state
process by linking them to the transition probabiltiy matrix. Most
commonly, this is done by specifying linear predictors \$\$
\eta\_{ij}^{(t)} = \beta^{(ij)}\_0 + \beta^{(ij)}\_1 z\_{t1} + \dots +
\beta^{(ij)}\_p z\_{tp} \$\$ for each off-diagonal element (\\i \neq
j\\) of the transition probability matrix and then applying the inverse
multinomial logistic link (also known as softmax) to each row. This
function efficiently calculates all transition probabilty matrices for a
given design matrix `Z` and parameter matrix `beta`.

## Usage

``` r
tpm_g2(Z, beta, byrow = FALSE, ad = NULL, report = TRUE, ref = NULL)
```

## Arguments

- Z:

  covariate design matrix with or without intercept column, i.e. of
  dimension c(n, p) or c(n, p+1)

  If `Z` has only p columns, an intercept column of ones will be added
  automatically.

  Can also be a list of N\*(N-1) design matrices with different number
  of columns but the same number of rows. In that case, no intercept
  column will be added.

- beta:

  matrix of coefficients for the off-diagonal elements of the transition
  probability matrix

  Needs to be of dimension c(N\*(N-1), p+1), where the first column
  contains the intercepts.

  If `Z` is a list, `beta` can also be a list of length N\*(N-1) with
  each entry being a vector or a (long) matrix of coefficients, each
  matching the dimension of the corresponding entry in `Z`.

- byrow:

  logical indicating if each transition probability matrix should be
  filled by row

  Defaults to `FALSE`, but should be set to `TRUE` if one wants to work
  with a matrix of beta parameters returned by popular HMM packages like
  `moveHMM`, `momentuHMM`, or `hmmTMB`.

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

- report:

  logical, indicating whether the coefficient matrix `beta` should be
  reported from the fitted model. Defaults to `TRUE`, but only works if
  `ad = TRUE`.

- ref:

  optional vector of length N with the reference state indices for each
  column of the transition probability matrix. Each row in the
  transition matrix corresponds to a multinomial regression, hence one
  state needs to be the reference category. Defaults to off-diagonal
  elements (`ref = 1:N`).

## Value

array of transition probability matrices of dimension c(N,N,n)

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_cont()`](https://janolefi.github.io/LaMa/reference/tpm_cont.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
Z = matrix(runif(200), ncol = 2)
beta = matrix(c(-1, 1, 2, -2, 1, -2), nrow = 2, byrow = TRUE)
Gamma = tpm_g(Z, beta)
```
