# Build all embedded transition probability matrices of an inhomogeneous HSMM

Hidden semi-Markov models are defined in terms of state durations and an
**embedded** transition probability matrix that contains the conditional
transition probabilities given that the **current state is left**. This
matrix necessarily has diagonal entries all equal to zero as
self-transitions are impossible. We can allow this matrix to vary with
covariates, which is the purpose of this function.

It builds all embedded/ conditional transition probability matrices
based on a design and parameter matrix. For each row of the matrix, the
inverse multinomial logistic link is applied.

For a matrix of dimension c(N,N), the number of free off-diagonal
elements is N\*(N-2) which determines the number of rows of the
parameter matrix.

Compatible with automatic differentiation by `RTMB`

## Usage

``` r
tpm_emb_g(Z, beta, report = TRUE)
```

## Arguments

- Z:

  covariate design matrix with or without intercept column, i.e. of
  dimension c(n, p) or c(n, p+1)

  If `Z` has only p columns, an intercept column of ones will be added
  automatically.

- beta:

  matrix of coefficients for the off-diagonal elements of the embedded
  transition probability matrix

  Needs to be of dimension c(N\*(N-2), p+1), where the first column
  contains the intercepts. p can be 0, in which case the model is
  homogeneous.

- report:

  logical, indicating whether the coefficient matrix beta should be
  reported from the fitted model. Defaults to `TRUE`.

## Value

array of embedded/ conditional transition probability matrices of
dimension c(N,N,n)

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_cont()`](https://janolefi.github.io/LaMa/reference/tpm_cont.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
## parameter matrix for 3-state HSMM
beta = matrix(c(rep(0, 3), -0.2, 0.2, 0.1), nrow = 3)
# no intercept
Z = rnorm(100)
omega = tpm_emb_g(Z, beta)
# intercept
Z = cbind(1, Z)
omega = tpm_emb_g(Z, beta)
```
