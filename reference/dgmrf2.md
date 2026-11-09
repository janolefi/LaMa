# Reparametrised multivariate Gaussian distribution

Density function of the multivariate Gaussian distribution
reparametrised in terms of its precision matrix (inverse variance). This
implementation is particularly useful for defining the **joint
log-likelihood** with penalised splines or i.i.d. random effects that
have a multivariate Gaussian distribution with fixed precision/ penalty
matrix \\\lambda S\\. As \\S\\ is fixed and only scaled by \\\lambda\\,
it is more efficient to precompute the determinant of \\S\\ (for the
normalisation constant) and only scale the quadratic form by \\\lambda\\
when multiple spline parameters/ random effects with different
\\\lambda\\'s but the same penalty matrix \\S\\ are evaluated.

## Usage

``` r
dgmrf2(x, mu = 0, S, lambda, logdetS = NULL, log = FALSE)
```

## Arguments

- x:

  density evaluation point, either a vector or a matrix

- mu:

  mean parameter. Either scalar or vector

- S:

  unscaled precision matrix

- lambda:

  precision scaling parameter

  Can be a vector if `x` is a matrix. Then each row of `x` is evaluated
  with the corresponding `lambda`. This is beneficial from an efficiency
  perspective because the determinant of `S` is only computed once.

- logdetS:

  Optional precomputed log determinant of the precision matrix `S`. If
  the precision matrix does not depend on parameters, it can be
  precomputed and passed to the function.

- log:

  logical; if `TRUE`, densities are returned on the log scale.

## Value

vector of density values

## Details

This implementation allows for automatic differentiation with `RTMB`.

## Examples

``` r
x = matrix(runif(30), nrow = 3)

# iid random effects
S = diag(10)
sigma = c(1, 2, 3) # random effect standard deviations
lambda = 1 / sigma^2
d = dgmrf2(x, 0, S, lambda)

# P-splines
L = diff(diag(10), diff = 2) # second-order difference matrix
S = t(L) %*% L
lambda = c(1,2,3)
d = dgmrf2(x, 0, S, lambda, log = TRUE)
```
