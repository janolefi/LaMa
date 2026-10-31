# Variance-covariance matrix of a `LaMaModel` object

Computes the approximate covariance matrix of the estimated parameters
`coef(object)`, as the inverse Hessian of the negative log-likelihood at
the estimate.

## Usage

``` r
# S3 method for class 'LaMaModel'
vcov(object, ...)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- ...:

  additional arguments (not used)

## Value

covariance matrix with rows and columns named like `coef(object)`

## Details

The Hessian is

- the exact Hessian `obj$he()` for models without random effects,

- for models with random effects, the Hessian of the
  Laplace-approximated marginal likelihood w.r.t. the fixed parameters,
  by finite differences of `obj$gr()` (as in `sdreport`),

- for models from
  [`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md), the
  Hessian of the penalised negative log-likelihood, including the
  penalised coefficients. The covariance is then conditional on the
  estimated smoothing parameters.

If the Hessian is not positive definite, which indicates that the model
has not converged, a warning is issued.

## See also

[`coef.LaMaModel`](https://janolefi.github.io/LaMa/reference/coef.LaMaModel.md)
