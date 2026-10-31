# Estimated parameters of a `LaMaModel` object

Returns the estimated fixed parameters of a model from
[`lama`](https://janolefi.github.io/LaMa/reference/lama.md) or
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md), on the
working scale and in the order of `obj$par`, i.e. without parameters
fixed via `map` and without random effects. Parameters with several
entries are named by their position among the estimated entries, e.g.
`logmu.1` and `logmu.2`, matrices in column-major order as in `obj$par`.

## Usage

``` r
# S3 method for class 'LaMaModel'
coef(object, ...)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- ...:

  additional arguments (not used)

## Value

named numeric vector of the estimated parameters

## See also

[`vcov.LaMaModel`](https://janolefi.github.io/LaMa/reference/vcov.LaMaModel.md)
