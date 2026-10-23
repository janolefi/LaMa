# Summary of a `qremlModel` object

Summarises a model from
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md): the
summary of
[`summary.LaMaModel`](https://janolefi.github.io/LaMa/reference/summary.LaMaModel.md),
plus the effective degrees of freedom and smoothing parameters of the
smooth terms and diagnostics of the smoothness selection.

## Usage

``` r
# S3 method for class 'qremlModel'
summary(object, ...)

# S3 method for class 'summary.qremlModel'
print(x, digits = max(3, getOption("digits") - 3), ...)
```

## Arguments

- object:

  model object of class `"qremlModel"`

- ...:

  additional arguments (not used)

- x:

  object of class `"summary.qremlModel"`

- digits:

  number of significant digits to print

## Value

[`summary()`](https://rdrr.io/r/base/summary.html) returns an object of
class `"summary.qremlModel"`, extending `"summary.LaMaModel"` by `edf`,
`lambda` and `selection`. [`print()`](https://rdrr.io/r/base/print.html)
returns `x` invisibly.
