# Plot a `LaMaModel` object

Plots the pseudo-residuals of a model from
[`report`](https://janolefi.github.io/LaMa/reference/report.md) or
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md), computed
by
[`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md)
from the model object alone, with one row of plots per data stream.

## Usage

``` r
# S3 method for class 'LaMaModel'
plot(x, ...)
```

## Arguments

- x:

  model object of class `"LaMaModel"`

- ...:

  further arguments passed to
  [`plot.LaMaResiduals`](https://janolefi.github.io/LaMa/reference/plot.LaMaResiduals.md)

## Value

the pseudo-residuals, invisibly
