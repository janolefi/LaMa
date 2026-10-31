# Pseudo-residuals of a `LaMaModel` object

Computes the pseudo-residuals of a model from
[`lama`](https://janolefi.github.io/LaMa/reference/lama.md) or
[`qreml`](https://janolefi.github.io/LaMa/reference/qreml.md) from the
model object alone, by calling `pseudo_res(mod = object)`.

## Usage

``` r
# S3 method for class 'LaMaModel'
residuals(object, ...)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- ...:

  additional arguments (not used)

## Value

pseudo-residuals of class `LaMaResiduals`, see
[`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md)

## See also

[`pseudo_res`](https://janolefi.github.io/LaMa/reference/pseudo_res.md),
[`plot.LaMaResiduals`](https://janolefi.github.io/LaMa/reference/plot.LaMaResiduals.md)
