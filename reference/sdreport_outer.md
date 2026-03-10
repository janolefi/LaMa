# Report uncertainty of the estimated smoothing parameters or variances

Computes standard deviations for the smoothing parameters of a model
object returned by `qreml` using the delta method.

## Usage

``` r
sdreport_outer(mod, invert = FALSE)
```

## Arguments

- mod:

  model objects as returned by
  [`qreml`](https://janolefi.github.io/reference/qreml.md)

- invert:

  optional logical; if `TRUE`, the inverse smoothing paramaters
  (variances) are returned along with the transformed standard
  deviations obtained via the delta method.

## Value

list containing `report` matrix summarising parameters and standard
deviations as well as the outer `Hessian` matrix.

## Details

The computations are based on the approximate gradient of the restricted
log likelihood. The outer Hessian is computed by finite differencing of
this gradient. If the inverse smoothing parameters are requested, the
standard deviations are transformed to the variances using the delta
method.

## Examples

``` r
## no examples
```
