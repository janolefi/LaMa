# Repeat initial penalty strengths for replicated smooths

Expands the `sp0` element of a
[`modmat`](https://janolefi.github.io/LaMa/reference/modmat.md) object
into a `sp` vector for
[`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md) or
[`penalty2`](https://janolefi.github.io/LaMa/reference/penalty2.md),
when each smooth is used more than once, for example one spline per
state or per off-diagonal element of the transition probability matrix.

## Usage

``` r
rep_sp(model_matrices, N)
```

## Arguments

- model_matrices:

  object of class `LaMa_matrices` as returned by
  [`modmat`](https://janolefi.github.io/LaMa/reference/modmat.md)

- N:

  number of replicates. Either a single number used for every smooth, or
  one number per smooth, i.e. per element of `model_matrices$S`

## Value

named numeric vector of initial penalty strengths, of the length and in
the order expected for `sp`

## Details

`sp` runs over replicates **within** each smooth, so for a tensor
product the margins vary fastest. A single
[`rep()`](https://rdrr.io/r/base/rep.html) over `sp0` therefore only
gives the right answer for some models, silently producing the right
length but the wrong values for others. This function repeats each
smooth's penalty strengths separately, which is always correct.

## See also

[`modmat`](https://janolefi.github.io/LaMa/reference/modmat.md), which
computes `sp0`

## Examples

``` r
modmat = modmat(~ s(x) + ti(x, y), data.frame(x = runif(100), y = runif(100)))
modmat$sp0         # one value per penalty
#>      s(x) ti(x,y).x ti(x,y).y 
#>  58.74017 144.93454 161.40577 
rep_sp(modmat, 2)  # ready to use as 'sp' for two replicates of each smooth
#>      s(x).1      s(x).2 ti(x,y).x.1 ti(x,y).y.1 ti(x,y).x.2 ti(x,y).y.2 
#>    58.74017    58.74017   144.93454   161.40577   144.93454   161.40577 
```
