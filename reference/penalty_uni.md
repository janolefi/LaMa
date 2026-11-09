# Penalty approximation of unimodality constraints for univariate smooths

Penalty approximation of unimodality constraints for univariate smooths

## Usage

``` r
penalty_uni(coef, m, kappa = 1000, concave = TRUE, rho = 20)
```

## Arguments

- coef:

  coefficient vector or matrix on which to apply the unimodality penalty

- m:

  vector of indices for the position of the coefficient mode. If `coef`
  is a vector, must be of length 1. Otherwise, must be of length equal
  to nrow(coef)

- kappa:

  global scaling factor for the penalty

- concave:

  logical; if `TRUE` (default), the penalty enforces increasing until
  the mode then decreasing. If the coefficients should decrease until
  the mode, then increase, set `concave = FALSE`.

- rho:

  control parameter for smooth approximation to `min(x, 0)` used
  internally. For large values, gets closer to true minimum function but
  less stable.

## Value

a numeric value of the penalty for the given coefficients

## Examples

``` r
## coefficient vector
coef <- c(1, 2, 3, 2, 1)
# mode at position 3
penalty_uni(coef, m = 3) # basically zero
#> [1] 5.152884e-07
#' # mode at position 2
penalty_uni(coef, m = 2) # large positive penalty
#> [1] 1000

## coefficient matrix
coef <- rbind(coef, coef)
m <- c(1, 4)
penalty_uni(coef, m)
#> [1] 3000
```
