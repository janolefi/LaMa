# Build a standardised P-Spline design matrix and the associated P-Spline penalty matrix

This function builds the B-spline design matrix for a given data vector.
Importantly, the B-spline basis functions are normalised such that the
integral of each basis function is 1, hence this basis can be used for
spline-based density estimation, when the basis functions are weighted
by non-negative weights summing to one.

## Usage

``` r
make_matrices_dens(
  x,
  k,
  type = "real",
  degree = 3,
  knots = NULL,
  diff_order = 2,
  pow = 0.5,
  npoints = 10000
)
```

## Arguments

- x:

  data vector

- k:

  number of basis functions

- type:

  type of the data, either `"real"` for data on the reals, `"positive"`
  for data on the positive reals or `"circular"` for circular data like
  angles.

- degree:

  degree of the B-spline basis functions, defaults to cubic B-splines

- knots:

  optional vector of knots (including the boundary knots) to be used for
  basis construction. If not provided, the knots are placed
  equidistantly for `"real"` and `"circular"` and using polynomial
  spacing for `"positive"`.

  For `"real"` and `"positive"` `k - degree + 1` knots are needed, for
  `"circular"` `k + 1` knots are needed.

- diff_order:

  order of differencing used for the P-Spline penalty matrix. Defaults
  to second-order differences.

- pow:

  power for polynomial knot spacing, only used for `type = "positive"`

- npoints:

  number of points used in the numerical integration for normalizing the
  B-spline basis functions

## Value

list containing the design matrix `Z`, the penalty matrix `S`, the
prediction design matrix `Z_predict`, the prediction grid `xseq`, and
details for the basis expansion.

The penalty matrix leaves out the last coefficient, which is fixed at
zero for identifiability.

## Examples

``` r
set.seed(1)
# real-valued
x <- rnorm(100)
modmat <- make_matrices_dens(x, k = 20)
# positive-continuouos
x <- rgamma2(100, mean = 5, sd = 2)
modmat <- make_matrices_dens(x, k = 20, type = "positive")
# circular
x <- rvm(100, mu = 0, kappa = 2)
modmat <- make_matrices_dens(x, k = 20, type = "circular")
# bounded in an interval
x <- rbeta(100, 1, 2)
modmat <- make_matrices_dens(x, k = 20)
```
