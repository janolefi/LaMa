# Build the design and penalty matrices for smooth density estimation

`smooth_dens_construct()` is deprecated, use
[`smoothdens`](https://janolefi.github.io/LaMa/reference/smoothdens.md)
instead. It returns the same objects, but ordered by component first and
data stream second.

## Usage

``` r
smooth_dens_construct(
  data,
  par,
  type = "real",
  k = 25,
  knots = NULL,
  degree = 3,
  diff_order = 2
)
```

## Arguments

- data:

  named data frame of 1 or multiple data streams

- par:

  list of initial means and sds (`"real"`, `"positive"`) or means and
  concentrations (`"circular"`, called `kappa` or `concentration`), one
  value per state. For several data streams, a named list of such lists.

- type:

  type of each data stream, either `"real"` for data on the reals,
  `"positive"` for data on the positive reals or `"circular"` for
  angular data.

- k:

  number of basis functions for each data stream

- knots:

  optional vector of knots (including the boundary knots) to be used for
  basis construction, for several data streams a named list of such
  vectors. If not provided, the knots are placed equidistantly for
  `"real"` and `"circular"` and using polynomial spacing for
  `"positive"`.

  For `"real"` and `"positive"` `k - degree + 1` knots are needed, for
  `"circular"` `k + 1` knots are needed.

- degree:

  degree of the B-spline basis functions for each data stream, defaults
  to cubic B-splines

- diff_order:

  order of differencing used for the P-Spline penalty matrix for each
  data stream. Defaults to second-order differences.

## Value

a nested list containing the design matrices `Z`, the penalty matrices
`S`, the initial coefficients `coef` the prediction design matrices
`Z_predict`, the prediction grids `xseq`, and details for the basis
expansion for each data stream.

## Examples

``` r
data = data.frame(x = rnorm(100))
# deprecated:
# SmoothDens = smooth_dens_construct(data, par = list(x = list(mean = 0, sd = 1)))
# use instead:
dens = smoothdens(data$x, par = list(mean = 0, sd = 1))
```
