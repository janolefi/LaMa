# Build the design and the penalty matrix for models involving penalised splines based on a formula and a data set

Build the design and the penalty matrix for models involving penalised
splines based on a formula and a data set

## Usage

``` r
make_matrices(formula, data, knots = NULL)
```

## Arguments

- formula:

  formula as used in `mgcv`. Formulas can be right-side only, or contain
  a response variable, which is just extracted for naming.

  Can also be a list of formulas, which are then processed separately.
  In that case, both a named list of right-side only formulas or a list
  of formulas with response variables can be provided.

- data:

  data frame containing all the variables on the right side of the
  formula(s)

- knots:

  optional list containing user specified knot values for each covariate
  to be used for basis construction. For most bases the user simply
  supplies the `knots` to be used, which must match up with the `k`
  value supplied (note that the number of knots is not always just `k`).
  See `mgcv` documentation for more details.

  If `formula` is a list, this needs to be a named (based on the
  response variables) list over such lists.

## Value

a list of class `LaMa_matrices` containing:

- `Z`:

  design matrix (or list of such matrices if `formula` is a list))

- `S`:

  list of penalty matrices (with names based on the response terms of
  the formulas as well as the smooth terms and covariates). For
  tensorproduct smooths, corresponding entries are themselves lists,
  containing the \\d\\ marginal penalty matrices if \\d\\ is the
  dimension of the tensor product)

- `pardim`:

  list of parameter dimensions (fixed and penalised separately) for each
  formula, for ease of setting up initial parameters

- `coef`:

  list of coefficient vectors filled with zeros of the correct length
  for each formula, for ease of setting up initial parameters

- `sp0`:

  named vector of **initial penalty strengths**, with one entry per
  penalty matrix and names matching `S` (a smooth with several
  penalties, e.g. a tensor product, contributes one entry each, named
  `<smooth>.<margin>`).

  These are
  [`initial.sp`](https://rdrr.io/pkg/mgcv/man/initial.sp.html)'s guesses
  times 100, capped at `1e6`: the inner optimisation is better behaved
  when it starts constrained, and the outer iteration relaxes the
  penalty from there.

  It can be passed straight to
  [`penalty`](https://janolefi.github.io/LaMa/reference/penalty.md) or
  [`penalty2`](https://janolefi.github.io/LaMa/reference/penalty2.md)
  when each smooth appears once. When a smooth is replicated, e.g. one
  spline per state, use
  [`rep_sp`](https://janolefi.github.io/LaMa/reference/rep_sp.md), as a
  single [`rep()`](https://rdrr.io/r/base/rep.html) over `sp0` is only
  correct for some models.

- `formula`:

  the formula(s) as supplied, in the same nesting as `Z`

- `data`:

  the data frame used for the model(s)

- `gam`:

  unfitted [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) object
  used for construction of `Z` and `S` (or list of such objects if
  `formula` is a list). Its `X` element is not kept, as it is returned
  as `Z`

- `gam0`:

  fitted [`mgcv::gam`](https://rdrr.io/pkg/mgcv/man/gam.html) which is
  used internally to create prediction design matrices (or list of such
  objects if `formula` is a list)

- `knots`:

  knot list used in the basis construction (or named list over such
  lists if `formula` is a list)

## See also

[`predict.LaMa_matrices`](https://janolefi.github.io/LaMa/reference/predict.LaMa_matrices.md)
for prediction design matrix construction based on the model matrices
object created by this function.

## Examples

``` r
data = data.frame(x = runif(100), 
                  y = runif(100),
                  g = factor(rep(1:10, each = 10)))

# unvariate thin plate regression spline
modmat = make_matrices(~ s(x), data)
# univariate P-spline
modmat = make_matrices(~ s(x, bs = "ps"), data)
# adding random intercept
modmat = make_matrices(~ s(g, bs = "re") + s(x, bs = "ps"), data)
# tensorproduct of x and y
modmat = make_matrices(~ s(x) + s(y) + ti(x,y), data)
# multiple formulas at once
modmat = make_matrices(list(mu ~ s(x) + y, sigma ~ s(g, bs = "re")), data = data)
```
