# Process and standardise formulas for the state process of hidden Markov models

Process and standardise formulas for the state process of hidden Markov
models

## Usage

``` r
process_hid_formulas(formulas, nStates, ref = NULL)
```

## Arguments

- formulas:

  formulas for the transition process of a hidden Markov model, either
  as a single formula, a list of formulas, or a matrix.

- nStates:

  number of states of the Markov chain

- ref:

  optional vector of reference categories for each state, defaults to
  `1:nStates`. If provided, must be of length `nStates` and contain
  valid state indices. If a formula matrix is provided, this cannot be
  specified because reference categories are specified by one `"."`
  entry in each row.

## Value

named list of formulas of length `nStates * (nStates - 1)`, where each
formula corresponds to a transition from state \\i\\ to state \\j\\,
excluding transitions from `i` to `ref[i]`.

## Examples

``` r
# single formula for all non-reference category elements
formulas = process_hid_formulas(~ s(x), nStates = 3)
# now a list of length 6 with names tr.ij, not including reference categories

# different reference categories
formulas = process_hid_formulas(~ s(x), nStates = 3, ref = c(1,1,1))

# different formulas for different entries (and only for 2 of 6)
formulas = list(tr.12 ~ s(x), tr.23 ~ s(y))
formulas = process_hid_formulas(formulas, nStates = 3, ref = c(1,1,1))
# also a list of length 6, remaining entries filled with tr.ij ~ 1

# matrix input with reference categories
formulas = matrix(c(".", "~ s(x)", "~ s(y)",
                    "~ g", ".", "~ I(x^2)",
                    "~ y", "~ 1", "."), 
                    nrow = 3, byrow = TRUE)
# dots define reference categories
formulas = process_hid_formulas(formulas, nStates = 3)
```
