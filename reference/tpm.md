# Build the transition probability matrix from unconstrained parameter vector

Markov chains are parametrised in terms of a transition probability
matrix \\\Gamma\\, for which each row contains a conditional probability
distribution of the next state given the current state. Hence, each row
has entries between 0 and 1 that need to sum to one.

For numerical optimisation, we parameterise in terms of unconstrained
parameters, thus this function computes said matrix from an
unconstrained parameter vector via the inverse multinomial logistic link
(also known as softmax) applied to each row.

## Usage

``` r
tpm(param, byrow = FALSE, ref = NULL)
```

## Arguments

- param:

  unconstrained parameter vector of length N\*(N-1) where N is the
  number of states of the Markov chain

- byrow:

  logical indicating if the transition probability matrix should be
  filled by row

  Defaults to `FALSE`, but should be set to `TRUE` if one wants to work
  with a matrix of beta parameters returned by popular HMM packages like
  `moveHMM`, `momentuHMM`, or `hmmTMB`.

- ref:

  Optional integer vector of length N giving, for each row, the column
  index of the reference state (its predictor is fixed to 0). Defaults
  to the diagonal (`ref = 1:N`).

## Value

Transition probability matrix of dimension c(N,N)

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/reference/generator.md),
[`tpm_cont()`](https://janolefi.github.io/reference/tpm_cont.md),
[`tpm_emb()`](https://janolefi.github.io/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/reference/tpm_p.md)

## Examples

``` r
# 2 states: 2 free off-diagonal elements
par1 = rep(-1, 2)
Gamma1 = tpm(par1)

# 3 states: 6 free off-diagonal elements
par2 = rep(-2, 6)
Gamma2 = tpm(par2)
```
