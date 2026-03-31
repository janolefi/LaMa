# Calculate continuous time transition probabilities

A continuous-time Markov chain is described by an infinitesimal
generator matrix \\Q\\. When observing data at time points \\t_1, \dots,
t_n\\ the transition probabilites between \\t_i\\ and \\t\_{i+1}\\ are
caluclated as \$\$\Gamma(\Delta t_i) = \exp(Q \Delta t_i),\$\$ where
\\\exp()\\ is the matrix exponential. The mapping \\\Gamma(\Delta t)\\
is also called the **Markov semigroup**. This function calculates all
transition matrices based on a given generator and time differences.

## Usage

``` r
tpm_ct(Q, timediff, rates = NULL, ad = NULL, report = TRUE)

tpm_cont(Q, timediff, rates = NULL, ad = NULL, report = TRUE)
```

## Arguments

- Q:

  infinitesimal generator matrix of the continuous-time Markov chain of
  dimension `c(nStates, nStates)`

- timediff:

  time differences between observations of length `nObs`-1 when based on
  `nObs` observations

- rates:

  optional vector of state-dependent rates for MM(M)PP fitting. For the
  MM(M)PP likelihood, the matrices needed in the forward algorithm are
  \\\exp((Q - \Lambda) \Delta t)\\, where \\\Lambda\\ is a diagonal
  matrix with the state-dependent rates on the diagonal.

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

- report:

  logical, indicating whether `Q` should be reported from the fitted
  model. Defaults to `TRUE`, but only works if `ad = TRUE`.

## Value

array of continuous-time transition matrices of dimension
`c(nStates, nStates, nObs-1)`

## See also

Other transition probability matrix functions:
[`generator()`](https://janolefi.github.io/LaMa/reference/generator.md),
[`generator_g()`](https://janolefi.github.io/LaMa/reference/generator_g.md),
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md),
[`tpm_emb()`](https://janolefi.github.io/LaMa/reference/tpm_emb.md),
[`tpm_emb_g()`](https://janolefi.github.io/LaMa/reference/tpm_emb_g.md),
[`tpm_g()`](https://janolefi.github.io/LaMa/reference/tpm_g.md),
[`tpm_g2()`](https://janolefi.github.io/LaMa/reference/tpm_g2.md),
[`tpm_p()`](https://janolefi.github.io/LaMa/reference/tpm_p.md)

## Examples

``` r
# building a Q matrix for a 3-state cont.-time Markov chain
Q = generator(rep(-2, 6))

# draw random time differences
timediff = rexp(100, 10)

# compute all transition matrices
Gamma = tpm_cont(Q, timediff)
```
