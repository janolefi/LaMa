# Periodically stationary distribution of a periodically inhomogeneous Markov chain

Computes the periodically stationary distribution of a periodically
inhomogeneous Markov chain.

## Usage

``` r
stationary_p(Gamma, t = NULL, ad = NULL)
```

## Arguments

- Gamma:

  array of transition probability matrices of dimension c(N,N,L)

- t:

  integer index of the time point in the cycle, for which to calculate
  the stationary distribution

  If `t` is not provided, the function calculates all stationary
  distributions for each time point in the cycle.

- ad:

  optional logical, indicating whether automatic differentiation with
  `RTMB` should be used. By default, the function determines this
  itself.

## Value

either the periodically stationary distribution at time t or all
periodically stationary distributions.

## Details

If the transition probability matrix of an inhomogeneous Markov chain
varies only periodically (with period length \\L\\), it converges to a
so-called periodically stationary distribution. This happens, because
the thinned Markov chain, which has a full cycle as each time step, has
homogeneous transition probability matrix \$\$\Gamma_t = \Gamma^{(t)}
\Gamma^{(t+1)} \dots \Gamma^{(t+L-1)}\$\$ for all \\t = 1, \dots, L.\\
The stationary distribution for time \\t\\ satifies \\\delta^{(t)}
\Gamma_t = \delta^{(t)}\\.

This function calculates said periodically stationary distribution.

## References

When using this function, please cite:

Koslik, J. O., Feldmann, C. C., Mews, S., Michels, R., & Langrock, R.
(2025). Inference on the state process of periodically inhomogeneous
hidden Markov models for animal behavior. The Annals of Applied
Statistics, 19(4), 2724–2737.

## See also

[`tpm_p`](https://janolefi.github.io/LaMa/reference/tpm_p.md) and
[`tpm_g`](https://janolefi.github.io/LaMa/reference/tpm_g.md) to create
multiple transition matrices based on a cyclic variable or design matrix

Other stationary distribution functions:
[`stationary()`](https://janolefi.github.io/LaMa/reference/stationary.md),
[`stationary_ct()`](https://janolefi.github.io/LaMa/reference/stationary_ct.md)

## Examples

``` r
# setting parameters for trigonometric link
beta = matrix(c(-1, 2, -1, -2, 1, -1), nrow = 2, byrow = TRUE)
Gamma = tpm_p(beta = beta, degree = 1)
# periodically stationary distribution for specific time point
delta = stationary_p(Gamma, 4)

# all periodically stationary distributions
Delta = stationary_p(Gamma)
```
