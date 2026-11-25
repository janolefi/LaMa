# Takeuchi information criterion

Computes the Takeuchi information criterion (TIC), a version of AIC
whose penalty remains valid if the model is misspecified, for one or
several fitted models.

## Usage

``` r
TIC(object, ...)

# S3 method for class 'LaMaModel'
TIC(
  object,
  ...,
  type = c("HAC", "sandwich", "cluster"),
  kernel = c("QS", "Bartlett"),
  bw = NULL,
  cluster = NULL,
  cr = c("CR1", "CR0", "CR3", "CR3onestep")
)
```

## Arguments

- object:

  model object of class `"LaMaModel"`

- ...:

  further models of class `"LaMaModel"`, compared as in
  [`AIC`](https://rdrr.io/r/stats/AIC.html)

- type:

  robust covariance used for the penalty: `"HAC"` (default),
  `"sandwich"` or `"cluster"`, see
  [`vcov.LaMaModel`](https://janolefi.github.io/LaMa/reference/vcov.LaMaModel.md)

- kernel, bw, cluster, cr:

  passed to
  [`vcov.LaMaModel`](https://janolefi.github.io/LaMa/reference/vcov.LaMaModel.md)

## Value

For one model, the TIC, with the effective number of parameters as
attribute `"edf"`. For several, a data frame with one row per model and
columns `df` (number of parameters), `edf` and `TIC`.

## Details

With \\H\\ the Hessian of the negative log-likelihood and \\J\\ the
variance of the score sum (the meat of the sandwich, see
[`vcov.LaMaModel`](https://janolefi.github.io/LaMa/reference/vcov.LaMaModel.md)),
\$\$\mathrm{TIC} = -2 \log L(\hat\theta) + 2 \\ \mathrm{edf}, \qquad
\mathrm{edf} = \mathrm{tr}(J H^{-1}).\$\$ The penalty of AIC, the number
of parameters \\p\\, assumes \\J = H\\, which holds if the model is
correctly specified. TIC estimates it instead: the effective number of
parameters `edf` approximates the optimism \\E\[(\hat\theta -
\theta^\*)^\top H (\hat\theta - \theta^\*)\] \approx \mathrm{tr}(H \\
\mathrm{Var}(\hat\theta))\\ of the maximised log-likelihood, so models
whose estimates are more variable than the model itself suggests are
penalised more.

It is computed from the robust covariance \\V = H^{-1} J H^{-1}\\ as
\\\mathrm{edf} = \mathrm{tr}(H V)\\, so every robust `type` of
[`vcov.LaMaModel`](https://janolefi.github.io/LaMa/reference/vcov.LaMaModel.md)
can be used, including the jackknives `cr = "CR3"` and
`cr = "CR3onestep"`. The default is the HAC sandwich, as the penalty
needs the correct variance of the estimator, which for serially
dependent scores is not given by the outer product of the scores
(`type = "sandwich"`). For several independent tracks,
`type = "cluster"` is an alternative.

Under correct specification, `edf` is close to \\p\\, and TIC close to
AIC. Values well above \\p\\ indicate misspecification, e.g. serial
dependence the model misses. Values below \\p\\ arise from noise in the
estimated \\J\\, from a model that overstates the variability of the
data, or from parameters that are nearly not identified, e.g. a
superfluous state. In the last case, neither AIC nor TIC is justified,
so neither should be used to choose the number of states. `edf` is
itself estimated and, especially for short series, noisy.

## References

Takeuchi, K. (1976). Distribution of informational statistics and a
criterion of model fitting. *Suri-Kagaku (Mathematical Sciences)*, 153,
12–18 (in Japanese).

Konishi, S. and Kitagawa, G. (2008). *Information Criteria and
Statistical Modeling*. Springer, New York.

## See also

[`vcov.LaMaModel`](https://janolefi.github.io/LaMa/reference/vcov.LaMaModel.md),
[`AIC`](https://rdrr.io/r/stats/AIC.html)
