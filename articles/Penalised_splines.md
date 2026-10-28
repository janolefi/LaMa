# 5 Penalised splines

> Before diving into this vignette, we recommend reading the vignettes
> [**Introduction to
> LaMa**](https://janolefi.github.io/LaMa/articles/Intro_to_LaMa.html),
> [**Automatic differentiation via
> RTMB**](https://janolefi.github.io/LaMa/articles/LaMa_and_RTMB.html)
> and
> [**Extensions**](https://janolefi.github.io/LaMa/articles/Extensions.html).

This vignette explores how `LaMa` can be used to fit models involving
nonparametric components, represented by **penalised splines**. The main
idea here is that it may be useful to represent some relationships in
our model by smooth functions for which the functional form is not
pre-specified, but flexibly estimated from the data. For HMM-like
models, this is particularly valuable, as the latent nature of the state
process makes modelling choices more difficult. For example, choosing an
appropriate parametric family for the state-dependent distributions may
be difficult as we cannot do state-specific EDA before fitting the
model. Also very difficult is the dependence of transition probabilities
on covariates as the transitions are not directly observed. Hence, the
obvious alternative is to model these kind of relationships flexibly
using splines but imposing a penalty on the smoothness of the estimated
functions. This leads us to **penalised splines**.

`LaMa` contains helper functions that build **design** and **penalty
matrices** for given formulas (using `mgcv` under the hood) and also
functions to estimate models involving penalised splines in a random
effects framework. For the latter to work, the **penalised negative
log-likelihood** needs to be compatible with the R package `RTMB` to
allow for **automatic differentiation (AD)**. For more information on
`RTMB`, see the vignette *LaMa and RTMB* or check out its
[documentation](https://CRAN.R-project.org/package=RTMB). For more
information on penalised splines, we recommend Wood
([2017](#ref-wood2017generalized)).

## Smooth transition probabilities

We will start by investigating a 2-state HMM for the `trex` data set,
containing hourly step lengths and turning angles of a Tyrannosaurus rex
living 66 million years ago. The transition probabilities are modelled
as smooth functions of the time of day using **cyclic P-Splines**. The
relationship can be summarised as

\text{logit}(\gamma\_{ij}^{(t)}) = \beta_0^{(ij)} + s\_{ij}(t), where
s\_{ij}(t) is a smooth periodic function of time of day. We model the
T-rex’s **step lengths** and **turning angles** using state-dependent
**gamma** and **von Mises** distributions.

To ease with model specification, `LaMa` provides the function
[`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md) which
creates **design** and **penalty** matrices for regression settings
based on the R package `mgcv`. The user only needs to specify the right
side of a formula using `mgcv` syntax and provide data. Here, we use
`s(tod, bs = "cp")` to create the matrices for cyclic P-splines (`cp`).
This results in a cubic B-Spline basis, that is wrapped at boundary of
the support (0 and 24). We then append both resulting matrices to the
`dat` list.

``` r

library(LaMa)
#> Loading required package: RTMB

head(trex)
#>   tod      step     angle state
#> 1   9 0.3252437        NA     1
#> 2  10 0.2458265  2.234562     1
#> 3  11 0.2173252 -2.262418     1
#> 4  12 0.5114665 -2.958732     1
#> 5  13 0.3828494  1.811840     1
#> 6  14 0.4220099  1.834668     1
```

``` r

modmat = modmat(~ s(tod, bs = "cp"), # formula
                       data = data.frame(tod = 1:24), # data
                       knots = list(tod = c(0, 24))) # where to wrap the cyclic basis
modmat # formula and initial penalty strength
#> Model matrices for
#> 
#>   ~ s(tod, bs = "cp")
#> 
#> ---
#> 
#> Initial penalty strengths:
#> 
#> s(tod) 
#>  242.8 
#> 
#> Use 'rep_sp()' if a smooth is used more than once.
Z = modmat$Z # spline design matrix
S = modmat$S # penalty matrix
```

We can now specify the **penalised negative log-likelihood function**.
The transition probability matrix can be computed the regular way using
[`tpm()`](https://janolefi.github.io/LaMa/reference/tpm.md). In the last
line we need to add the curvature penalty based on `S`, which we can
conveniently do using
[`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md).

``` r

pnll = function(par) {
  getAll(par, dat)
  # cbinding intercept and spline coefs, because intercept is not penalised
  Gamma = tpm(cbind(beta0, betaSpline), Z)
  # computing all periodically stationary distributions for easy access later
  Delta = stationary_p(Gamma)[c(2:24, 1), ]; REPORT(Delta) # row k for tod k
  # parameter transformations
  mu = exp(logmu); REPORT(mu)
  sigma = exp(logsigma); REPORT(sigma)
  kappa = exp(logkappa); REPORT(kappa)
  # calculating all state-dependent densities
  allprobs = matrix(1, nrow = length(step), ncol = N)
  ind = which(!is.na(step) & !is.na(angle)) # only for non-NA obs.
  for(j in 1:N){
    allprobs[ind,j] = dgamma2(step[ind],mu[j],sigma[j]) * dvm(angle[ind],0,kappa[j])
  }
  -forward(Delta[tod[1],], Gamma, allprobs, tod = tod) + # regular forward algorithm
    penalty(betaSpline, S, sp) # this does all the penalisation work
}
```

We also have to append an `sp` vector to our `dat` list which is the
**initial penalty strength** parameter vector. In this case it is of
length two because our coefficient matrix has two rows.
[`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md)
provides one initial value per penalty in `modmat$sp0`, based on
`mgcv`’s initial guess, and
[`rep_sp()`](https://janolefi.github.io/LaMa/reference/rep_sp.md)
repeats it for each row of the coefficient matrix.

If you are wondering why `sp` is not added to the `par` list, this is
because for penalised likelihood estimation, it is a **hyperparameter**,
hence not a true parameter in the sense of the other parameters in
`par`. One could, at his point, just use the above penalised likelihood
function to do penalised ML for a **fixed** penalty strength `sp`.

``` r

par = list(logmu = log(c(0.3, 2.5)), # state-dependent mean step
           logsigma = log(c(0.3, 2)), # state-dependent sd step
           logkappa = log(c(0.1, 2)), # state-dependent concentration angle
           beta0 = c(-2, -2), # state process intercepts
           betaSpline = matrix(0, nrow = 2, ncol = ncol(Z)-1)) # spline coefs

dat = list(step = trex$step, # observed steps
           angle = trex$angle, # observed angle
           N = 2, # number of states
           tod = trex$tod, # time of day (used for indexing)
           Z = Z, # spline design matrix
           S = S, # penalty matrix
           sp = rep_sp(modmat, 2)) # initial penalty strength, one per row of betaSpline
```

The model fit with automatic smoothness selection can then be conducted
by using the
[`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) function
contained in `LaMa`. The **quasi restricted maximum likelihood**
algorithm finds a good penalty strength parameter `sp` by treating the
spline coefficients as random effects. Under the hood,
[`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) also
constructs an AD function with `RTMB` but uses the **qREML** algorithm
([Koslik 2024](#ref-koslik2024efficient)) to fit the model. We have to
tell the [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
function which parameters are spline coefficients by providing the name
of the corresponding list element of `par`.

There are some rules to follow when using
[`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md):

1.  The likelihood function needs to be `RTMB`-compatible, i.e. have the
    same structure as the likelihood functions in the vignette *LaMa and
    RTMB*. Most importantly, it should be a function of the parameter
    list only.

2.  The penalty strength vector `sp` needs its length to correspond to
    the *total* number of spline coefficient vectors used. In our case,
    this is the number of rows of `betaSpline`, but if we additionally
    had a different spline coefficient (with a different name) in our
    parameter list, possibly with a different length and a different
    penalty matrix, we would have needed more elements in `sp`.

3.  The
    [`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md)
    function can only be called *once* in the likelihood. If several
    spline coefficients are penalised,
    [`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md)
    expects a list of coefficient matrices or vectors and a list of
    penalty matrices. This is shown in the third example in this
    vignette.

&nbsp;

4.  [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
    looks for the initial penalty strengths in `dat` under the name
    `sp`, and falls back to `lambda` (the former default, which still
    works). If `dat` contains `sp`, the name `lambda` is free for other
    things, e.g. the mean of a Poisson distribution. If you want to use
    a different name for the penalty strengths, you have to specify it
    as a character string in the
    [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md) call
    using the `spname` argument.

``` r

system.time(
  mod1 <- qreml(pnll, par, dat, random = "betaSpline")
)
#> outer 1 - sp: 24.279 24.279 
#> outer 2 - sp: 2.428 2.428 
#> outer 3 - sp: 0.243 0.243 
#> outer 4 - sp: 0.378 0.076 
#> outer 5 - sp: 0.312 0.104 
#> outer 6 - sp: 0.311 0.107 
#> outer 7 - sp: 0.309 0.11 
#> outer 8 - sp: 0.309 0.11 
#> outer 9 - sp: 0.309 0.11 
#> outer 10 - sp: 0.309 0.11 
#> outer 11 - sp: 0.309 0.11
#>    user  system elapsed 
#>   7.834   0.082   7.916
```

The `mod` object is now a list that contains everything that is reported
by the likelihood function, but also the `RTMB` object created in the
process. After fitting the model, we can use
[`predict()`](https://rdrr.io/r/stats/predict.html) with the `modmat`
object we created earlier to build a new interpolating design matrix
using the exact same basis expansion specified above. This allows us to
plot the estimated transition probabilities as a smooth function of time
of day. We use
[`MCreport()`](https://janolefi.github.io/LaMa/reference/MCreport.md) to
draw samples from the approximate normal distribution of the MLE and
propagate uncertainty into the plotted quantities.

``` r

color = LaMaColors(2)

tod_seq = seq(0, 24, length = 100)
Z_p = predict(modmat, data.frame(tod = tod_seq))

# MC samples: with report = TRUE, REPORT()ed quantities are sampled directly
samples = MCreport(mod1$obj, report = TRUE)
#> Evaluating Hessian...
#> Computing reported quantities...

# transition probability samples on the plotting grid
g12 = sapply(samples$beta, function(b) tpm(b, Z_p)[1,2,])
g21 = sapply(samples$beta, function(b) tpm(b, Z_p)[2,1,])
g12.q = apply(g12, 1, quantile, probs = c(0.025, 0.975))
g21.q = apply(g21, 1, quantile, probs = c(0.025, 0.975))

# Delta[,2] samples available directly via report = TRUE
Delta2.samples = sapply(samples$Delta, function(d) d[,2])
Delta2.q = apply(Delta2.samples, 1, quantile, probs = c(0.025, 0.975))

# MLE quantities
Gamma_plot = tpm(mod1$beta, Z_p)
Delta = mod1$Delta

oldpar = par(mfrow = c(1, 2))

# transition probabilities with 95% confidence bands
plot(NA, bty = "n", ylim = c(0, 1), xlim = c(0, 24),
     xlab = "time of day", ylab = "transition probability")
polygon(c(tod_seq, rev(tod_seq)), c(g12.q[1,], rev(g12.q[2,])),
        col = adjustcolor(color[1], 0.3), border = NA)
lines(tod_seq, Gamma_plot[1,2,], lwd = 2, col = color[1])
polygon(c(tod_seq, rev(tod_seq)), c(g21.q[1,], rev(g21.q[2,])),
        col = adjustcolor(color[2], 0.3), border = NA)
lines(tod_seq, Gamma_plot[2,1,], lwd = 2, col = color[2])
legend("topleft", lwd = 2, col = color, bty = "n",
       legend = c(expression(gamma[12]^(t)), expression(gamma[21]^(t))))

# periodically stationary distribution with 95% confidence band
plot(Delta[,2], type = "b", lwd = 2, pch = 16, ylim = c(0, 1),
     xlab = "time of day", ylab = "Pr(active)",
     col = "deepskyblue", bty = "n", xaxt = "n")
polygon(c(1:24, 24:1), c(Delta2.q[1,], rev(Delta2.q[2,])),
        col = adjustcolor("deepskyblue", 0.3), border = NA)
axis(1, at = seq(0, 24, by = 4), labels = seq(0, 24, by = 4))
```

![](Penalised_splines_files/figure-html/results%20qreml-1.png)

``` r


par(oldpar)
```

## Smooth density estimation

To demonstrate nonparametric estimation of the state-dependent
densities, we will consider the `nessi` data set. It contains
acceleration data of the Loch Ness Monster, specifically the **overall
dynamic body acceleration (ODBA)**. ODBA is strictly positive with some
very extreme values, making direct analysis difficult. Hence, for our
analysis we consider the logarithm of ODBA as our observed process.

``` r

head(nessi)
#>         ODBA   logODBA state
#> 1 0.03775025 -3.276763     2
#> 2 0.05417830 -2.915475     2
#> 3 0.03625247 -3.317248     2
#> 4 0.01310802 -4.334531     1
#> 5 0.05402441 -2.918319     3
#> 6 0.06133794 -2.791357     3
hist(nessi$logODBA, prob = TRUE, breaks = 50, bor = "white", 
     main = "", xlab = "log(ODBA)")
```

![](Penalised_splines_files/figure-html/shark_data-1.png)

Clearly, there are at least three behavioural states in the data, and we
start by fitting a simple 3-state Gaussian HMM with likelihood function:

``` r

nll = function(par){
    getAll(par, dat)
    sigma = exp(logsigma) # exp because strictly positive
    REPORT(mu); REPORT(sigma)
    Gamma = tpm(eta) # multinomial logit link
    delta = stationary(Gamma) # stationary dist of the homogeneous Markov chain
    allprobs = matrix(1, length(logODBA), N)
    ind = which(!is.na(logODBA))
    for(j in 1:N) allprobs[ind,j] = dnorm(logODBA[ind], mu[j], sigma[j])
    -forward(delta, Gamma, allprobs)
}
```

We then fit the model as explained in the vignette *LaMa and RTMB*.

``` r

# initial parameter list
par = list(mu = c(-4.5, -3.5, -2.5),
           logsigma = log(rep(0.5, 3)),
           eta = rep(-2, 6))

# data and hyperparameters
dat = list(logODBA = nessi$logODBA, N = 3)

# creating automatically differentiable objective function
obj = MakeADFun(nll, par, silent = TRUE)

# fitting the model
opt = nlminb(obj$par, obj$fn, obj$gr)

# reporting to get calculated quantities
mod = lama(obj)

# visualising the results
color = LaMaColors(3)
hist(nessi$logODBA, prob = TRUE, breaks = 50, bor = "white",
     main = "", xlab = "log(ODBA)")
for(j in 1:3) curve(mod$delta[j] * dnorm(x, mod$mu[j], mod$sigma[j]), 
                    add = TRUE, col = color[j], lwd = 2, n = 500)
```

![](Penalised_splines_files/figure-html/Gaussian%20fit-1.png)

We see a clear lack-of-fit due to the inflexibility of the Gaussian
state-dependent densities. Thus, we now fit a model with state-dependent
densities based on P-Splines.

In a first step, this requires us to prepare the **design** and
**penalty matrices** needed using
[`smoothdens()`](https://janolefi.github.io/LaMa/reference/smoothdens.md).
This function takes the data and a set of initial parameters (specifying
initial means and standard deviations for each state). It then builds
the **P-Spline** design and penalty matrices as well as a matrix of
**initial spline coefficients** based on the provided parameters. The
basis functions are standardised such that they integrate to one, which
is needed for density estimation. For several data streams,
[`smoothdens()`](https://janolefi.github.io/LaMa/reference/smoothdens.md)
also takes a named data frame or list and returns one such object per
data stream.

``` r

# providing initial means and sds to initialise spline coefficients
par0 = list(mean = c(-4, -3.3, -2.8), sd = c(0.3, 0.2, 0.5))

# construct the smooth density object
dens = smoothdens(nessi$logODBA, par = par0)

# objects for model fitting
S = dens$S # penalty matrix for logODBA
beta = dens$coef # initial spline coefficients
```

Then, we can specify the penalised negative log-likelihood function. The
coefficient matrix `beta` provided by
[`smoothdens()`](https://janolefi.github.io/LaMa/reference/smoothdens.md)
has one column less than the number of basis functions. This is because
the last column, i.e. the last coefficient for each state, needs to be
fixed to zero for **identifiability**. `dens$density(beta)` takes care
of this: it adds the zero column and transforms the unconstrained
parameter matrix to non-negative weights that sum to one (called
`alpha`) for each state using the inverse multinomial logistic link
(softmax). The columns of the `allprobs` matrix are then computed as
linear combinations of the columns of the design matrix `dens$Z` and the
weights `alpha`, with rows of missing observations equal to one. Lastly,
we penalise the unconstrained coefficients `beta` (not the constrained
`alpha`’s) using the
[`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md)
function.

``` r

pnll = function(par){
  getAll(par, dat)
  # regular stationary HMM stuff
  Gamma = tpm(eta)
  delta = stationary(Gamma)
  # smooth state-dependent densities
  allprobs = dens$density(beta)
  # forward algorithm + penalty
  -forward(delta, Gamma, allprobs) + 
    penalty(beta, S, sp)
}
```

Now we specify the initial parameter and data list and fit the model. In
this case, we actually don’t need to add the observations to our `dat`
list, as all the information is contained in `dens`.

``` r

par = list(beta = beta, # spline coefficients prepared by smoothdens()
           eta = rep(-2, 6)) # initial transition matrix on logit scale

dat = list(S = S, # spline penalty matrix
           sp = rep(10, 3)) # initial penalty strength vector

# fitting the model using qREML
system.time(
  mod2 <- qreml(pnll, par, dat, random = "beta")
)
#> Creating AD function
#> Initialising with sp: 10 10 10
#> outer 1 - sp: 1.733 1.189 2.623 
#> outer 2 - sp: 0.836 0.885 1.46 
#> outer 3 - sp: 0.956 0.994 1.653 
#> outer 4 - sp: 0.962 0.996 1.661 
#> outer 5 - sp: 0.962 0.996 1.661 
#> outer 6 - sp: 0.962 0.996 1.661 
#> outer 7 - sp: 0.962 0.996 1.661
#> Converged
#> Final model fit with sp: 0.962 0.996 1.661
#>    user  system elapsed 
#>   9.058   2.057   8.634
```

After fitting the model, we can easily visualise the smooth densities:
[`predict()`](https://rdrr.io/r/stats/predict.html) evaluates them for
the estimated coefficients on a grid and returns the grid `x` together
with the densities. We already have access to all reported quantities
because [`qreml()`](https://janolefi.github.io/LaMa/reference/qreml.md)
automatically runs the reporting after model fitting.

``` r

pred = predict(dens, coef = mod2$par$beta)
sDens = pred$density # all three state-dependent densities on a grid
xseq = pred$x

hist(nessi$logODBA, prob = TRUE, breaks = 50, bor = "white", main = "", xlab = "log(ODBA)")
for(j in 1:3) lines(xseq, mod2$delta[j] * sDens[,j], col = color[j], lwd = 2)
lines(xseq, colSums(mod2$delta * t(sDens)), col = "black", lwd = 2, lty = 2)
```

![](Penalised_splines_files/figure-html/shark_smooth_results-1.png)

The P-Spline model results in a very good fit to the empirical
distribution. This is because the first state has a skewed distribution,
the second state has a high kurtosis and the third state has a funny
right tail. The P-Spline model can capture all of these features where
the parametric model failed to do so.

## Markov-switching GAMLSS

Lastly, we want to demonstrate how one can easily fit Markov-switching
regression models where the state-dependent means and potentially other
parameters depend on covariates via smooth functions. For this, we
consider the `energy` data set contained in the R package `MSwM`. It
comprises 1784 daily observations of energy prices (in Cents per kWh) in
Spain which we want to explain using the daily oil prices (in Euros per
barrel) also provided in the data. Specifically, we consider a 2-state
MS-GAMLSS defined by \text{price}\_t \mid \\ S_t = i \\ \sim N
\bigl(\mu_t^{(i)}, (\sigma_t^{(i)})^2 \bigr), \mu_t^{(i)} =
\beta\_{0,\mu}^{(i)} + s\_{\mu}^{(i)}(\text{oil}\_t), \quad
\text{log}(\sigma_t^{(i)}) = \beta\_{0, \sigma}^{(i)} +
s\_{\sigma}^{(i)}(\text{oil}\_t), \quad i = 1,2, not covering other
potential explanatory covariates for the sake of simplicity.

``` r

data(energy, package = "MSwM")
head(energy)
#>      Price      Oil      Gas     Coal   EurDol Ibex35   Demand
#> 1 3.188083 22.43277 14.40099 38.35157 1.134687 8.3976 477.3856
#> 2 4.953667 22.27263 19.02747 38.35157 1.106439 8.3771 609.1261
#> 3 4.730917 22.65383 18.48417 38.35157 1.106684 8.5547 650.3715
#> 4 4.531000 23.67657 18.30143 38.35157 1.116819 8.4631 647.0499
#> 5 5.141875 23.67209 14.55602 38.35157 1.122965 8.1773 627.9698
#> 6 6.322083 23.60534 15.22485 38.35157 1.122460 8.1866 693.2467
```

Similar to the first example, we can prepare the model matrices using
[`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md):

``` r

modmat = modmat(~ s(Oil, k = 12, bs = "ps"), energy)
Z = modmat$Z # design matrix
S = modmat$S # penalty matrix (list)
```

Then, we specify the penalised negative log-likelihood function. It
differs from the first example as the state-dependent distributions, as
opposed to the state process parameters, depend on the covariate.
Additionally, we now have two separate spline-coefficient matrices
called `betaSpline` and `alphaSpline` for the state-dependent means and
standard deviations respectively. Thus, we need to pass them as a list
to the
[`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md)
function.

We also pass the penalty matrix list `S` that is provided by
[`modmat()`](https://janolefi.github.io/LaMa/reference/modmat.md). This
could potentially be a list of length two if the two spline coefficient
matrices were penalised differently (e.g. by us using a different spline
basis). In this case, however, they are the same and we only pass the
list of length one. It does not matter to
[`penalty()`](https://janolefi.github.io/LaMa/reference/penalty.md) if
we pass a list of length one or just one matrix.

``` r

pnll = function(par) {
  getAll(par, dat)
  Gamma = tpm(eta) # computing the tpm
  delta = stationary(Gamma) # stationary distribution

  # regression parameters for mean and sd
  beta = cbind(beta0, betaSpline); REPORT(beta) # mean parameter matrix
  alpha = cbind(alpha0, alphaSpline); REPORT(alpha) # sd parameter matrix

  # calculating all covariate-dependent means and sds
  Mu = Z %*% t(beta) # mean
  Sigma = exp(Z %*% t(alpha)) # sd

  allprobs = cbind(dnorm(price, Mu[,1], Sigma[,1]),
                   dnorm(price, Mu[,2], Sigma[,2])) # state-dependent densities
  
  - forward(delta, Gamma, allprobs) +
    penalty(list(betaSpline, alphaSpline), S, sp)
}
```

From this point on, the model fit is now basically identical to the
previous two examples. We specify initial parameters and include an
intial penalty strength parameter in the `dat` list.

``` r

# initial parameter list
par = list(eta = rep(-4, 2), # state process intercepts
           beta0 = c(2, 5), # state-dependent mean intercepts
           betaSpline = matrix(0, nrow = 2, ncol = ncol(Z)-1), # mean spline coef
           alpha0 = c(0, 0), # state-dependent sd intercepts
           alphaSpline = matrix(0, nrow = 2, ncol = ncol(Z)-1)) # sd spline coef

# data, model matrices and initial penalty strength
dat = list(price = energy$Price, 
           Z = Z, 
           S = S, 
           sp = rep_sp(modmat, 4)) # 2 rows each of betaSpline and alphaSpline

# model fit
system.time(
  mod3 <- qreml(pnll, par, dat, random = c("betaSpline", "alphaSpline"))
)
#> Creating AD function
#> Initialising with sp: 15317.756 15317.756 15317.756 15317.756
#> outer 1 - sp: 3024.686 1531.776 2056.01 1531.776 
#> outer 2 - sp: 302.469 153.178 205.601 153.178 
#> outer 3 - sp: 30.247 15.318 20.56 15.318 
#> outer 4 - sp: 21.818 4.441 8.246 2.187 
#> outer 5 - sp: 22.318 10.491 8.229 6.273 
#> outer 6 - sp: 22.947 5.581 8.315 3.311 
#> outer 7 - sp: 22.352 8.746 8.224 4.794 
#> outer 8 - sp: 22.747 6.295 8.29 3.844 
#> outer 9 - sp: 22.457 7.975 8.244 4.375 
#> outer 10 - sp: 22.651 6.711 8.274 4.053 
#> outer 11 - sp: 22.512 7.603 8.255 4.239 
#> outer 12 - sp: 22.609 6.943 8.267 4.128 
#> outer 13 - sp: 22.528 7.415 8.257 4.193 
#> outer 14 - sp: 22.596 7.068 8.266 4.154 
#> outer 15 - sp: 22.595 7.072 8.266 4.154 
#> outer 16 - sp: 22.594 7.08 8.265 4.155 
#> outer 17 - sp: 22.591 7.094 8.265 4.156 
#> outer 18 - sp: 22.591 7.094 8.265 4.156
#> Converged
#> Final model fit with sp: 22.596 7.068 8.266 4.154
#>    user  system elapsed 
#>  12.251   0.015  12.267
```

Having fitted the model, we can visualise the results. We first decode
the most probable state sequence and then plot the estimated
state-dependent densities as a function of the oil price, as well as the
decoded time series. For the former, we again create a fine grid of oil
price values and use [`predict()`](https://rdrr.io/r/stats/predict.html)
to build the associated interpolating design matrix.

``` r

xseq = seq(min(energy$Oil), max(energy$Oil), length = 200) # sequence for prediction
Z_p = predict(modmat, newdata = data.frame(Oil = xseq)) # prediction design matrix

energy$states = viterbi(mod = mod3) # decoding most probable state sequence

Mu_plot = Z_p %*% t(mod3$beta)
Sigma_plot = exp(Z_p %*% t(mod3$alpha))

library(scales) # to make colors semi-transparent

oldpar = par(mfrow = c(1,2))

# state-dependent distribution as a function of oil price
plot(energy$Oil, energy$Price, pch = 20, bty = "n", col = alpha(color[energy$states], 0.1),
     xlab = "oil price", ylab = "energy price")
for(j in 1:2) lines(xseq, Mu_plot[,j], col = color[j], lwd = 3) # means
qseq = qnorm(seq(0.5, 0.95, length = 4)) # sequence of quantiles
for(i in qseq){ for(j in 1:2){
  lines(xseq, Mu_plot[,j] + i * Sigma_plot[,j], col = alpha(color[j], 0.7), lty = 2)
  lines(xseq, Mu_plot[,j] - i * Sigma_plot[,j], col = alpha(color[j], 0.7), lty = 2)
}}
legend("topright", bty = "n", legend = paste("state", 1:2), col = color, lwd = 3)

# decoded time series
plot(NA, xlim = c(0, nrow(energy)), ylim = c(1,10), bty = "n",
     xlab = "time", ylab = "energy price")
segments(x0 = 1:(nrow(energy)-1), x1 = 2:nrow(energy),
         y0 = energy$Price[-nrow(energy)], y1 = energy$Price[-1],
         col = color[energy$states[-1]], lwd = 0.5)
```

![](Penalised_splines_files/figure-html/energy_results-1.png)

``` r

par(oldpar)
```

## References

Koslik, Jan-Ole. 2024. *Efficient Smoothness Selection for Nonparametric
Markov-Switching Models via Quasi Restricted Maximum Likelihood*.
<https://arxiv.org/abs/2411.11498>.

Wood, Simon. 2017. *Generalized Additive Models: An Introduction with
r*. Chapman; hall/CRC.
