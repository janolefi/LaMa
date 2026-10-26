# qreml(): what to improve, and why

Notes taken while reading `R/qreml_functions.R` and `R/helper_functions.R`
alongside an independent implementation of the same algorithm (the `efs`
engine in gamRTMB, `R/efs.R`). Nothing here has been changed; this is a list
for a later session, ordered roughly by how much it matters.

Two things to keep in mind while reading. First, `qreml` solves a strictly
harder problem than a GAM fitter does: the likelihood is the user's, written
by hand, and the penalty structure is declared by convention rather than
constructed. Several of the awkward bits below are the price of that, and the
fix is a clearer convention rather than less flexibility. Second, the
algorithm itself is right — the two implementations agree on the smoothing
parameters when both converge. What follows is about correctness at the
edges, cost, and how easy the thing is to reason about.

Notation here follows the code: `J` is the penalized Hessian (the Hessian of
`pnll`, which already contains the penalty), and `H = J - bigS` is the data
Hessian. Note that this is the opposite of the usual convention, where `H` is
the penalized one; the comment at the assignment says "penalised Hessian"
about the unpenalized matrix, which is worth fixing on its own.

---

## 1. The tensor-product update is off by a factor of lambda

This is a real bug, not a style point.

`penalty()` and the simple branch of `penalty2()` report the bare quadratic
form:

```r
quadform <- rowSums(current_re %*% S[[i]] * current_re)
Pen[[i]] <- quadform                       # b' S b
```

but the tensor-product branch of `penalty2()` reports the *weighted* one:

```r
subthispen[k] <- thislambda[k] * (t(this_tp[j,]) %*% thisS[[k]] %*% this_tp[j,])
thispen[[j]]  <- subthispen                # lambda_k b' S_k b
```

The update in `qreml` is `lambda_new <- this_edoF / this_pen` in both cases.
For a simple smooth that is exactly Wood & Fasiolo:

```
edoF = rank - lambda tr(J^-1 S) = lambda [tr(S_lambda^- S) - tr(J^-1 S)]
pen  = b' S b
=> lambda_new = lambda * [tr(S_lambda^- S) - tr(J^-1 S)] / (b' S b)     correct
```

For a tensor product both numerator and denominator carry a `lambda_k`, and it
cancels:

```
edoF_k = lambda_k [tr(S_lambda^- S_k) - tr(J^-1 S_k)]
pen_k  = lambda_k b' S_k b
=> lambda_new = [tr(...) - tr(...)] / (b' S_k b) = lambda_new_correct / lambda_k
```

So every tensor-product smooth takes a step that is the correct one divided by
its own current `lambda`. Near `lambda = 1` this is invisible; away from it,
it converges to a different fixed point. The outer gradient reported by
`outer_gr` has the same asymmetry, so the gradient convergence criterion does
not catch it.

**Fix**: make `penalty2()` report `b' S_k b` in both branches, and let the
caller apply `lambda`. Two callers change; both are in this file.

A regression test worth having: a one-smooth model where `S` is passed once as
a matrix (simple branch) and once as a one-element list (tp branch). The two
should give the same `lambda`.

## 2. Tensor products and `L` matrices need a different step, not a different loop

Related to the above, and the deeper version of it. The multiplicative update
is derived for penalty weights that can move one at a time. For a tensor
product they can — `lambda_1 S_1 + lambda_2 S_2` with both free — so the fix
in §1 is enough there.

It is *not* enough for a smooth whose weights are tied through mgcv's `L`
convention, `lambda = exp(L theta)`, which is how an SPDE or an adaptive
smooth is written. There the weight-by-weight update, projected back onto
`theta`, has a fixed point where the projected log-ratios vanish, which is not
where the gradient vanishes. Measured in the gamRTMB implementation: it
converged smoothly and confidently to a point 1.3 nats short of the optimum
with `|dV/dtheta| = 3.9`, and every diagnostic said "converged".

`qreml` does not support `L` today, so this is only a note for if it ever
does. What works is a step in the metric the multiplicative update induces:
linearising `log(n/d)` about `n = d` gives `2(n-d)/(n+d)`, so the exact update
is to first order `-M^{-1} g` with `M = L' diag((n_i + d_i)/4) L`. `M` is
positive definite, so the step is always downhill, it reproduces the
multiplicative update on the free-weight case, and its fixed point is `g = 0`.

## 3. No step-size control

`alpha` applies fixed exponential smoothing to every proposed `lambda`:

```r
lambda_mapped[i] <- (1 - alpha) * lambda_new + alpha * lambda_mapped[i]
```

unconditionally, whether the step improves the criterion or not. Wood &
Fasiolo's safeguard is different in kind: halve the step in `log lambda` until
the criterion does not increase. That is what makes the method monotone, and
monotonicity is what lets you trust `llk_r` as a convergence signal.

Everything needed is already computed. `llk_r[k]` is the restricted likelihood
at the current `lambda`; it is stored and then never read. Turning it into a
line search is a loop around the existing inner solve:

```
propose step; a <- 1
repeat { evaluate criterion at lambda * exp(a * step)
         if it improved, accept; else a <- a/2 }
```

with a cap on `max_halve`. The inner solve is warm-started, so a rejected
candidate is cheap after the first one.

Note also that damping with `alpha` in the *linear* `lambda` scale is not the
same as damping in `log lambda`, and the multiplicative update lives in the
log scale. `alpha = 0.3` on a `lambda` that wants to move from 1 to 1e6 is a
very different amount of damping than on one moving from 1 to 2.

## 4. The full inverse is computed where a few blocks are needed

```r
J_inv <- safe_chol_inv(J_pd)     # the entire p x p inverse
```

and then every use is `J_inv[idx, idx]` for one smooth's indices, or
`rowSums(J_inv * (J - bigS))` for the EDF. Nothing ever needs an off-block
entry of the inverse. This is the single biggest obstacle to scaling: it is
`O(p^3)` time and `O(p^2)` memory every outer iteration, on a matrix that is
often sparse.

Two cheaper routes, both using the Cholesky that `safe_chol_inv` already
forms:

- Per smooth, solve the factor against that smooth's `q_k` indicator columns.
  Peak memory `p * max(q_k)` instead of `p^2`, and the total work is one
  triangular solve per penalized coefficient.
- Or a Takahashi / inverse-subset recursion, which gives the inverse exactly
  on the factor's sparsity pattern — a superset of what any of these traces
  needs.

The first is a few lines and is what gamRTMB's engine does; the second is
better asymptotically and is more work.

Same point for `leading_diag <- rowSums(J_inv * (J - bigS))` at the end: that
is `diag(J^-1 (J - bigS))`, which only needs the diagonal, and the comment
already says so — but it is computed from a dense inverse anyway.

## 5. `MASS::ginv` where a Cholesky would do, and where the structure is known

Three places:

- The final `J_inv <- MASS::ginv(J)`, after the loop, even though
  `safe_chol_inv` exists and was used inside it. `ginv` is an SVD: `O(p^3)`
  with a large constant, and it silently pseudo-inverts a singular `J`
  rather than reporting it.
- `thisS_inv <- MASS::ginv(thisS)` per tensor-product smooth per outer
  iteration. The null space of `sum_k lambda_k S_k` does not depend on the
  `lambda_k` (for positive weights and PSD `S_k` it is the intersection of the
  null spaces), so it can be found once at setup — or removed entirely by
  corner-constraining the block, which is what makes the penalty proper.
- `mod$outer_gr` repeats both.

## 6. PD repair: right idea, three inconsistencies

The repair targets `H = J - bigS`, which is correct — Wood's Theorem 1 is a
condition on the data Hessian, and a test on `J` passes routinely while
`J - bigS` is indefinite. But:

- The sparse branch adds `1e-8 * mean(diag(H))` **once** and does not re-check;
  if that was not enough, the code proceeds with a still-indefinite matrix.
  The dense branch calls `Matrix::nearPD`. Two different repairs for the same
  condition, one of which can silently fail.
- `mean(Matrix::diag(H))` can be negative, in which case the "jitter" is
  negative and makes things worse. Use `max(abs(diag(H)))`, or a fraction of
  the largest eigenvalue.
- `llk_r[k] <- restr_llk(-opt$value, bigS, J)` uses the **unrepaired** `J`,
  while the traces in the same iteration use the repaired one. The two halves
  of the iteration are then describing different matrices.
- The final EDF are computed from `ginv(J)` on the unrepaired `J`, so a fit
  that needed repair reports effective degrees of freedom that need not lie in
  `[0, 1]`.

Worth distinguishing, when reporting this to the user, between a repair at an
intermediate iteration — routine, the starting values are often a place with
negative curvature and the iteration walks away from it — and a repair at the
final point, which means the reported criterion and EDF are those of a
different problem. Only the second deserves a warning.

## 7. The criterion is computed expensively and then not used

```r
restr_llk <- function(lp_opt, bigS, J) lp_opt + gdeterminant(bigS)/2 - gdeterminant(J)/2
```

`gdeterminant` is an **eigendecomposition**. Two of them per outer iteration,
on `p x p` matrices, for a quantity that is stored in `llk_r` and read by
nothing. Either delete it or make it earn its keep (§3). If it stays:

- `log|J|` comes free from the Cholesky that is already formed.
- `log|S_lambda|_+ = sum_k rank_k log(lambda_k) + sum_k log|S_k|_+` for
  distinct blocks, and the second sum does not depend on `lambda` at all, so
  it never needs recomputing — and for a convergence test it can be dropped
  entirely as an additive constant.

## 8. Convergence criteria

Three separate things decide convergence and none of them is the natural one:

- `conv_crit = "relchange"` (the default) on `lambda`. Scale-free, but
  insensitive exactly where `lambda` is large, which is where terms are being
  shrunk away.
- `opt$counts[2] < 3` — "the inner optimiser needed fewer than three gradient
  evaluations" — declares outer convergence. That is a statement about the
  warm start, not about the smoothing parameters, and it will fire early on an
  easy inner problem.
- `k > 3` forces at least four outer iterations even when the first one was
  already at the optimum.

The natural criterion is the outer gradient, which is already computed and
stored in `outer_gr`: `-0.5 * pen + edoF / (2 * lambda)`. With §3 in place, a
relative change in the criterion is the other good one. Both are cheap and
both mean what they say.

## 9. Two thresholds for "effectively infinite"

`lambda` is clamped at `1e8`, excluded from the convergence check above `1e6`,
and only after `k > 3`. Three magic numbers for one idea. Worse, exclusion is
permanent in effect: `convInd <- which(lambda_mapped <= 1e6)` is recomputed
each iteration, so a `lambda` can come back — but its gradient was not
examined while it was out, so nothing was driving it back.

A cleaner statement: bound `log lambda` to a range (`lambda` in `[1e-10,
1e10]`, say), let the step control handle the approach to the boundary, and
report which parameters ended up on it.

## 10. The inner optimiser is BFGS, and that is correct — do not change it

**This item said the opposite when it was first written.** It argued that
since the Hessian is formed anyway for the update, a Newton inner solve is
nearly free and would be more reliable. That reasoning is wrong, it was
implemented in gamRTMB's EFS engine on the strength of it, and the engine was
twenty times slower than `qreml` on `bcpe`/`film90` and did not converge at
all. The correction is recorded here because the argument is seductive.

The error is in "the Hessian is formed anyway". It is formed **once per outer
iteration**, at the converged inner point, because that is where the update's
traces and the criterion's log determinant are needed. A Newton inner loop
forms and factorises it once per inner *step* — on that model, a hundred per
outer iteration instead of one.

Three further things BFGS gets right that a Newton loop has to be built to
get right:

- It never forms the Hessian during the solve, so an indefinite one costs
  nothing. `bcpe`'s inner Hessian **is** indefinite, by the family's
  parameterisation rather than by the starting values, so a Newton step needs
  the Hessian repaired before it means anything, and the repair changes the
  step's length to whatever the correction happened to give.
- Its direction comes from a positive definite approximation, so the scale is
  always sensible.
- `optim`'s line search enforces a curvature condition as well as sufficient
  decrease. A hand-rolled backtracking loop that accepts any decrease at all
  can creep indefinitely — the symptom in gamRTMB was a two-cycle taking
  0.001 off the objective per iteration for a hundred iterations.

Measured, same model, same machine: `qreml` 56s, EFS with a BFGS inner solve
54s and converged, EFS with a Newton inner solve still running after twenty
minutes. So `reltol = 1e-10` and the advice not to touch it are not symptoms
of anything; leave them alone.

Where a Newton inner solve *does* win is a large structured term whose inner
problem is well conditioned, because a sparse Hessian beats BFGS's dense
`p x p` approximation: on a Markov random field over 196 regions, 0.28s
against 0.71s. Worth offering as an option, never as the default.

Related, and unaffected by the above: `spHess = FALSE` by default, so
`optimHess` finite-differences the AD gradient — `p` gradient calls and a
dense `p x p` result, once per outer iteration. The sparse AD Hessian is
already implemented and is opt-in. It should at least be the default when the
Hessian turns out to be sparse, which can be checked once at setup.

## 11. `obj$fn` is not a pure function of its argument

`lambda` reaches the tape through

```r
getLambda <- function(x) lambda
dat[[psname]] <- DataEval(getLambda, rep(advector(1), 0))
```

reading a variable in the enclosing frame that the outer loop reassigns. It
works, but the consequences are worth knowing:

- `obj$fn(p)` after the loop uses whatever `lambda` was last, so the returned
  object's `obj` is not the object for the returned `lambda` unless the final
  refit happens to have left it there.
- `mod$outer_gr` re-enters the same environment (`environment(mod$outer_gr) <-
  environment()`) and mutates `newpar` as a side effect of being called.
- `?Tape` documents `DataEval` as caching until `F$force.update()` is
  called, and nothing in `qreml` calls it — which looks alarming. **It is
  not a bug**: checked directly on the exact `qreml` construction, changing
  `lambda` and re-evaluating at a fresh parameter vector gives the right
  answer to machine precision for `lambda` = 1, 100 and 10000. An
  un-reordered tape re-reads on every pass; the caching in the help page's
  example only starts after `F$reorder()`.

  What *does* bite, and what makes this look broken when you go to check it,
  is that TMB caches `obj$fn(p)` for a repeated `p`. Evaluating at the same
  parameter vector after changing `lambda` returns the old value and looks
  exactly like a stale `DataEval`. Nudge the parameter vector when testing
  this.

An alternative worth considering is to make `lambda` an ordinary parameter of
the objective and simply hold it fixed during the inner solve, which is what
it is mathematically. Then the tape has no hidden state, and the inner
gradient and Hessian are sub-blocks of the joint ones.

While on RTMB sharp edges, one more, found the hard way: `RTMB::OBS(x)` keys
on the **deparsed name** of its argument in a registry global to the package.
`MakeADFun` resets that registry per object, but `MakeTape` does not. So a
tape built with `MakeTape` that calls `OBS(y)` leaves `y` behind under that
name, and the *next* model's `OBS(y)` — on a tape or in plain R — silently
returns the previous model's response. It presents as a second fit that
converges confidently, gradient `1e-10`, to the wrong answer. If `qreml` ever
moves to `MakeTape`, or if a user's `pnll` calls `OBS`, this is waiting.

## 12. Conventions that are checked by comment rather than by code

- The ordering of `random` must match the order in which coefficients were
  passed to `penalty()`. The docs say so in bold; a mismatch silently fits a
  different model. It is checkable: `length(Pen[[i]])` against
  `nrow(re_inds[[i]])`, and the dimensions of each `S[[i]]` against the
  corresponding parameter's.
- `dat` must be called `dat`: `argname_dat <- as.character(substitute(dat))`
  and then `assign(argname_dat, dat, envir = environment())`. This works only
  because the user's `pnll` refers to it by the same name. The code comment
  says "(hopefully)". Either document it as a hard requirement and check it,
  or pass the data explicitly.
- `argname_par <- as.character(substitute(par))` makes the returned list's
  element name depend on what the caller happened to name their variable —
  `mod$par` or `mod$mypars`. Downstream code has to read `mod$parname` to find
  it. A fixed name plus a recorded original would be simpler for everyone.
- `penalty()` guesses coefficient orientation:

  ```r
  if (ncol(re_coef[[i]]) != nrow(S[[i]])) t(re_coef[[i]])
  else if (nrow(...) != nrow(S[[i]])) re_coef[[i]]
  else stop(...)
  ```

  so a square coefficient matrix whose dimension matches the penalty is an
  error, and a non-square one is transposed without comment. Orientation
  should be declared once (an argument, or "always one smooth per row") rather
  than inferred.

## 13. Smaller things

- `build_bigS` allocates a dense `p x p` matrix every outer iteration, and
  again at the end, even when every `S` is sparse.
- `ranks <- sapply(S, function(x) ... Matrix::rankMatrix(x))` is an SVD per
  penalty at setup. Once, so it matters less, but the ranks of the standard
  difference penalties are known analytically, and `rankMatrix` on a large `S`
  is not cheap.
- `gc()` before returning. Nothing in the function warrants forcing a
  collection; it just adds a pause.
- `smoothing` is applied after convergence and then the model is refitted, so
  `llk_restricted` describes a trajectory that does not end at the returned
  fit. Worth saying in the docs.
- The "no penalty parameters will be estimated" case sets `maxiter <- 1` but
  still runs the whole update machinery for one iteration.
- `allmods` is accumulated when `saveall = TRUE` and holds a full report list
  per iteration; on a large model that is a lot of memory for a debugging
  aid.

---

## What a rewrite would look like

If more than the bug in §1 is worth fixing, the shape that came out of doing
this independently is:

1. Build the objective once, with the penalty weights as ordinary parameters
   held fixed during the inner solve.
2. Keep the BFGS inner solve (see §10), warm-started, and form the sparse AD
   Hessian once per outer iteration at the converged point.
3. Per outer iteration: form `J` and `H = J - bigS`; check `H` is PSD and
   repair it if not; factor `J` once and reuse the factor for the log
   determinant, the block inverses, and the Newton step.
4. Propose the multiplicative step for free weights and the metric step of §2
   for tied ones; cap it; halve it against the real criterion until it
   improves.
5. Converge on the outer gradient or the relative change in the criterion.

That is about 250 lines including the documentation, and it is the same
algorithm — the differences are all in §§1, 3, 4, 6 and 8.
