Cleans up `forward()`, `forward_g()` and `forward_p()` before adding per-observation scores. Best reviewed commit by commit.

## Bug fixes
- **Unsorted track IDs (AD only):** track-specific `delta`/`Gamma` were matched to the wrong track whenever the IDs were not in sorted order, e.g. `"ID1"`, …, `"ID10"`. In a realistic fit, `obj$fn` and `nll(par)` differed by 39 at the same parameters.
- Tracks of length one, and `nObs = 1`, failed with AD.
- `logspace = TRUE` added `1e-7` to `delta`, which could change the likelihood by ~1000 in edge cases. Without AD it underflowed to `NaN`. It now uses an exact row-max shift.
- One-state models returned `NaN` with `logspace`, and skipped input checks and `REPORT`.
- `forward_p()` failed if `tod` did not contain every value 1, …, L.
- `forward_g()` failed without AD for two observations.
- Interleaved (non-contiguous) tracks now give an informative error instead of `Mat::row(): index out of bounds`.

## Structure
- `prepare_forward()` does all input checks. It returns `gidx`, the `Gamma` slice used for the transition into each time point (homogeneous, track-specific, time-varying or periodic).
- There is one R kernel for AD (`forward_ad()`, exact and banded) and one C++ kernel (`forward_cpp()`), replacing four C++ functions and four copies of the R recursion.
- `forward_g()` calls `forward()`. `forward_p()` indexes `Gamma` by `tod` and only expands it for `REPORT`, because the decoding functions expect that.
- Without AD, `NA` in `allprobs` is an error. The TapeConfig tip is shown once per session.

## Tests
- New `tests/testthat/test-forward.R` compares the AD and C++ versions with an exact log-scale reference, including gradients.
- `R CMD check` (without vignettes, since pandoc is missing locally) shows one NOTE, which was already there before (`tail` in `summary.qremlModel`).
- The R code of all 9 vignettes runs without errors or new warnings.

🤖 Generated with [Claude Code](https://claude.com/claude-code)
