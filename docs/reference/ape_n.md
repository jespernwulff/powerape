# Required sample size for a target power

Search plus confirmation. A pilot run estimates the SE scale (SE ~
c/sqrt(n)), a normal-approximation formula proposes n, and simulation at
the proposal verifies and adjusts until the estimated power is within
Monte Carlo precision of the goal. By default a **confirmation stage**
then re-measures the candidate at `nsim_confirm` replications and
accepts only if the confirmed power is no more than
`max(0.005, 1.5 * MCSE)` below the goal, pushing n upward otherwise (up
to three rounds). The confirmation never trims n on overshoot – slight
overshoot is the conservative direction for sample-size planning. The
returned power and MCSE come from the confirmation run. Works for APE
and AIE designs alike.

## Usage

``` r
ape_n(
  dgp,
  power = 0.9,
  claim = c("minimum", "detect", "equivalence"),
  sesoi = NULL,
  alpha = 0.05,
  conf = NULL,
  nsim = 1500,
  seed = NULL,
  n_range = c(30, 2e+06),
  max_iter = 5,
  confirm = TRUE,
  nsim_confirm = 4 * nsim,
  se = c("model", "robust"),
  separation = c("fail", "keep")
)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
  or
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md).

- power:

  Target power (default 0.90).

- claim:

  `"minimum"` (default), `"detect"`, or `"equivalence"`.

- sesoi:

  Smallest effect size of interest, in APE units. Required for
  `"minimum"` and `"equivalence"`; optional for `"detect"` (if supplied,
  the outcome table is still broken out against it).

- alpha:

  The claim's error rate (default 0.05): two-sided for `"detect"`,
  one-sided for `"minimum"`, TOST for `"equivalence"`. May also be a
  **sample-size rule**: a function of `n` returning the error rate to
  use at that `n`, for example `function(n) alphaN::alphaN(n, BF = 3)`,
  the Bayes-factor calibration of Wulff and Taylor (2024) that lowers
  alpha as `n` grows. The searching functions evaluate the rule at every
  candidate `n`, so `ape_n()` designs jointly over the pair (alpha(n),
  n); results store the rule and the level it produced.

- conf:

  Optional override: the interval level used for the claim itself
  (`1 - alpha` for detection, `1 - 2 alpha` otherwise). Supplying
  `conf = 0.95` for a minimum-effect or equivalence claim reproduces the
  pre-1.8.0 behavior (one-sided error rate 2.5%).

- nsim:

  Replications per search step.

- seed:

  Optional seed (the caller's RNG state is preserved).

- n_range:

  Admissible range for the answer; hitting the ceiling warns.

- max_iter:

  Maximum search iterations.

- confirm:

  Run the high-precision confirmation stage (default TRUE). Set to FALSE
  for speed in exploratory loops; the returned power then carries
  search-stage precision only.

- nsim_confirm:

  Replications per confirmation round (default `4 * nsim`).

- se:

  Standard errors for the exogenous cross-sectional routes: `"model"`
  (default, expected-information ML) or `"robust"`
  (heteroskedasticity-robust HC0 sandwich, as in the sandwich package).
  Panel designs always use unit-clustered SEs and IV designs the stacked
  method-of-moments robust sandwich; `se` is ignored there.

- separation:

  How simulated studies with separation in an identifying cell are
  scored: `"fail"` (default; counted as failed, as the field's reference
  analysis refuses them) or `"keep"` (retained with their degenerate
  Wald intervals). See the section 'Separation and sparse cells'.

## Value

A `powerape_n` object: the required `n`, the confirmed (or, with
`confirm = FALSE`, search-stage) power and MCSE at that `n`, the search
history (`stage` column distinguishes search and confirmation rounds),
the standard-error type and separation convention used, the separation
diagnostics of the final run (`separated`, `min_cell`; see
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)),
and the embedded DGP.

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
ape_n(d, power = 0.8, claim = "detect", nsim = 600, seed = 1)
#> powerape required sample size -- detect claim
#>   n = 706 for 80% target power (confirmed 0.805, MCSE 0.008)
#>   assumed true APE +0.1000, sesoi -, 95% CI (alpha = 0.05), probit
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 2400.
# }
```
