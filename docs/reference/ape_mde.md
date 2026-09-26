# Minimum detectable APE/AIE at a given sample size

The inverse of
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md):
the sample size is fixed – an archive with so many firms, a panel with
so many units and waves – and the question is the smallest effect the
design can reliably conclude. For `claim = "detect"` and
`claim = "minimum"` the function searches over the assumed true effect,
re-running the APE/AIE inversion at every candidate, until simulated
power at `n` hits the goal; for `claim = "equivalence"` the target is
held at the DGP's pinned value (typically 0) and the search is over the
margin instead, returning the tightest equivalence bounds the design can
expect to establish.

## Usage

``` r
ape_mde(
  dgp,
  n,
  power = 0.8,
  claim = c("detect", "minimum", "equivalence"),
  sesoi = NULL,
  alpha = 0.05,
  conf = NULL,
  nsim = 1000,
  seed = NULL,
  main_focal = NULL,
  main_moderator = NULL,
  max_iter = 6,
  confirm = TRUE,
  nsim_confirm = 4 * nsim,
  se = c("model", "robust"),
  direction = NULL,
  separation = c("fail", "keep")
)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
  or
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md).

- n:

  Sample size of the planned study (units for panel designs).

- power:

  Target power for the claim (default 0.80).

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
  candidate `n`, so
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  designs jointly over the pair (alpha(n), n); results store the rule
  and the level it produced.

- conf:

  Optional override: the interval level used for the claim itself
  (`1 - alpha` for detection, `1 - 2 alpha` otherwise). Supplying
  `conf = 0.95` for a minimum-effect or equivalence claim reproduces the
  pre-1.8.0 behavior (one-sided error rate 2.5%).

- nsim:

  Number of simulation replications.

- seed:

  Optional seed (the caller's RNG state is preserved).

- main_focal, main_moderator:

  Main-effect anchors for AIE designs; defaults to the values stored by
  a previous
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
  call.

- max_iter:

  Maximum search refinements.

- confirm, nsim_confirm:

  Confirmation stage as in
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md).

- se:

  Standard errors for the exogenous cross-sectional routes: `"model"`
  (default, expected-information ML) or `"robust"`
  (heteroskedasticity-robust HC0 sandwich, as in the sandwich package).
  Panel designs always use unit-clustered SEs and IV designs the stacked
  method-of-moments robust sandwich; `se` is ignored there.

- direction:

  `"positive"` (an increase) or `"negative"` (a decrease) for the detect
  and minimum-effect claims; default: the sign of the DGP's pinned
  effect, `"positive"` when none is pinned. Ignored for equivalence.

- separation:

  How simulated studies with separation in an identifying cell are
  scored: `"fail"` (default; counted as failed, as the field's reference
  analysis refuses them) or `"keep"` (retained with their degenerate
  Wald intervals). See the section 'Separation and sparse cells'.

## Value

A `powerape_mde` object: `mde` (the minimum detectable effect, signed,
or for equivalence the smallest establishable margin), the confirmed
`power` and `mcse` at that effect, the search `history`, the
`direction`, the separation diagnostics of the final run (see
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)),
and the DGP re-pinned at the answer (so the object feeds
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)).

## Details

**Direction.** For the detect and minimum-effect claims the search runs
in one direction: `direction = "positive"` searches increases,
`"negative"` decreases (a drop in a rare event, an attenuating
interaction). With a binary outcome the two differ whenever the baseline
is not 0.5, because the Bernoulli variance moves with the rate: above a
0.5 baseline a decrease is harder to detect than an increase of the same
size, below it easier. The default follows the sign of the DGP's pinned
effect (`"positive"` for an unpinned DGP), and the returned `mde`
carries the sign.

**Search.** A pilot simulation measures the standard error at a
reference effect (the pinned one, or half the distance from the baseline
to the nearer bound), a normal approximation proposes a candidate, and
the standard error is re-measured at that candidate before the search
starts: the standard error of an APE moves with the effect through the
Bernoulli variances, so a standard error taken far from the answer would
bias the first candidate (by 2-4% at large `n` from a 0.15 reference at
a 0.30 baseline, more at rarer baselines). By default the answer is then
verified the way
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
verifies its n: a high-precision confirmation stage re-measures power at
the candidate and pushes the effect outward (never inward) if it falls
short, so the reported MDE errs on the conservative side.

For a DGP with a moderator the searched effect is the AIE; the two
conditional-at-reference main-effect anchors are taken from the pinned
DGP (or passed via `main_focal`/`main_moderator`) and held fixed across
candidates.

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
ape_mde(d, n = 712, claim = "detect", nsim = 600, seed = 1)
#> powerape minimum detectable APE -- detect claim
#>   MDE = 0.0996 at n = 712 for 80% target power (confirmed 0.804, MCSE 0.008)
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 2400.
## a decrease from the same baseline
ape_mde(d, n = 712, claim = "detect", nsim = 600, seed = 1,
        direction = "negative")
#> powerape minimum detectable APE -- detect claim, decreases
#>   MDE = -0.0908 at n = 712 for 80% target power (confirmed 0.801, MCSE 0.008)
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 2400.
# }
```
