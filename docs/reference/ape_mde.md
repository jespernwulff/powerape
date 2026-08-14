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
  conf = 0.95,
  nsim = 1000,
  seed = NULL,
  main_focal = NULL,
  main_moderator = NULL,
  max_iter = 6,
  confirm = TRUE,
  nsim_confirm = 4 * nsim,
  se = c("model", "robust")
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

- conf:

  CI level (default 0.95).

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

## Value

A `powerape_mde` object: `mde` (the minimum detectable effect, or for
equivalence the smallest establishable margin), the confirmed `power`
and `mcse` at that effect, the search `history`, and the DGP re-pinned
at the answer (so the object feeds
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)).

## Details

By default the answer is verified the way
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
verifies its n: a high-precision confirmation stage re-measures power at
the candidate and pushes the effect upward (never downward) if it falls
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
#>   MDE = 0.1003 at n = 712 for 80% target power (confirmed 0.808, MCSE 0.008)
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 2400.
# }
```
