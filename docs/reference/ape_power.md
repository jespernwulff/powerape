# Simulated power for an APE/AIE claim at a given sample size

Simulates studies from the pinned DGP, fits the index model by maximum
likelihood, computes the APE (or, for moderated DGPs, the AIE) with a
delta-method Wald CI, and evaluates the requested confidence-interval
claim (Riesthuis, 2024): `"minimum"` (CI lower bound above the SESOI –
the default and the package's point), `"detect"` (CI excludes 0,
directional), or `"equivalence"` (CI within +/- SESOI). Also reports the
full outcome distribution.

## Usage

``` r
ape_power(
  dgp,
  n,
  claim = c("minimum", "detect", "equivalence"),
  sesoi = NULL,
  alpha = 0.05,
  conf = NULL,
  nsim = 1000,
  seed = NULL,
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

  Total sample size of the simulated study.

- claim:

  `"minimum"` (default), `"detect"`, or `"equivalence"`.

- sesoi:

  Smallest effect size of interest, in APE units. Required for
  `"minimum"` and `"equivalence"`; optional for `"detect"` (if supplied,
  the outcome table is still broken out against it).

- alpha:

  The claim's error rate (default 0.05): two-sided for `"detect"`,
  one-sided for `"minimum"`, TOST for `"equivalence"`.

- conf:

  Optional override: the interval level used for the claim itself
  (`1 - alpha` for detection, `1 - 2 alpha` otherwise). Supplying
  `conf = 0.95` for a minimum-effect or equivalence claim reproduces the
  pre-1.8.0 behavior (one-sided error rate 2.5%).

- nsim:

  Number of simulation replications.

- seed:

  Optional seed (the caller's RNG state is preserved).

- se:

  Standard errors for the exogenous cross-sectional routes: `"model"`
  (default, expected-information ML) or `"robust"`
  (heteroskedasticity-robust HC0 sandwich, as in the sandwich package).
  Panel designs always use unit-clustered SEs and IV designs the stacked
  method-of-moments robust sandwich; `se` is ignored there.

## Value

A `powerape_power` object: power, Monte Carlo standard error, outcome
distribution, the failed-fit count (failures count against power,
conservatively), and the embedded DGP spec for reproducibility and
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md).

## Details

**Error-rate convention.** Every claim is tested at `alpha` (default
0.05) in its conventional form. Detection is the two-sided test at
`alpha`, read off the `1 - alpha` (95%) interval (only correctly signed
rejections count, which costs no power). The minimum-effect claim is a
one-sided test at `alpha` and equivalence a two-one-sided-tests (TOST)
procedure at `alpha`; both read off the `1 - 2 alpha` (90%) interval,
the convention of Lakens (2017), Lakens et al. (2018), TOSTER, and
Riesthuis (2024) for equivalence. The outcome distribution uses both
intervals accordingly. To reproduce the more conservative choice of a
95% interval for the minimum-effect test (Riesthuis, 2024) pass
`conf = 0.95`, which sets that claim's one-sided error rate to 2.5%.

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
ape_power(d, n = 800, claim = "minimum", sesoi = 0.03, nsim = 500, seed = 1)
#> powerape -- minimum-effect claim (CI lower bound > 0.030)
#>   probit, n = 800, assumed true APE +0.1000, 90% CI (alpha = 0.05), nsim = 500
#>   power = 0.678 (MCSE 0.021)
#>   outcomes: minimum 0.678 | detect-only 0.176 | inconclusive 0.146 | equivalence 0.000 | failed 0.000 
# }
```
