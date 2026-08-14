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
  conf = 0.95,
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

- conf:

  CI level (default 0.95).

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

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
ape_power(d, n = 800, claim = "minimum", sesoi = 0.03, nsim = 500, seed = 1)
#> powerape -- minimum-effect claim (CI lower bound > 0.030)
#>   probit, n = 800, assumed true APE +0.1000, 95% CI, nsim = 500
#>   power = 0.554 (MCSE 0.022)
#>   outcomes: minimum 0.554 | detect-only 0.300 | inconclusive 0.146 | equivalence 0.000 | failed 0.000 
# }
```
