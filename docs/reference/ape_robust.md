# Robustness sweep over contextual assumptions (worst-case power and n_max)

Power computed from a single guess at the contextual inputs – baseline
rate, nuisance signal, covariate dependence – is fragile (Hancock &
Feng, 2025). `ape_robust()` re-runs the power analysis over a scenario
grid of those inputs and reports per-scenario power, the worst case, and
(optionally) `n_max`: the sample size that reaches the target power in
the worst scenario – the insurance-premium n.

## Usage

``` r
ape_robust(
  dgp,
  n,
  claim = c("minimum", "detect", "equivalence"),
  sesoi = NULL,
  conf = 0.95,
  nsim = 800,
  seed = NULL,
  vary,
  pin = c("ape", "coefficients"),
  grid_points = 3,
  nmax = TRUE,
  nmax_power = 0.9
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

- vary:

  Named list of scenario values for contextual inputs. Parametric route:
  `baseline`, `signal`, `correlation` (scalar); empirical route:
  `baseline`, `signal`. A length-2 element is expanded to `grid_points`
  equally spaced values; longer vectors are used as-is. Scenarios are
  the full factorial grid.

- pin:

  `"ape"` (re-invert per scenario) or `"coefficients"` (hold
  coefficients, let the implied effect drift).

- grid_points:

  Grid size used to expand length-2 `vary` elements.

- nmax:

  Run the n_max search in the worst scenario (default TRUE).

- nmax_power:

  Target power for the n_max search (default 0.90).

## Value

A `powerape_robust` object: `scenarios` (inputs, implied effect, power,
MCSE), the worst scenario, marginal mean power per input, and the n_max
result.

## Details

Two pinning modes (DESIGN.md section 3.3):

- `pin = "ape"` (default): the inversion is re-run in every scenario, so
  the true effect stays at the planning value and only the *precision*
  channel varies. Scenarios where the target is infeasible are reported
  as such, not dropped silently.

- `pin = "coefficients"`: the effect coefficients stay at their
  base-case values while the contextual calibration is redone; the
  implied true effect drifts and is reported per scenario.
  Claim-boundary guards are relaxed here on purpose – showing that a
  scenario collapses power is the point.

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
ape_robust(d, n = 700, claim = "detect",
           vary = list(baseline = c(0.20, 0.40)),
           nsim = 300, seed = 1, nmax = FALSE)
#> powerape robustness sweep -- detect claim, APE, pin = ape
#>   n = 700, nsim = 300 per scenario, 3 scenario(s) over: baseline
#>   power range [0.773, 0.890]; worst scenario:
#>  baseline implied_effect     power       mcse
#>       0.4            0.1 0.7733333 0.02417222
#>   scenarios:
#>  baseline implied_effect     power       mcse
#>       0.2            0.1 0.8900000 0.01806470
#>       0.3            0.1 0.8466667 0.02080242
#>       0.4            0.1 0.7733333 0.02417222
# }
```
