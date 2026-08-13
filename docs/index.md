# powerape

**Power analysis for average partial effects (APEs) and average
interaction effects (AIEs) from probit and logit models.**

In nonlinear models the coefficient is not the quantity anyone
interprets — conclusions are stated in probability points, i.e. the APE.
But the APE depends on *all* coefficients and the covariate
distribution, so coefficient- or Cohen’s-d-based power tools cannot
answer the question practitioners actually have: *“How many observations
do I need to detect a treatment effect of 5 percentage points?”*

`powerape` answers it by simulation, with targets specified **in APE
units**, built around the confidence-interval claims of Riesthuis
(2024):

| Claim | CI criterion | Conclusion it buys |
|----|----|----|
| `minimum` (default) | lower bound \> SESOI | the effect is *meaningfully* large |
| `detect` | CI excludes 0 | the effect exists (directional) |
| `equivalence` | CI inside ±SESOI | the effect is too small to matter |

## Installation

``` r

# install.packages("remotes")
remotes::install_github("jespernwulff/powerape")
```

## Quick start

``` r

library(powerape)

## 1. Describe the world: baseline rate + covariate signal, no coefficients
d <- ape_dgp(
  model       = "probit",
  focal       = pa_var("treat", "binary", p = 0.5),
  covariates  = list(pa_var("age", "normal", mean = 45, sd = 12),
                     pa_var("female", "binary", p = 0.55)),
  correlation = 0.2,
  baseline    = 0.30,
  signal      = 0.15
)

## 2. Pin the true APE (inversion + feasibility check)
d <- set_ape(d, target = 0.10)

## 3. Power the minimum-effect claim against a 5-pp SESOI
ape_power(d, n = 1500, claim = "minimum", sesoi = 0.05, nsim = 1000)

## 4. Required n, power curves, robustness, reporting
ape_n(d, power = 0.80, claim = "minimum", sesoi = 0.05)
plot(ape_curve(d, n = seq(500, 3000, 500), claim = "minimum", sesoi = 0.05),
     target_power = 0.80)
ape_robust(d, n = 2200, claim = "minimum", sesoi = 0.05,
           vary = list(baseline = c(0.20, 0.40), signal = c(0, 0.30)))
power_statement(ape_power(d, n = 2200, claim = "minimum", sesoi = 0.05))
```

## Features

- **Binary and continuous focal variables** (discrete-change and
  average-derivative APEs), with plateau-aware inversion and informative
  feasibility errors.
- **Interaction designs**: pin the average interaction effect — the
  [`ginteff`](https://github.com/jespernwulff/ginteff) estimand — with
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md);
  the internal AIE is unit-tested to reproduce `ginteff` exactly for all
  four variable-type pairs.
- **Three covariate routes**: parametric marginals with a Gaussian
  copula, pilot-data resampling
  ([`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md)),
  or a fitted pilot model
  ([`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md)).
- **Robustness sweeps**
  ([`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)):
  worst-case power and the insurance-premium sample size n_max over
  contextual assumptions, in the spirit of Hancock & Feng (2025).
- **Citable output**:
  [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  renders any result as a self-contained methods paragraph for grants
  and preregistrations.
- **Monte Carlo honesty**: every power estimate carries its MCSE; failed
  fits count against the claim, never silently dropped.

## Validation

The engine is anchored against independent implementations: the
two-proportion case reproduces
[`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html) and
Stata `power twoproportions`; APE estimates and delta-method SEs match
`marginaleffects` and Stata `margins, dydx()` to ~1e-6 on fixed
datasets; the AIE matches `ginteff` exactly. See the package tests and
the design document for details.

## Learn more

Four vignettes cover the workflows:
[`vignette("minimum-effect")`](https://jespernwulff.github.io/powerape/articles/minimum-effect.md),
[`vignette("classic-sesoi-detection")`](https://jespernwulff.github.io/powerape/articles/classic-sesoi-detection.md),
[`vignette("equivalence")`](https://jespernwulff.github.io/powerape/articles/equivalence.md),
and
[`vignette("interaction-effects")`](https://jespernwulff.github.io/powerape/articles/interaction-effects.md).
