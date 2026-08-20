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
## ...or, when n is fixed by the archive, the smallest effect it can find:
ape_mde(d, n = 1500, claim = "minimum", sesoi = 0.05)
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
- **Panel designs**
  ([`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)):
  correlated random effects probit — Mundlak means, unit-clustered SEs,
  per-variable within-unit persistence via `pa_var(icc = )`, and sample
  sizes counted in **units (clusters)**. The target stays the ASF-based
  APE, which the pooled CRE estimator recovers; binary-by-binary panel
  AIE designs included. Probit-only by design: a pooled CRE logit would
  estimate the ASF only approximately, so `model = "logit"` is refused
  with a teaching error.
- **IV designs**
  ([`ape_dgp_iv()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_iv.md)):
  endogenous focal variables via control-function probit — the estimator
  of Stata’s `cfprobit`, replicated exactly with its stacked
  no-bootstrap standard errors (robust and clustered). Friendly knobs:
  `endogeneity` (error correlation) and `iv_strength` (first-stage share
  of the focal’s variance); endogenous-interaction AIE designs include
  the instrument-by-moderator first-stage terms.
- **Robust standard errors** for the standard routes via `se = "robust"`
  (HC0 sandwich); panel and IV routes carry their own clustered /
  stacked inference.
- **Three covariate routes**: parametric marginals with a Gaussian
  copula, pilot-data resampling
  ([`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md)),
  or a fitted pilot model
  ([`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md)).
- **Minimum detectable effects**
  ([`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)):
  when n is fixed — an archive of so many firms, a panel of so many
  units and waves — find the smallest APE or AIE the design reliably
  concludes, with the same conservative confirmation stage as
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md);
  the equivalence variant returns the tightest establishable margin.
- **Robustness sweeps**
  ([`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)):
  worst-case power and the insurance-premium sample size n_max over
  contextual assumptions, in the spirit of Hancock & Feng (2025) — and
  `mode = "mde"` for the minimum detectable effect under the least
  favorable assumptions.
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
datasets; the AIE matches `ginteff` exactly; simulated power matches
exact finite-sample enumeration in the saturated case and TOSTER under
matched conventions. Panel designs match Stata’s clustered `margins`
(APE and AIE double difference) to all printed digits on fixed panels,
the clustered sandwich matches
[`sandwich::vcovCL`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
exactly, and the engine reproduces the Donner-Klar cluster design
effect. The control-function route reproduces Stata’s `cfprobit` across
all four first-stage models — coefficients, stacked robust/clustered
SEs, and `margins` ASF effects to ~1e-6 — including the
endogenous-interaction design run with Stata’s own syntax. See the
package tests and the validation battery for details.

## Learn more

Six vignettes cover the workflows:
[`vignette("minimum-effect")`](https://jespernwulff.github.io/powerape/articles/minimum-effect.md),
[`vignette("classic-sesoi-detection")`](https://jespernwulff.github.io/powerape/articles/classic-sesoi-detection.md),
[`vignette("equivalence")`](https://jespernwulff.github.io/powerape/articles/equivalence.md),
[`vignette("interaction-effects")`](https://jespernwulff.github.io/powerape/articles/interaction-effects.md),
[`vignette("panel-designs")`](https://jespernwulff.github.io/powerape/articles/panel-designs.md),
and
[`vignette("iv-designs")`](https://jespernwulff.github.io/powerape/articles/iv-designs.md).
