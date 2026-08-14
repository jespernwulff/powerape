# Package index

## Specify the world

Describe the data-generating process the study will sample from.

- [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  : Declare a variable of the data-generating process
- [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
  : Specify the data-generating process for an APE/AIE power analysis
- [`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md)
  : Specify a DGP with covariates resampled from pilot data
- [`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md)
  : Specify a DGP from a fitted pilot model
- [`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)
  : Specify a panel DGP: correlated random effects probit/logit
- [`ape_dgp_iv()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_iv.md)
  : Specify an IV DGP: endogenous focal, control-function probit

## Pin the effect

Solve coefficients from targets in APE/AIE units.

- [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
  : Pin the true APE of a DGP (inversion step)
- [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
  : Pin the true AIE of a moderated DGP (inversion step)
- [`true_ape()`](https://jespernwulff.github.io/powerape/reference/true_ape.md)
  : Recompute the true APE implied by a pinned DGP
- [`true_aie()`](https://jespernwulff.github.io/powerape/reference/true_aie.md)
  : Recompute the true AIE implied by a pinned DGP

## Power

Simulated power, curves, and required sample sizes.

- [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
  : Simulated power for an APE/AIE claim at a given sample size
- [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md)
  [`plot(`*`<powerape_curve>`*`)`](https://jespernwulff.github.io/powerape/reference/ape_curve.md)
  : Power curve over a grid of sample sizes
- [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  : Required sample size for a target power

## Robustness and reporting

- [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
  : Robustness sweep over contextual assumptions (worst-case power and
  n_max)
- [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  : Render a power analysis as a citable methods paragraph
