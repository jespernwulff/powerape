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
  : Specify a panel DGP: correlated random effects probit
- [`ape_dgp_iv()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_iv.md)
  : Specify an IV DGP: endogenous focal, control-function probit
- [`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md)
  : Specify a difference-in-differences design with a binary outcome
- [`ape_dgp_did_pilot()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did_pilot.md)
  : A difference-in-differences design from pre-treatment pilot data
- [`set_analysis()`](https://jespernwulff.github.io/powerape/reference/set_analysis.md)
  : Choose the analysis the simulated studies receive

## Pin the effect

Solve coefficients from targets in APE/AIE/ATT units.

- [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
  : Pin the true APE of a DGP (inversion step)
- [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
  : Pin the true AIE of a moderated DGP (inversion step)
- [`set_aie3()`](https://jespernwulff.github.io/powerape/reference/set_aie3.md)
  : Pin the three-way interaction effect of a design with two moderators
- [`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md)
  : Pin the true effect on the treated of a difference-in-differences
  design
- [`true_ape()`](https://jespernwulff.github.io/powerape/reference/true_ape.md)
  : Recompute the true APE implied by a pinned DGP
- [`true_aie()`](https://jespernwulff.github.io/powerape/reference/true_aie.md)
  : Recompute the true AIE implied by a pinned DGP
- [`true_aie3()`](https://jespernwulff.github.io/powerape/reference/true_aie3.md)
  : Recompute the true three-way interaction effect implied by a pinned
  DGP
- [`true_att()`](https://jespernwulff.github.io/powerape/reference/true_att.md)
  : Recompute the true effect on the treated implied by a pinned design

## Power

Simulated power, curves, required sample sizes, detectable effects.

- [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
  : Simulated power for an APE/AIE claim at a given sample size
- [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md)
  [`plot(`*`<powerape_curve>`*`)`](https://jespernwulff.github.io/powerape/reference/ape_curve.md)
  : Power curve over a grid of sample sizes
- [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  : Required sample size for a target power
- [`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
  : Minimum detectable APE/AIE at a given sample size

## Justified error rates

Alpha by sample-size rule or by error costs, from one stored simulation.

- [`power_at()`](https://jespernwulff.github.io/powerape/reference/power_at.md)
  : Power at other error rates from one stored simulation
- [`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
  [`plot(`*`<powerape_alpha>`*`)`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
  : Justify the error rate: minimize or balance weighted Type I and Type
  II errors
- [`alpha_frontier()`](https://jespernwulff.github.io/powerape/reference/alpha_frontier.md)
  : Closed-form justified alpha at the detectability frontier

## Robustness and reporting

- [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
  : Robustness sweep over contextual assumptions (worst-case power and
  n_max)
- [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  : Render a power analysis as a citable methods paragraph
