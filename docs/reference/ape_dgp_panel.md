# Specify a panel DGP: correlated random effects probit/logit

Panel power analysis under the most general applied panel
binary-response model: a correlated random effects (CRE) specification
estimated by pooled probit/logit with Mundlak means and unit-clustered
standard errors (Wooldridge, 2010, ch. 15). Unobserved unit
heterogeneity
`c_i = xi'(centered observed unit means of time-varying regressors) + a_i`
may correlate with the regressors through the Mundlak part; the target
effect remains the ASF-based average partial effect, which the pooled
CRE estimator recovers even though its coefficients are attenuated.

## Usage

``` r
ape_dgp_panel(
  model = c("probit", "logit"),
  focal,
  moderator = NULL,
  covariates = list(),
  n_periods,
  rho = 0,
  cre_share = 0,
  correlation = NULL,
  baseline,
  signal = 0,
  n_int = 1e+05,
  seed_int = 20260814L
)
```

## Arguments

- model:

  `"probit"` (default) or `"logit"`.

- focal, covariates:

  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  objects; use their `icc` argument to set within-unit persistence
  (default 0; `icc = 1` = time-constant).

- moderator:

  Optional binary
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  for panel AIE designs: the index gains the moderator and a
  focal-by-moderator interaction, the estimand becomes the AIE (pin it
  with
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)),
  and the estimating model gains the interaction's own Mundlak mean
  whenever the product is time-varying. Currently binary focal x binary
  moderator only.

- n_periods:

  Number of periods per unit (balanced; \>= 2).

- rho:

  Latent unit-effect share `Var(c) / (Var(c) + 1)` — the xtprobit-style
  rho. 0 = no unobserved heterogeneity.

- cre_share:

  Fraction of `Var(c)` loaded (with equal standardized weights) on the
  centered observed unit means of the time-varying regressors; the
  remainder is the independent component `a_i`. `cre_share = 0` gives a
  pure random-effects world; the estimator includes Mundlak means either
  way.

- correlation:

  `NULL` (independence), a scalar exchangeable correlation, or a full
  correlation matrix over `(focal, moderator, covariates)` –
  Gaussian-copula dependence.

- baseline:

  Baseline outcome rate in (0, 1): the average probability `P(Y = 1)`
  with the focal (and moderator, if any) set to their reference values
  (0 for binary, the mean for continuous), covariates at their
  population distribution. The intercept is calibrated to reproduce it.

- signal:

  Nuisance-signal strength: the McKelvey-Zavoina-style latent
  pseudo-R-squared of the nuisance index, i.e.
  `Var(z'gamma) / (Var(z'gamma) + 1)`. Equal standardized weights,
  scaled to hit `signal`. Ignored when `gamma` is supplied.

- n_int, seed_int:

  Size (in **units**, each contributing `n_periods` rows) and seed of
  the deterministic integration draw used for calibration and inversion.

## Value

An object of class `powerape_dgp` (route `"panel"`). Pin the effect with
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
(no moderator) or
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
(with moderator); the empirical/pilot routes are not yet available for
panel DGPs.

## Details

Sample-size arguments of
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
and
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
refer to the **number of units (clusters)** for panel DGPs; each unit
contributes `n_periods` observations.

## Examples

``` r
# \donttest{
d <- ape_dgp_panel(
  focal = pa_var("treat", "binary", p = 0.5, icc = 1),  # unit-level treatment
  covariates = list(pa_var("size", "normal", icc = 0.6)),
  n_periods = 4, rho = 0.3, cre_share = 0.5,
  baseline = 0.30, signal = 0.10
)
d <- set_ape(d, target = 0.08)
ape_power(d, n = 150, claim = "detect", nsim = 300, seed = 1)  # n = units
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 150 units x 4 periods (600 obs), assumed true APE +0.0800, 95% CI, nsim = 300
#>   power = 0.490 (MCSE 0.029)
#>   outcomes: detect 0.490 | inconclusive 0.510 | failed 0.000 
# }
```
