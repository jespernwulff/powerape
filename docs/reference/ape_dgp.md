# Specify the data-generating process for an APE/AIE power analysis

Describes the world the study will sample from: an index model
(probit/logit), a binary or continuous (normal) focal variable, an
optional moderator (for average-interaction-effect designs), covariates
with a Gaussian-copula dependence structure, and an outcome anchored on
the friendly scale – the baseline outcome rate plus the strength of the
nuisance covariates. Effect coefficients are set later, in effect-size
units, by
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
(no moderator) or
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
(with moderator).

## Usage

``` r
ape_dgp(
  model = c("probit", "logit"),
  focal,
  moderator = NULL,
  covariates = list(),
  correlation = NULL,
  baseline,
  signal = 0,
  gamma = NULL,
  n_int = 4e+05,
  seed_int = 20260812L
)
```

## Arguments

- model:

  `"probit"` (default) or `"logit"`.

- focal:

  A
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md),
  `"binary"` or `"normal"`. For a continuous focal the APE is the
  average derivative, and reference = its mean.

- moderator:

  Optional
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md).
  When present, the index includes the moderator and a
  focal-by-moderator interaction, the estimand is the average
  interaction effect (AIE, as in the `ginteff` package), and the DGP
  must be pinned with
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md).

- covariates:

  List of
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  objects (may be empty).

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

- gamma:

  Optional explicit nuisance coefficients (length =
  `length(covariates)`); overrides `signal`.

- n_int, seed_int:

  Size and seed of the deterministic Monte Carlo integration draw used
  for calibration and inversion.

## Value

An object of class `powerape_dgp`. Pin the effect with
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
/
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md),
then use
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
or
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md).

## Examples

``` r
d <- ape_dgp(
  model = "probit",
  focal = pa_var("treat", "binary", p = 0.5),
  covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
  baseline = 0.30, signal = 0.15
)
d
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50)
#>   covariates: age (normal, mean 45.00, sd 12.00) 
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.150
#>   target effect not set yet - call set_ape() or set_aie().
```
