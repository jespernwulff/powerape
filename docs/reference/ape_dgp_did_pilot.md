# A difference-in-differences design from pre-treatment pilot data

Estimates, from data observed **before** the treatment date, the inputs
of
[`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md)
that such data can inform, and builds the design:

## Usage

``` r
ape_dgp_did_pilot(data, outcome, group, id = NULL, period = NULL, trend, ...)
```

## Arguments

- data:

  A data frame of pre-treatment observations.

- outcome, group:

  Names of the binary outcome and of the indicator of the group that
  will be treated (0/1 or logical).

- id:

  Optional name of the unit identifier. With it the design is a panel
  and `rho` is estimated (units need two or more periods); without it
  the design is repeated cross-sections.

- period:

  Optional name of the period variable; used to report the control
  group's change between pre-treatment periods.

- trend:

  The control group's change in the outcome rate from before to after
  treatment, in probability points. Required.

- ...:

  Further arguments of
  [`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md),
  such as `model`, `n_pre`, `n_post`, `trends`, `analysis`,
  `covariates`. Naming `baseline`, `gap`, `p_treated`, `rho` or `design`
  here overrides the estimate. A `subgroup` is not supported.

## Value

A `powerape_dgp` as from
[`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md),
with the estimates and the pilot's size in `$pilot`. Pin the effect with
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md).

## Details

- `baseline`: the control group's outcome rate;

- `gap`: the treated group's rate minus the control group's;

- `p_treated`: the share of units (rows, without `id`) in the treated
  group;

- `rho`, for a panel: the latent within-unit correlation of the route's
  Gaussian copula, by the method of moments (the value at which the
  model's probability of the outcome in two periods of one unit equals
  the observed share of such pairs).

All pre-treatment periods are pooled, because the route holds rates flat
before the treatment date. Pre-treatment data cannot tell how the
control group will change from before to after treatment, so `trend`
must be supplied; when the pilot has several periods, the control
group's average change from one to the next is reported as a guide.
Covariates are not taken from the pilot: pass them as
[`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
objects.

The estimates are used as if they were known. They are not: vary them
with
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
(`baseline`, `gap`, `rho`, `p_treated`), in particular for a small
pilot.

## Examples

``` r
## a pilot panel: 300 units observed in two pre-treatment periods
set.seed(1)
u <- rnorm(300)
treated <- rbinom(300, 1, 0.5)
pilot <- data.frame(
  firm = rep(1:300, 2), year = rep(1:2, each = 300),
  treated = rep(treated, 2),
  y = rbinom(600, 1, pnorm(-0.7 + 0.3 * rep(treated, 2) + 0.8 * rep(u, 2))))
d <- suppressWarnings(
  ape_dgp_did_pilot(pilot, outcome = "y", group = "treated", id = "firm",
                    period = "year", trend = 0.05, model = "logit"))
d <- set_att(d, 0.08)
d
#> powerape DGP -- logit, difference-in-differences (two periods)
#>   data: panel, both periods per unit; rho = 0.40 (latent within-unit correlation); unit-clustered SEs
#>   analysis: pooled logit; ATT = average partial effect of the treatment indicator over treated post-period rows
#>   treated group: 49% of units
#>   covariates: none
#>   untreated trends: parallel on the logit index (the scale the logit analysis assumes)
#>   rates before -> after: control 0.322 -> 0.372; treated 0.385 -> 0.519 (0.439 without treatment)
#>   pilot: baseline, gap, treated share, rho estimated from 600 pre-treatment rows (300 units); control group's change between pilot periods +0.000
#>   true ATT pinned at +0.0800 (treatment coefficient = 0.3213)
#>   the same four rates read under linear-parallel trends: +0.0834 (their cross difference)
```
