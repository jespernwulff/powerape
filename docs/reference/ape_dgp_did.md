# Specify a difference-in-differences design with a binary outcome

Describes a difference-in-differences study with one treatment date: a
treated and a control group observed before and after the treated group
receives the intervention, in one period each by default or in `n_pre`
periods before and `n_post` periods after. The target is the **average
treatment effect on the treated (ATT)** after treatment, in probability
points, pinned with
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md).
The outcome is regressed on the group indicator, period indicators, the
treatment indicator (treated group after the treatment date), and any
covariates. Under `analysis = "index"` (default) that regression is a
pooled probit or logit and the ATT is the average partial effect of the
treatment indicator over the treated post-period observations
(Wooldridge, 2023), with a delta-method standard error. Under
`analysis = "linear"` it is a linear probability model fitted by
ordinary least squares and the ATT is the coefficient of the treatment
indicator, the standard difference-in-differences estimate.

## Usage

``` r
ape_dgp_did(
  model = c("probit", "logit"),
  baseline,
  gap = 0,
  trend = 0,
  p_treated = 0.5,
  design = c("panel", "cross-sections"),
  rho = 0,
  p_post = NULL,
  trends = c("index", "linear"),
  analysis = c("index", "linear"),
  n_pre = 1,
  n_post = 1,
  subgroup = NULL,
  covariates = list(),
  correlation = NULL,
  signal = 0,
  gamma = NULL,
  n_int = 4e+05,
  seed_int = 20260812L
)
```

## Arguments

- model:

  `"probit"` (default) or `"logit"`.

- baseline:

  Outcome rate of the control group before treatment, in (0, 1).

- gap:

  Treated-group rate minus control-group rate before treatment
  (probability points; default 0).

- trend:

  Change in the control-group rate from before to after the treatment
  date (probability points; default 0).

- p_treated:

  Share of units (or observations) in the treated group.

- design:

  `"panel"` (default) or `"cross-sections"`.

- rho:

  Panel designs: latent correlation between any two outcomes of a unit,
  in `[0, 1)`; with covariates, the correlation given the covariates
  (section 'Data structure').

- p_post:

  Cross-section designs: share of observations drawn after the treatment
  date (default `n_post / (n_pre + n_post)`).

- trends:

  Scale on which the untreated trends are parallel in the simulated
  world: `"index"` (default) or `"linear"`.

- analysis:

  The analysis the simulated studies receive: `"index"` (default; pooled
  probit/logit, ATT as an average partial effect) or `"linear"` (linear
  probability model, ATT as the coefficient of the treatment indicator,
  robust standard errors).

- n_pre, n_post:

  Number of periods before and after the treatment date (default one
  each).

- subgroup:

  Optional binary
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  splitting both groups into two subgroups, e.g.
  `pa_var("female", "binary", p = 0.4)`; see the section 'Comparing two
  subgroups'.

- covariates:

  List of
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  objects (may be empty); time-constant.

- correlation:

  `NULL`, a scalar, or a correlation matrix over `(group, covariates)`,
  or over `(group, subgroup, covariates)` when there is a subgroup –
  Gaussian-copula dependence, as in
  [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md).

- signal, gamma:

  Nuisance-signal strength or explicit nuisance coefficients, as in
  [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md).

- n_int, seed_int:

  Size and seed of the deterministic integration draw used for
  calibration when there are covariates.

## Value

An object of class `powerape_dgp` (route `"did"`). Pin the effect with
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md),
then use
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
or
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md).

## Details

The world is described on the probability scale by three rates:
`baseline` (control group, before), `gap` (treated minus control,
before), and `trend` (the control group's change from before to after).
With several periods the rates are the same in every period before and
in every period after the treatment date, and so is the effect. What the
treated group's rate would have been after the intervention date
*without* treatment then depends on the scale on which the two groups'
untreated trends are parallel, set by `trends`:

- `"index"` (default): parallel on the probit/logit index, the
  assumption under which the probit/logit analysis identifies the ATT
  (Wooldridge, 2023; Puhani, 2012).

- `"linear"`: parallel in probabilities, the assumption under which the
  linear analysis identifies the ATT. With covariates the rates are
  parallel on average over each group's covariates.

Each analysis is unbiased only in the world whose scale it assumes. When
`trends` and `analysis` name the same scale,
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
reports power. When they differ, the analysis centres on another number
than the ATT:
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
and
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md)
then report a **rejection rate** together with the mean of the
estimates, and
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
and
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
refuse the design. One case is exempt: with `gap = 0` or `trend = 0` and
no covariates that move the outcome, parallel trends mean the same on
both scales, the two worlds are the same world, and either analysis is
treated as valid. With two periods and no covariates the data cannot
tell the two scales apart (Wooldridge, 2023), and the same four cell
rates imply different effects under each: the printed design shows both
readings. The linear analysis always centres on the cross difference of
the four rates.

## Data structure

`design = "panel"` follows the same units over all periods (a balanced
panel); `n` in the power functions counts units, and standard errors are
clustered by unit. Outcomes of one unit are dependent through a Gaussian
copula with latent correlation `rho` that leaves the marginal rates
exact under either link; under probit this is the random-effects probit
with `rho = Var(c) / (Var(c) + 1)`. With covariates, `rho` is the
dependence that is left **given the covariates**, as in a random-effects
probit with regressors: covariates are constant over time, so they make
a unit's outcomes more alike on top of `rho`. Under probit with normal
covariates the raw latent correlation between two outcomes of a unit is
`(rho + s2) / (1 + s2)` with `s2 = signal / (1 - signal)`; to plan from
a raw within-unit correlation `r`, use
`rho = (r - signal) / (1 - signal)`. Under logit the covariates add less
and there is no closed form;
[`ape_dgp_did_pilot()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did_pilot.md)
solves for the value under either link. `design = "cross-sections"`
draws independent observations, each after the treatment date with
probability `p_post` and then in one of that block's periods with equal
probability; `n` counts all observations, and `se` in
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
chooses model-based or robust standard errors for the probit/logit
analysis.

## More than two periods

The analysis includes an indicator for every period but the first and
one treatment indicator, so it estimates one effect, not one per period.
In a panel the dependence between two outcomes of a unit is the same
however far apart the periods are (one shared unit component). That is
the favourable case for adding periods: averaging over periods then
removes period-specific noise only, and the gain is larger than under
dependence that fades with the distance between periods, which this
version does not offer. As an example, with rates of .20 and .30 in the
control group, .40 and .60 in the treated group and `rho = 0.5`, two
periods on each side instead of one cut the variance of the linear panel
estimate to 51%; if the latent correlation halved with every period of
distance it would be 71%. Treat the precision bought by extra periods as
an upper bound.

## Comparing two subgroups

With `subgroup`, both subgroups are treated in the treated group and the
target is the **difference between the two subgroups' effects on the
treated**, in probability points (the estimand of a triple difference
that compares, say, the effect on women with the effect on men; Caron,
2025). `baseline`, `gap` and `trend` may then hold two numbers, for the
reference subgroup (`subgroup = 0`) and for `subgroup = 1`, and
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md)
takes the two effects. The probit/logit analysis is fully interacted
with the subgroup and the estimate is the difference of the two
subgroups' average partial effects of treatment; the linear analysis
estimates the coefficient of the subgroup-by- treatment product. Trends
are parallel within each subgroup, on the scale `trends` names. Three
things differ from a single effect. The same eight rates imply different
subgroup differences on the two scales, and the printed design shows
both. A difference of zero in probability points is not a zero
interaction coefficient when the subgroups have different rates. And
with equal subgroups a difference as large as the effect itself needs
about four times the units the effect needs.
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
searches the difference, holding the reference subgroup's effect at its
pinned value.

## The linear analysis

A linear probability model is heteroskedastic by construction, so its
standard errors are always robust and the `se` argument of the power
functions is ignored: HC1 in cross-sections (as Stata's
`regress, vce(robust)`) and clustered by unit in panels (with the
small-sample factor of Stata's `regress, vce(cluster)`). Ordinary least
squares has no separation, so no simulated study is discarded for an
all-zero or all-one cell; the sparse-cell warning still applies. In a
balanced panel a time-constant covariate changes neither the estimate
nor its clustered standard error beyond the degrees of freedom (its
contribution cancels in each unit's change), so adjusting for covariates
buys the linear panel analysis no precision; in cross-sections it does.
A panel design with a larger `signal` at the same `rho` still has more
power under either analysis, because `rho` is the dependence given the
covariates and the design's units are then more alike over time (section
'Data structure'), not because of the adjustment.

## Interval coverage in small samples

Under the index analysis the effect on the treated is a nonlinear
function of the estimated rates, and its delta-method interval (the one
`margins` reports in Stata for the same model) covers slightly less than
its nominal level when the group-by-period cells are small. In the
design without covariates the coverage can be computed exactly. For a
95% interval at a true effect of zero, with rates between 25% and 44%,
it is 94.0% with 125 observations per cell, 94.5% with 250, 94.7% with
500 and 94.9% with 1,000 under logit (probit is up to 0.2 points closer
to 95%): the shortfall halves when the cells double. A detection claim
is therefore made for a true zero effect in 6.0% of studies at 125 per
cell and in 5.3% at 500, not in 5%. With covariates and in panels the
simulated coverage at about 500 per cell is similar (94.4% to 94.7%).
The simulated studies use this interval, so the power reported is the
power of the analysis a study will run and needs no correction; with
small cells, part of it is this excess. The interval of the linear
analysis is close to nominal at the same sizes (94.9% to 95.0% for the
same rates).

A second effect concerns both analyses. Robust and cluster-robust
standard errors are too small when one group has few independent units,
and more periods per unit do not help, because the units are what is
independent. With 10% of the units treated and a true effect of zero,
the 5% test of the linear analysis rejected in 12%, 10%, 7% and 6% of
simulated two-period panels with 4, 10, 20 and 40 treated units, and in
17%, 10% and 6% of panels with four periods on each side and 4, 10 and
40 treated units; the probit/logit analysis behaved alike (9%, 8% and 6%
with 10, 20 and 40 treated units). With groups of equal size the excess
of the linear analysis is small (5.7% with 20 units a group in a
two-period panel, 6.9% with four periods on each side, 5% to 6% from 50
a group). The power functions warn when fewer than 30 units are expected
in the smaller group; in repeated cross-sections, when fewer than 30
observations are expected in the smallest group-by-before/after cell;
with a subgroup, in the smallest group-by-subgroup cell. Power reported
for such a design contains the excess.

## Limits of this version

One treated group, one treatment date, treatment assigned at the unit
level, balanced panels. Staggered adoption, effects that differ by
period, pre-trend tests, rates that drift within the before or after
block, and treatment assigned to whole clusters are not available. All
power functions accept these designs, including
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
(which can also vary the scale on which trends are parallel),
[`power_at()`](https://jespernwulff.github.io/powerape/reference/power_at.md),
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
and
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md).
Functions that need power stop for a design whose analysis assumes
another scale than its world, because its result is a rejection rate.
[`ape_dgp_did_pilot()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did_pilot.md)
builds the design from data observed before the treatment date.

## References

Caron, L. (2025). Triple difference designs with heterogeneous treatment
effects. arXiv:2502.19620.

Puhani, P. A. (2012). The treatment effect, the cross difference, and
the interaction term in nonlinear "difference-in-differences" models.
*Economics Letters, 115*(1), 85-87.

Wooldridge, J. M. (2023). Simple approaches to nonlinear
difference-in-differences with panel data. *The Econometrics Journal,
26*(3), C31-C66.

## Examples

``` r
d <- ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20,
                 trend = 0.10, rho = 0.4)
d <- set_att(d, target = 0.07)
d
#> powerape DGP -- logit, difference-in-differences (two periods)
#>   data: panel, both periods per unit; rho = 0.40 (latent within-unit correlation); unit-clustered SEs
#>   analysis: pooled logit; ATT = average partial effect of the treatment indicator over treated post-period rows
#>   treated group: 50% of units
#>   covariates: none
#>   untreated trends: parallel on the logit index (the scale the logit analysis assumes)
#>   rates before -> after: control 0.200 -> 0.300; treated 0.400 -> 0.603 (0.533 without treatment)
#>   true ATT pinned at +0.0700 (treatment coefficient = 0.2858)
#>   the same four rates read under linear-parallel trends: +0.1033 (their cross difference)

## the linear probability analysis in a world with linear-parallel trends
dl <- ape_dgp_did(baseline = 0.20, gap = 0.20, trend = 0.10,
                  trends = "linear", analysis = "linear")
set_att(dl, target = 0.10)
#> powerape DGP -- probit, difference-in-differences (two periods)
#>   data: panel, both periods per unit; rho = 0.00 (latent within-unit correlation); unit-clustered SEs
#>   analysis: linear probability model (OLS); ATT = coefficient of the treatment indicator; unit-clustered SEs
#>   treated group: 50% of units
#>   covariates: none
#>   untreated trends: parallel in probabilities (linear; the scale the linear probability analysis assumes)
#>   rates before -> after: control 0.200 -> 0.300; treated 0.400 -> 0.600 (0.500 without treatment)
#>   true ATT pinned at +0.1000 (treatment coefficient = 0.2533)
#>   the same four rates read on the probit index: +0.0745 (covariate-free reading)
```
