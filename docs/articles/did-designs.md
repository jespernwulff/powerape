# Difference-in-differences designs

## The question

A treated group and a control group are observed before and after one
treatment date, and the outcome is binary. How many units does the study
need to show that the treatment changed the outcome rate of the treated?

The target here is the **average treatment effect on the treated (ATT)
in probability points**: the treated group’s outcome rate after
treatment, minus the rate it would have had without treatment.
[`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md)
describes such a study,
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md)
pins the effect, and the usual power functions do the rest.

## Describing the study in rates

The design is stated in three rates and one share: the control group’s
rate before treatment (`baseline`), the gap between the groups before
treatment (`gap`), the control group’s change from before to after
(`trend`), and the share of treated units (`p_treated`).

- R
- Stata

``` r

library(powerape)
```

``` r

d <- ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10,
                 design = "panel", rho = 0.5)
d <- set_att(d, 0.067)
d
#> powerape DGP -- logit, difference-in-differences (two periods)
#>   data: panel, both periods per unit; rho = 0.50 (latent within-unit correlation); unit-clustered SEs
#>   analysis: pooled logit; ATT = average partial effect of the treatment indicator over treated post-period rows
#>   treated group: 50% of units
#>   covariates: none
#>   untreated trends: parallel on the logit index (the scale the logit analysis assumes)
#>   rates before -> after: control 0.200 -> 0.300; treated 0.400 -> 0.600 (0.533 without treatment)
#>   true ATT pinned at +0.0670 (treatment coefficient = 0.2733)
#>   the same four rates read under linear-parallel trends: +0.1003 (their cross difference)
```

``` stata
powerape dgp d, did model(logit) baseline(0.20) gap(0.20) trend(0.10)      ///
    design(panel) rho(0.5)

powerape setatt d, target(0.067)
```

``` stata-output
powerape DGP -- logit, difference-in-differences (two periods)
  data: panel, both periods per unit; rho = 0.50 (latent within-unit correlation); unit-clustered SEs
  analysis: pooled logit; ATT = average partial effect of the treatment indicator over treated
    post-period rows
  treated group: 50% of units
  covariates: none
  untreated trends: parallel on the logit index (the scale the logit analysis assumes)
  rates before -> after: control 0.200 -> 0.300; treated 0.400 -> 0.600 (0.533 without treatment)
  true ATT pinned at +0.0670 (treatment coefficient = 0.2733)
  the same four rates read under linear-parallel trends: +0.1003 (their cross difference)
```

`design = "panel"` follows the same units in both periods; `n` then
counts units, standard errors are clustered by unit, and `rho` is the
dependence between a unit’s outcomes (in a design with covariates: the
dependence that is left given the covariates).
`design = "cross-sections"` draws new observations in each period.

## Which effect do four rates show?

The print ends with a second reading of the same four rates. A control
group that goes from 20% to 30% and a treated group that goes from 40%
to 60% show an effect of 10 points if the groups would have moved in
parallel **in probabilities**, and of 6.7 points if they would have
moved in parallel **on the logit index** (7.4 points on the probit
index). With a binary outcome the two cannot both hold unless the groups
start at the same rate or the control group does not change, and two
periods of data cannot tell them apart.

The package therefore keeps two things separate:

- `trends`: the scale on which the untreated trends are parallel **in
  the simulated world**: `"index"` (the default) or `"linear"`;
- `analysis`: what each simulated study is analysed with: `"index"`, a
  pooled logit or probit whose effect on the treated is the average
  partial effect of the treatment indicator over the treated
  post-treatment observations (Wooldridge, 2023; Puhani, 2012), or
  `"linear"`, a linear probability model whose effect is the coefficient
  of the treatment indicator.

An analysis is unbiased in the world whose scale it assumes. When the
two agree, the functions report power:

- R
- Stata

``` r

ape_power(d, n = 3000, claim = "detect", nsim = 300, seed = 1)
#> powerape -- detection claim (CI excludes 0, directional)
#>   logit, n = 3000 units x 2 periods (6000 obs), assumed true ATT +0.0670, 95% CI (alpha = 0.05), nsim = 300
#>   power = 0.797 (MCSE 0.023)
#>   outcomes: detect 0.797 | inconclusive 0.203 | failed 0.000
```

``` r

ape_n(d, power = 0.80, claim = "detect", nsim = 1000, seed = 2)
```

``` stata
powerape power d, n(3000) claim(detect) nsim(300) seed(1)
```

``` stata-output
powerape -- detection claim (CI excludes 0, directional)
  logit, n = 3000 units x 2 periods (6000 obs), assumed true ATT + 0.0670, 95% CI (alpha = .05),
    nsim = 300
  fit: logit, vce(cluster) via engine(stata), effects(closed)
  power = 0.820 (MCSE 0.022)
  outcomes: detect 0.820 | inconclusive 0.180 | failed 0.000
```

``` stata
powerape n d, power(0.80) claim(detect) nsim(1000) seed(2) engine(mata)
```

When they differ, the result is labelled a rejection rate, and the mean
of the estimates shows what the analysis centres on:

- R
- Stata

``` r

d_mis <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20,
                             trend = 0.10, rho = 0.5,
                             trends = "linear", analysis = "index"), 0.10)
ape_power(d_mis, n = 1000, claim = "detect", nsim = 300, seed = 3)
#> powerape -- detection claim (CI excludes 0, directional)
#>   logit, n = 1000 units x 2 periods (2000 obs), assumed true ATT +0.1000, 95% CI (alpha = 0.05), nsim = 300
#>   rejection rate = 0.367 (MCSE 0.028) -- not power: the simulated world has linear-parallel trends,
#>     the logit analysis assumes index-parallel trends; its estimates average +0.0658, the true ATT is +0.1000
#>   outcomes: detect 0.367 | inconclusive 0.633 | failed 0.000
```

``` stata
powerape dgp d_mis, did model(logit) baseline(0.20) gap(0.20) trend(0.10)  ///
    rho(0.5) trends(linear) analysis(index)
powerape setatt d_mis, target(0.10)

powerape power d_mis, n(1000) claim(detect) nsim(300) seed(3)
```

``` stata-output
powerape -- detection claim (CI excludes 0, directional)
  logit, n = 1000 units x 2 periods (2000 obs), assumed true ATT + 0.1000, 95% CI (alpha = .05),
    nsim = 300
  fit: logit, vce(cluster) via engine(stata), effects(closed)
  rejection rate = 0.393 (MCSE 0.028) -- not power: the simulated world has linear-parallel trends,
    the logit analysis assumes index-parallel trends; its estimates average +0.0675, the true ATT is +0.1000
  outcomes: detect 0.393 | inconclusive 0.607 | failed 0.000
```

Here the world has an effect of 10 points, and the logit analysis
averages about 6.7, the logit reading of the same rates.
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
and
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
stop for such a design, because no sample size delivers the claim for
the right reason.

## What if the trends are parallel on the other scale?

[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
re-runs the analysis over a grid of the design’s inputs. For these
designs the grid can include the scale itself:

- R
- Stata

``` r

ape_robust(d, n = 3000, claim = "detect", nsim = 200, seed = 4,
           vary = list(rho = c(0.3, 0.7), trends = c("index", "linear")),
           grid_points = 2, nmax = FALSE)
#> powerape robustness sweep -- detect claim, ATT, pin = ape
#>   n = 3000, nsim = 200 per scenario, 4 scenario(s) over: rho, trends
#>   note: in 2 scenario(s) (valid = FALSE) the analysis assumes another parallel-trends scale than the world has: `power` is a rejection rate there
#>     and `est_mean` shows where the estimates centre (the true ATT is +0.0670); they are left out of the range, the worst case and n_max.
#>   power range [0.825, 0.870]; worst scenario:
#>  rho trends implied_effect power       mcse separated   est_mean valid
#>  0.3  index          0.067 0.825 0.02686773         0 0.07174983  TRUE
#>   scenarios:
#>  rho trends implied_effect power       mcse separated   est_mean valid
#>  0.3  index          0.067 0.825 0.02686773         0 0.07174983  TRUE
#>  0.7  index          0.067 0.870 0.02378024         0 0.06493598  TRUE
#>  0.3 linear          0.067 0.230 0.02975735         0 0.03095584 FALSE
#>  0.7 linear          0.067 0.370 0.03413942         0 0.03440064 FALSE
```

``` stata
powerape robust d, n(3000) claim(detect) nsim(200) seed(4)                 ///
    vary(rho(0.3 0.7) trends(index linear)) gridpoints(2) nonmax engine(mata)
```

``` stata-output
powerape robustness sweep -- detect claim, ATT, pin = ape
  n = 3000, nsim = 200 per scenario, 4 scenario(s) over: rho, trends
  note: in 2 scenario(s) (valid = 0) the analysis assumes another parallel-trends scale than the
    world has: power is a rejection rate there
    and est_mean shows where the estimates centre (the true ATT is +0.0670); they are left out of
      the range, the worst case and n_max.
  power range [0.775, 0.915]; worst scenario: row 1
  fit: Mata engine (Fisher scoring)
 row          rho       trends    implied_eff    power     mcse separated  est_mean  valid
   1           .3        index         0.0670   0.7750   0.0295    0.0000    0.0696      1
   2           .7        index         0.0670   0.9150   0.0197    0.0000    0.0676      1
   3           .3       linear         0.0670   0.2500   0.0306    0.0000    0.0348      0
   4           .7       linear         0.0670   0.3750   0.0342    0.0000    0.0344      0
```

In every scenario the effect on the treated is pinned at the planning
value. Rows with `valid = FALSE` are worlds in which the analysis is
biased: their `power` is a rejection rate, `est_mean` shows where the
estimates centre, and they are left out of the worst case and of the
`n_max` search. The other inputs that can be varied are `baseline`,
`gap`, `trend`, `p_treated` and `signal`.

## More than two periods

`n_pre` and `n_post` give the number of periods before and after the
treatment date. The analysis then has an indicator for every period but
the first and one treatment indicator, so it estimates one effect.

- R
- Stata

``` r

d4 <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20,
                          trend = 0.10, rho = 0.5, n_pre = 2, n_post = 2),
              0.067)
ape_power(d4, n = 1500, claim = "detect", nsim = 300, seed = 5)
#> powerape -- detection claim (CI excludes 0, directional)
#>   logit, n = 1500 units x 4 periods (6000 obs), assumed true ATT +0.0670, 95% CI (alpha = 0.05), nsim = 300
#>   power = 0.807 (MCSE 0.023)
#>   outcomes: detect 0.807 | inconclusive 0.193 | failed 0.000
```

``` stata
powerape dgp d4, did model(logit) baseline(0.20) gap(0.20) trend(0.10)     ///
    rho(0.5) npre(2) npost(2)
powerape setatt d4, target(0.067)

powerape power d4, n(1500) claim(detect) nsim(300) seed(5)
```

``` stata-output
powerape -- detection claim (CI excludes 0, directional)
  logit, n = 1500 units x 4 periods (6000 obs), assumed true ATT + 0.0670, 95% CI (alpha = .05),
    nsim = 300
  fit: logit, vce(cluster) via engine(stata), effects(closed)
  power = 0.810 (MCSE 0.023)
  outcomes: detect 0.810 | inconclusive 0.190 | failed 0.000
```

Two limits matter when reading this. Rates are the same in every period
before and in every period after the treatment date, and so is the
effect. And the dependence between two outcomes of a unit is the same
however far apart the periods are. That is the favourable case for
adding periods, so the precision extra periods buy here is an upper
bound.

## Comparing the effect between two subgroups

With `subgroup`, both groups are split into two subgroups that are both
treated, and the target becomes the difference between their effects on
the treated.

- R
- Stata

``` r

ds <- ape_dgp_did(model = "logit", baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                  trend = c(0.10, 0.05), rho = 0.5,
                  subgroup = pa_var("female", "binary", p = 0.4))
ds <- set_att(ds, c(0.10, 0.03))
ds
#> powerape DGP -- logit, difference-in-differences with a subgroup (two periods)
#>   data: panel, both periods per unit; rho = 0.50 (latent within-unit correlation); unit-clustered SEs
#>   analysis: pooled logit, fully interacted with female; estimate = difference of the two subgroups' average partial effects of treatment
#>   treated group: 50% of units; subgroup female = 1: 40%
#>   covariates: none
#>   untreated trends: parallel within each subgroup, on the logit index
#>   rates before -> after, female = 0: control 0.300 -> 0.400; treated 0.400 -> 0.609 (0.509 without treatment)
#>   rates before -> after, female = 1: control 0.200 -> 0.250; treated 0.350 -> 0.448 (0.418 without treatment)
#>   effects on the treated: +0.1000 (female = 0), +0.0300 (female = 1); true difference (DATT) pinned at -0.0700
#>   the same eight rates read under linear-parallel trends: difference -0.0612 (their triple cross difference)
```

``` stata
powerape dgp ds, did model(logit) baseline(0.30 0.20) gap(0.10 0.15)       ///
    trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4))

powerape setatt ds, target(0.10 0.03)
```

``` stata-output
powerape DGP -- logit, difference-in-differences with a subgroup (two periods)
  data: panel, both periods per unit; rho = 0.50 (latent within-unit correlation); unit-clustered SEs
  analysis: pooled logit, fully interacted with female; estimate = difference of the two subgroups'
    average partial effects of treatment
  treated group: 50% of units; subgroup female = 1: 40%
  covariates: none
  untreated trends: parallel within each subgroup, on the logit index
  rates before -> after, female = 0: control 0.300 -> 0.400; treated 0.400 -> 0.609 (0.509 without treatment)
  rates before -> after, female = 1: control 0.200 -> 0.250; treated 0.350 -> 0.448 (0.418 without treatment)
  effects on the treated: +0.1000 (female = 0), +0.0300 (female = 1); true difference (DATT) pinned at -0.0700
  the same eight rates read under linear-parallel trends: difference -0.0612 (their triple cross difference)
```

A difference between two effects is far more expensive than one effect:
with equal subgroups and a difference as large as the effect itself, it
needs about four times the units.

## Starting from pilot data

Data from before the treatment date inform the control group’s rate, the
gap between the groups, the treated share and, in a panel, the
dependence within units.
[`ape_dgp_did_pilot()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did_pilot.md)
estimates these and builds the design. The control group’s change from
before to after is not in such data, so `trend` has to be supplied.

- R
- Stata

``` r

set.seed(1)
u <- rnorm(600)
treated <- rbinom(600, 1, 0.5)
pilot <- data.frame(
  firm = rep(1:600, 2), year = rep(1:2, each = 600), treated = rep(treated, 2),
  y = rbinom(1200, 1, pnorm(-0.7 + 0.3 * rep(treated, 2) + 0.8 * rep(u, 2))))
dp <- ape_dgp_did_pilot(pilot, outcome = "y", group = "treated", id = "firm",
                        period = "year", trend = 0.05, model = "logit")
set_att(dp, 0.08)
#> powerape DGP -- logit, difference-in-differences (two periods)
#>   data: panel, both periods per unit; rho = 0.43 (latent within-unit correlation); unit-clustered SEs
#>   analysis: pooled logit; ATT = average partial effect of the treatment indicator over treated post-period rows
#>   treated group: 48% of units
#>   covariates: none
#>   untreated trends: parallel on the logit index (the scale the logit analysis assumes)
#>   rates before -> after: control 0.320 -> 0.370; treated 0.341 -> 0.472 (0.392 without treatment)
#>   pilot: baseline, gap, treated share, rho estimated from 1200 pre-treatment rows (600 units); control group's change between pilot periods +0.010
#>   true ATT pinned at +0.0800 (treatment coefficient = 0.3268)
#>   the same four rates read under linear-parallel trends: +0.0814 (their cross difference)
```

``` stata
* The pilot rows are the ones the R code draws (Stata cannot draw R's random
* numbers): precompute/stata/did-pilot.csv in the R package's repository.
import delimited using "did-pilot.csv", encoding(utf-8) clear
powerape dgp dp, didpilot outcome(y) group(treated) id(firm) period(year)  ///
    trend(0.05) model(logit)

powerape setatt dp, target(0.08)
```

``` stata-output
powerape DGP -- logit, difference-in-differences (two periods)
  data: panel, both periods per unit; rho = 0.43 (latent within-unit correlation); unit-clustered SEs
  analysis: pooled logit; ATT = average partial effect of the treatment indicator over treated
    post-period rows
  treated group: 48% of units
  covariates: none
  untreated trends: parallel on the logit index (the scale the logit analysis assumes)
  rates before -> after: control 0.320 -> 0.370; treated 0.341 -> 0.472 (0.392 without treatment)
  pilot: baseline, gap, treated share, rho estimated from 1200 pre-treatment rows (600 units);
    control group's change between pilot periods +0.010
  true ATT pinned at +0.0800 (treatment coefficient = 0.3268)
  the same four rates read under linear-parallel trends: +0.0814 (their cross difference)
```

The estimates are used as if they were known. Sweep them with
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
before relying on a sample size. When covariates are added to the call,
the function returns the `rho` that is left given those covariates,
which is lower than the raw correlation in the pilot: the covariates are
constant over time and account for part of it.

## Reporting

[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
writes the methods paragraph: the groups and periods, the rates, the
scale on which trends are parallel, the analysis and its standard
errors.

- R
- Stata

``` r

pw <- ape_power(d, n = 3000, claim = "minimum", sesoi = 0.02, nsim = 300, seed = 6)
power_statement(pw)
#> We conducted a simulation-based power analysis for the average treatment
#> effect on the treated (ATT) in a difference-in-differences design, in
#> probability points, using the powerape package (version 1.12.2), following
#> the confidence-interval approach of Riesthuis (2024). The assumed
#> data-generating process was a difference-in-differences design with a treated
#> group (50% of the units) and a control group, observed in one period before
#> and one period after one common treatment date as a panel of the same units,
#> with a latent within-unit correlation of 0.50 shared by all periods (Gaussian
#> copula); outcome rates before and after treatment of 0.200 and 0.300 in the
#> control group and 0.400 and 0.600 in the treated group (0.533 without
#> treatment); untreated trends parallel on the logit index; no covariates. The
#> assumed true ATT was 0.067 (6.7 percentage points). At a sample size of n =
#> 3000 units observed in 2 periods (6,000 unit-period observations), simulated
#> power for the minimum-effect claim (the 90% confidence interval's lower bound
#> exceeding the smallest effect size of interest, 0.020; a one-sided test at
#> alpha = 5%) was 0.617 (Monte Carlo SE 0.028; 300 replications). Across
#> replications, the probability of concluding a meaningful effect was 0.617, of
#> detection without meaningfulness 0.200, of an inconclusive result 0.183, and
#> of equivalence 0.000. Each simulated study was analysed by a pooled logit
#> regression of the outcome on group, period indicators and the treatment
#> indicator, with the effect on the treated computed as the average partial
#> effect of the treatment indicator over the treated post-treatment
#> observations (Wooldridge, 2023), and delta-method Wald confidence intervals
#> from unit-clustered standard errors; replications without a
#> maximum-likelihood estimate, from a failed fit or from separation (a
#> group-by-before/after cell with no events or no non-events), counted against
#> the claim (0.0% of replications).
```

``` stata
powerape power d, n(3000) claim(minimum) sesoi(0.02) nsim(300) seed(6)

powerape statement
```

``` stata-output
We conducted a simulation-based power analysis for the average treatment
effect on the treated (ATT) in a difference-in-differences design, in
probability points, using the powerape package for Stata (version 0.2.2 (Stata
port of powerape 1.12.2)), following the confidence-interval approach of
Riesthuis (2024). The assumed data-generating process was a
difference-in-differences design with a treated group (50% of the units) and a
control group, observed in one period before and one period after one common
treatment date as a panel of the same units, with a latent within-unit
correlation of 0.50 shared by all periods (Gaussian copula); outcome rates
before and after treatment of 0.200 and 0.300 in the control group and 0.400
and 0.600 in the treated group (0.533 without treatment); untreated trends
parallel on the logit index; no covariates. The assumed true ATT was 0.067
(6.7 percentage points). At a sample size of n = 3000 units observed in 2
periods (6,000 unit-period observations), simulated power for the
minimum-effect claim (the 90% confidence interval's lower bound exceeding the
smallest effect size of interest, 0.020; a one-sided test at alpha = 5%) was
0.673 (Monte Carlo SE 0.027; 300 replications). Across replications, the
probability of concluding a meaningful effect was 0.673, of detection without
meaningfulness 0.153, of an inconclusive result 0.173, and of equivalence
0.000. Each simulated study was analysed by a pooled logit regression of the
outcome on group, period indicators and the treatment indicator, with the
effect on the treated computed as the average partial effect of the treatment
indicator over the treated post-treatment observations (Wooldridge, 2023), and
delta-method Wald confidence intervals from unit-clustered standard errors;
replications without a maximum-likelihood estimate, from a failed fit or from
separation (a group-by-before/after cell with no events or no non-events),
counted against the claim (0.0% of replications).
```

## Limits

One treated group and one treatment date; treatment assigned to units,
not to clusters of units; balanced panels; covariates that do not change
over time and are not interacted with treatment. Staggered adoption,
effects that differ by period and pre-trend tests are not available. The
Stata command covers these designs from its version 0.2.0
(`powerape dgp NAME, did`); on the package’s website this vignette shows
the Stata code next to the R code.

## References

Ai, C., & Norton, E. C. (2003). Interaction terms in logit and probit
models. *Economics Letters, 80*(1), 123-129.

Caron, L. (2025). Triple difference designs with heterogeneous treatment
effects. arXiv:2502.19620.

Puhani, P. A. (2012). The treatment effect, the cross difference, and
the interaction term in nonlinear “difference-in-differences” models.
*Economics Letters, 115*(1), 85-87.

Roth, J., & Sant’Anna, P. H. C. (2023). When is parallel trends
sensitive to functional form? *Econometrica, 91*(2), 737-747.

Wooldridge, J. M. (2023). Simple approaches to nonlinear
difference-in-differences with panel data. *The Econometrics Journal,
26*(3), C31-C66.
