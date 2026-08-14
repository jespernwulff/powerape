# Panel designs: correlated random effects

## Why a panel module?

Much applied binary-outcome work is panel work: firm-years,
person-waves, plant-quarters. The workhorse model for that setting is a
**correlated random effects (CRE) probit or logit** – pooled estimation
with Mundlak means and unit-clustered standard errors (Wooldridge, 2010,
ch. 15). It nests the pure random-effects model, allows the unit
heterogeneity to correlate with the regressors through their observed
unit means, and is what `xtprobit`-style advice converges to when the
heterogeneity cannot credibly be assumed independent.

The good news for design is that **the estimand does not change**. The
target is still the APE – here the average-structural-function (ASF)
version: the effect on P(Y = 1) integrating over the unit effects, with
the Mundlak means held at their observed values. The pooled CRE
estimator recovers exactly that quantity even though its *coefficients*
are attenuated. Designing in APE units is right for panels too; that is
the point of the package.

What changes is the **accounting**: for panel DGPs, every sample-size
argument
([`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md))
counts **units (clusters)**, not observations – clusters are what drive
power.

## Step 1: describe the panel world

``` r

library(powerape)

d <- ape_dgp_panel(
  model      = "probit",
  focal      = pa_var("treat", "binary", p = 0.5, icc = 1),  # unit-level treatment
  covariates = list(pa_var("size", "normal", icc = 0.6)),
  n_periods  = 4,
  rho        = 0.30,
  cre_share  = 0.50,
  baseline   = 0.30,
  signal     = 0.10
)
d
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50, icc = 1.00)
#>   covariates: size (normal, mean 0.00, sd 1.00, icc = 0.60) 
#>   panel: 4 periods per unit; rho = 0.30 (latent unit share), cre_share = 0.50; CRE probit w/ Mundlak means, unit-clustered SEs
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.100
#>   target effect not set yet - call set_ape() or set_aie().
```

Three knobs are panel-specific, and each is a quantity you can defend in
a grant application:

- `n_periods` – balanced waves per unit.
- `rho` – the share of latent variance coming from the unit effect,
  `Var(c) / (Var(c) + 1)`; the same convention `xtprobit` reports. Zero
  means no unobserved heterogeneity.
- `cre_share` – how much of that unit effect is *correlated* with the
  regressors (loaded on the centered unit means of the time-varying
  variables); the remainder is a pure random effect. The estimator
  includes Mundlak means either way.

Per-variable persistence comes from `pa_var(..., icc = )`: `icc = 1` is
a time-constant, unit-level variable (a one-shot treatment, a sector);
values in (0, 1) mix unit-level and wave-level variation.

## Step 2: pin the APE

``` r

d <- set_ape(d, target = 0.10)
```

Identical to the cross-sectional workflow – the inversion solves the
focal coefficient so the true ASF-APE is exactly 0.10, integrating the
unit effects (in closed form for probit, by Gauss-Hermite quadrature for
logit), and errors informatively if the target is infeasible.

## Step 3: power, in units

``` r

ape_power(d, n = 300, claim = "minimum", sesoi = 0.05,
          nsim = 400, seed = 1)
#> powerape -- minimum-effect claim (CI lower bound > 0.050)
#>   probit, n = 300 units x 4 periods (1200 obs), assumed true APE +0.1000, 95% CI, nsim = 400
#>   power = 0.442 (MCSE 0.025)
#>   outcomes: minimum 0.442 | detect-only 0.505 | inconclusive 0.052 | equivalence 0.000 | failed 0.000
```

`n = 300` means 300 units observed for 4 periods – 1,200 unit-period
observations, fitted by pooled probit with Mundlak means and
unit-clustered SEs in every replication. (Below 30 clusters the package
warns: cluster-robust inference is not trustworthy there.)

## Panels are not free observations

With a time-constant treatment, T waves of the same units add much less
than T-fold information: the outcome is correlated within units, so the
effective sample is closer to the number of clusters than to the number
of rows. Compare the panel above with a cross-section of the same
*total* size:

``` r

cs <- ape_dgp(
  model      = "probit",
  focal      = pa_var("treat", "binary", p = 0.5),
  covariates = list(pa_var("size", "normal")),
  baseline   = 0.30,
  signal     = 0.10
)
cs <- set_ape(cs, target = 0.10)
ape_power(cs, n = 1200, claim = "minimum", sesoi = 0.05,
          nsim = 400, seed = 2)
#> powerape -- minimum-effect claim (CI lower bound > 0.050)
#>   probit, n = 1200, assumed true APE +0.1000, 95% CI, nsim = 400
#>   power = 0.465 (MCSE 0.025)
#>   outcomes: minimum 0.465 | detect-only 0.502 | inconclusive 0.032 | equivalence 0.000 | failed 0.000
```

How the two compare is a *computable property of the design*, not a rule
of thumb. Here — with half of the firm heterogeneity loading on
observable means that the Mundlak terms absorb — the panel’s rows turn
out to be worth nearly as much as independent ones (high-replication
runs put the two powers within a point of each other). Push `rho` up, or
set `cre_share` toward zero so the heterogeneity is unobservable, and
the panel discount grows sharply. Reviewers who either count panel rows
as independent observations *or* apply a crude cluster-trial discount
are both guessing; this pair of calls replaces the guess.

Since `rho` is usually a guess, sweep it – with `pin = "ape"` the true
effect is re-pinned in every scenario, so only the precision channel
moves:

``` r

ape_robust(d, n = 300, claim = "minimum", sesoi = 0.05,
           vary = list(rho = c(0.1, 0.5)), grid_points = 3,
           nsim = 300, seed = 3, nmax = FALSE)
#> powerape robustness sweep -- minimum claim, APE, pin = ape
#>   n = 300, nsim = 300 per scenario, 3 scenario(s) over: rho
#>   power range [0.390, 0.477]; worst scenario:
#>  rho implied_effect power       mcse
#>  0.5            0.1  0.39 0.02816026
#>   scenarios:
#>  rho implied_effect     power       mcse
#>  0.1            0.1 0.4766667 0.02883606
#>  0.3            0.1 0.4300000 0.02858321
#>  0.5            0.1 0.3900000 0.02816026
```

## Panel interactions: the AIE under CRE

Moderation questions survive the move to panels: does a governance
reform bite harder in public firms? The estimand is the **average
interaction effect (AIE)** – the double difference of counterfactual
cell probabilities – and for panel DGPs it is defined on the same ASF as
the APE: unit effects integrated out, Mundlak means held fixed. It is
the same population quantity a cross-sectional AIE design targets.

``` r

di <- ape_dgp_panel(
  model      = "probit",
  focal      = pa_var("treat", "binary", p = 0.5, icc = 1),
  moderator  = pa_var("public", "binary", p = 0.4, icc = 0.7),
  covariates = list(pa_var("size", "normal", icc = 0.6)),
  n_periods  = 4,
  rho        = 0.30,
  cre_share  = 0.50,
  baseline   = 0.30,
  signal     = 0.10
)
di <- set_aie(di, target = 0.06, main_focal = 0.10, main_moderator = 0.05)
di
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50, icc = 1.00)
#>   moderator: public (binary, p = 0.40, icc = 0.70) (focal-by-moderator interaction included)
#>   covariates: size (normal, mean 0.00, sd 1.00, icc = 0.60) 
#>   panel: 4 periods per unit; rho = 0.30 (latent unit share), cre_share = 0.50; CRE probit w/ Mundlak means, unit-clustered SEs
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.100
#>   true AIE pinned at +0.0600 (beta_int = 0.1830)
#>   main APEs at reference: treat +0.100 (beta = 0.3564), public +0.050 (beta = 0.1829)
```

The three anchors are the familiar ones: the two
conditional-at-reference main-effect APEs plus the target AIE. One
panel-specific detail is handled for you: whenever the interaction term
varies within units, the estimating model includes the interaction’s
*own* Mundlak mean (the mean of a product is not the product of the
means).

Interactions are expensive, and panels make them more so – the
interaction tax and the clustering penalty compound:

``` r

ape_power(di, n = 400, claim = "detect", nsim = 250, seed = 4)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 400 units x 4 periods (1600 obs), assumed true AIE +0.0600, 95% CI, nsim = 250
#>   power = 0.112 (MCSE 0.020)
#>   outcomes: detect 0.112 | inconclusive 0.888 | failed 0.000
```

Even 400 units (1,600 observations) leave detection power low for a
6-point AIE. Budget accordingly –
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
with a larger `nsim` gives the required number of units, and
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
renders any of these results as a citable methods paragraph that counts
units and periods correctly.

## What stands behind the numbers

The panel machinery is validated the same way as the rest of the
package: the ASF-based APE and the AIE double difference match Stata
(`probit ..., vce(cluster id)` + `margins`) to all printed digits on
fixed panels, with the small clustered-SE difference fully attributed to
Stata’s observed-Hessian bread; the clustered sandwich matches
[`sandwich::vcovCL`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
exactly; the AIE reproduces `ginteff` on the identical fit; and the
engine reduces to exact enumeration and to the Donner-Klar cluster
design effect in the cases where those are available. See the validation
battery in the repository for the full dossier.

Current scope: balanced panels, the parametric covariate route, and
binary-by-binary panel interactions; continuous panel pairs, unbalanced
panels, and time effects are on the roadmap.

## References

- Mundlak, Y. (1978). On the pooling of time series and cross section
  data. *Econometrica*, 46(1), 69-85.
- Riesthuis, P. (2024). Simulation-based power analyses for the smallest
  effect size of interest: A confidence-interval approach for
  minimum-effect and equivalence testing. *AMPPS*, 7(2).
- Wooldridge, J. M. (2010). *Econometric Analysis of Cross Section and
  Panel Data* (2nd ed.). MIT Press. Ch. 15.
