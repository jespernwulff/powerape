# Choose the analysis the simulated studies receive

By default every simulated study is analysed with the design's index
model: a probit or logit fit, with the average partial effect (or
average interaction effect) and its delta-method standard error.
`set_analysis(dgp, "linear")` returns a copy of the design that is
analysed with a **linear probability model** instead: ordinary least
squares on the same regressors, where the estimate is the coefficient of
the focal variable (APE designs), of the focal-by-moderator product (AIE
designs), of the triple product (three-way designs), or of the treatment
indicator (difference-in-differences designs). One design can so be
priced both ways, for comparison.

## Usage

``` r
set_analysis(dgp, analysis = c("index", "linear"))
```

## Arguments

- dgp:

  A `powerape_dgp`.

- analysis:

  `"index"` (the design's probit or logit model, the default of every
  design) or `"linear"` (linear probability model).

## Value

A copy of `dgp` with the analysis set.

## Details

The world and the pinned effect do not change: the truth remains the
average partial effect (or AIE, or ATT) set with
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md),
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md),
or
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md).
The linear coefficient is another estimator of that effect and need not
centre on it. In large samples it equals the average partial effect when
a binary focal variable is independent of the covariates (as under
randomization) and when the focal variable and the covariates are
jointly normal; it equals the average interaction effect when a binary
focal variable and a binary moderator are independent of the covariates.
Otherwise it differs, and it can differ substantially when a binary
focal variable is related to covariates that also move the outcome: the
linear control then leaves part of the nonlinear confounding in the
coefficient. (With a 30% treated group, a 15% baseline, one covariate
with `signal = 0.3` and a latent correlation of 0.3 between the two, the
linear coefficient centres about 12% above the average partial effect.)
The printed design shows the value the linear analysis centres on in
large samples next to the truth, and
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
reports the mean of the simulated estimates.

## Standard errors and separation

A linear probability model is heteroskedastic by construction, so its
standard errors are always robust and the `se` argument of the power
functions is ignored: HC1 for independent observations (as Stata's
`regress, vce(robust)`) and clustered by unit for panel designs (with
the small-sample factor of Stata's `regress, vce(cluster)`). Ordinary
least squares has no separation, so no simulated study is discarded for
a cell without events or without non-events; the sparse-cell warning
still applies.

## Limits

Instrumental-variable designs are refused: their linear counterpart is
two-stage least squares, which is not available.
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
describes a linear analysis for difference-in-differences designs only.
Before searching,
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
checks that the value the linear analysis centres on lies on the same
side of the claim's boundary as the true effect, and stops otherwise.

## Examples

``` r
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
             covariates = list(pa_var("age", "normal")),
             baseline = 0.30, signal = 0.15)
d <- set_ape(d, target = 0.10)
dl <- set_analysis(d, "linear")
dl
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50)
#>   covariates: age (normal, mean 0.00, sd 1.00) 
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.150
#>   true APE pinned at +0.1000 (beta_focal = 0.2940, implied P(Y=1 | focal=1) = 0.400)
#>   analysis: linear probability model (OLS); estimate = coefficient of treat; robust (HC1) SEs
#>   in large samples that coefficient is +0.1000 (the true APE is +0.1000)
# \donttest{
ape_power(d, n = 800, claim = "detect", nsim = 300, seed = 1)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 800, assumed true APE +0.1000, 95% CI (alpha = 0.05), nsim = 300
#>   power = 0.867 (MCSE 0.020)
#>   outcomes: detect 0.867 | inconclusive 0.133 | failed 0.000 
ape_power(dl, n = 800, claim = "detect", nsim = 300, seed = 1)
#> powerape -- detection claim (CI excludes 0, directional)
#>   linear probability (OLS), n = 800, assumed true APE +0.1000, 95% CI (alpha = 0.05), nsim = 300, robust (HC1) SEs
#>   power = 0.863 (MCSE 0.020)
#>   the linear estimates average +0.1005 (MCSE 0.0019); the true APE is +0.1000
#>   outcomes: detect 0.863 | inconclusive 0.137 | failed 0.000 
# }
```
