# IV designs: endogenous focal variables via control functions

## Why an IV module?

Many focal variables of interest are not assigned: training take-up,
board structure, R&D intensity. When the focal variable is endogenous –
correlated with the structural error – a plain probit APE is biased, and
the modern flexible answer is the **control-function (CF) approach**
(Wooldridge, 2015): model the endogenous variable in a first stage on
instruments, and carry its (generalized) residual into the outcome
probit as an extra regressor. Stata 18.5 ships this as `cfprobit`, with
a clever stacked method-of-moments variance that avoids bootstrapping
entirely.

`powerape` replicates that estimator exactly – point estimates,
`margins`-style ASF effects, and the stacked robust and cluster-robust
standard errors, verified against Stata to about 1e-7 – and then wraps
the usual power machinery around it. The estimand is unchanged: the
ASF-based APE (or AIE) of the focal variable, which is what `margins`
reports after `cfprobit`.

## Step 1: describe the world

Two knobs carry the IV structure, and both are quantities you can
defend:

- `endogeneity` – the correlation between the structural and first-stage
  errors. This is *why* you need instruments; 0 means plain probit would
  have been fine.
- `iv_strength` – the share of the focal’s variance (latent variance for
  a binary focal) the instruments explain. This is the first-stage
  strength you would defend with an F statistic in the paper.

``` r

library(powerape)

d <- ape_dgp_iv(
  focal       = pa_var("training", "binary", p = 0.4),   # probit first stage
  covariates  = list(pa_var("size", "normal")),
  instruments = pa_var("subsidy", "binary", p = 0.5),
  endogeneity = 0.4,
  iv_strength = 0.2,
  baseline    = 0.30,
  signal      = 0.10
)
d <- set_ape(d, target = 0.10)
d
#> powerape DGP -- probit
#>   focal: training (binary, p = 0.40)
#>   covariates: size (normal, mean 0.00, sd 1.00) 
#>   iv: endogenous focal (probit first stage; instruments: subsidy); endogeneity rho = 0.40, first-stage strength = 0.20; control-function probit, stacked robust SEs
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.100
#>   true APE pinned at +0.1000 (beta_focal = 0.2857, implied P(Y=1 | focal=1) = 0.400)
```

A continuous focal (`pa_var("rd", "normal", ...)`) gets a linear first
stage instead; the pinning step is identical.

## Step 2: power – and the IV price

``` r

ape_power(d, n = 1200, claim = "detect", nsim = 400, seed = 1)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 1200, assumed true APE +0.1000, 95% CI, nsim = 400, stacked robust SEs
#>   power = 0.210 (MCSE 0.020)
#>   outcomes: detect 0.210 | inconclusive 0.790 | failed 0.000
```

Compare that with the same world *without* endogeneity, estimated by
plain probit:

``` r

ds <- ape_dgp(
  focal      = pa_var("training", "binary", p = 0.4),
  covariates = list(pa_var("size", "normal")),
  baseline   = 0.30,
  signal     = 0.10
)
ds <- set_ape(ds, target = 0.10)
ape_power(ds, n = 1200, claim = "detect", nsim = 400, seed = 2)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 1200, assumed true APE +0.1000, 95% CI, nsim = 400
#>   power = 0.955 (MCSE 0.010)
#>   outcomes: detect 0.955 | inconclusive 0.045 | failed 0.000
```

The gap is **the IV price**: the CF estimator identifies the effect from
the instrumented variation only, so its standard error is larger by
roughly `sqrt(1 / iv_strength)` – the familiar 2SLS variance logic.
Studies that budget their sample size for an exogenous analysis and then
switch to IV at review time are, mechanically, badly underpowered. This
pair of calls prices that in advance.

Since `endogeneity` and `iv_strength` are guesses, sweep them:

``` r

ape_robust(d, n = 1200, claim = "detect",
           vary = list(iv_strength = c(0.1, 0.3)), grid_points = 3,
           nsim = 300, seed = 3, nmax = FALSE)
#> powerape robustness sweep -- detect claim, APE, pin = ape
#>   n = 1200, nsim = 300 per scenario, 3 scenario(s) over: iv_strength
#>   power range [0.133, 0.353]; worst scenario:
#>  iv_strength implied_effect     power       mcse
#>          0.1            0.1 0.1333333 0.01962614
#>   scenarios:
#>  iv_strength implied_effect     power       mcse
#>          0.1            0.1 0.1333333 0.01962614
#>          0.2            0.1 0.2200000 0.02391652
#>          0.3            0.1 0.3533333 0.02759764
```

## Endogenous interactions

Moderation with an endogenous focal is the fully general case: does a
(self-selected) training program work differently in public firms? The
estimating model gains the focal-by-moderator term, the control
function, and the control-function-by-moderator interaction (Stata:
`interact()` + `mainonly()`), and the first stage includes the moderator
**and the instrument-by-moderator interactions** – as instruments must
in interaction models.

``` r

di <- ape_dgp_iv(
  focal       = pa_var("training", "binary", p = 0.4),
  moderator   = pa_var("public", "binary", p = 0.5),
  covariates  = list(pa_var("size", "normal")),
  instruments = pa_var("subsidy", "binary", p = 0.5),
  endogeneity = 0.4,
  iv_strength = 0.25,
  baseline    = 0.30,
  signal      = 0.10
)
di <- set_aie(di, target = 0.06, main_focal = 0.10, main_moderator = 0.05)
ape_power(di, n = 2000, claim = "detect", nsim = 250, seed = 4)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 2000, assumed true AIE +0.0600, 95% CI, nsim = 250, stacked robust SEs
#>   power = 0.064 (MCSE 0.015)
#>   outcomes: detect 0.064 | inconclusive 0.936 | failed 0.000
```

The interaction tax and the IV price compound. Interaction claims from
instrumented designs need very large samples, and it is far better to
learn that at the design stage than from a null result.

## What stands behind the numbers

On fixed datasets, powerape’s CF estimator reproduces Stata’s `cfprobit`
across first stages (linear, probit, fractional probit, Poisson):
coefficients and stacked SEs to ~1e-7 (robust and clustered), and the
`margins` ASF effects exactly, including the interaction design run with
Stata’s own syntax. One caveat we discovered in the process: for
*fractional* first stages, Stata’s `margins` rebuilds the control
function with a binary-outcome formula that ignores the fractional value
– an internal inconsistency with its own estimator, reproduced to 7
decimals in our validation battery. powerape uses the
estimation-consistent generalized residuals.

Current scope: one endogenous focal variable; linear and probit first
stages in the DGP route (fractional probit and Poisson are supported at
the estimator level); binary exogenous moderators for AIE designs;
robust (stacked sandwich) inference.

## References

- Rivers, D., & Vuong, Q. H. (1988). Limited information estimators and
  exogeneity tests for simultaneous probit models. *Journal of
  Econometrics*, 39(3), 347-366.
- Wooldridge, J. M. (2015). Control function methods in applied
  econometrics. *Journal of Human Resources*, 50(2), 420-445.
