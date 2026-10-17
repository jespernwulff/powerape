# Power for average interaction effects (AIEs)

## The estimand

“Does the treatment work differently for women?” is a question about the
**average interaction effect (AIE)** – in a nonlinear model *not* the
interaction coefficient (Ai & Norton, 2003), but the averaged double
difference / cross-partial that the
[`ginteff`](https://github.com/jespernwulff/ginteff) package estimates.
`powerape` reimplements the AIE internally for speed, unit-tested to
reproduce `ginteff`’s estimates and delta-method standard errors
exactly.

## Specifying and pinning an interaction design

Add a `moderator` to the DGP; the index then includes both mains and the
interaction. Pinning happens in effect-size units via
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md),
which takes three explicit anchors:

- `main_focal`, `main_moderator` – the *conditional-at-reference APEs*:
  the APE of each variable with the other held at its reference value (0
  for binary variables);
- `target` – the AIE itself.

For a binary-by-binary design this is exactly a bijection with the four
counterfactual cell rates: baseline p00, p10 = p00 + main_focal, p01 =
p00 + main_moderator, p11 = p00 + main_focal + main_moderator + AIE.

- R
- Stata

``` r

library(powerape)

d <- ape_dgp(model = "probit",
             focal = pa_var("treat", "binary", p = 0.5),
             moderator = pa_var("female", "binary", p = 0.55),
             covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
             baseline = 0.30, signal = 0.10)
d <- set_aie(d, target = 0.08, main_focal = 0.10, main_moderator = 0.03)
d
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50)
#>   moderator: female (binary, p = 0.55) (focal-by-moderator interaction included)
#>   covariates: age (normal, mean 45.00, sd 12.00) 
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.100
#>   true AIE pinned at +0.0800 (beta_int = 0.2044)
#>   main APEs at reference: treat +0.100 (beta = 0.2857), female +0.030 (beta = 0.0891)
```

``` stata
powerape dgp d, focal(treat, binary p(0.5)) moderator(female, binary p(0.55)) ///
    covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.10)

powerape setaie d, target(0.08) mainfocal(0.10) mainmoderator(0.03)
```

``` stata-output
powerape DGP -- probit
  focal: treat (binary, p = 0.50)
  moderator: female (binary, p = 0.55) (focal-by-moderator interaction included)
  covariates: age (normal, mean 45.00, sd 12.00)
  baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.100
  true AIE pinned at + 0.0800 (beta_int =    0.2044)
  main APEs at reference: treat  0.100 (beta =    0.2857), female  0.030 (beta =    0.0891)
```

Feasibility is checked at every step (all four cell rates must stay
inside (0, 1); with continuous variables the AIE plateaus in the
interaction coefficient), with errors reporting the attainable range.

## Interactions are expensive

The estimand rides on the DGP, so the usual functions power the AIE
directly:

- R
- Stata

``` r

ape_power(d, n = 3000, claim = "detect", nsim = 500, seed = 1)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 3000, assumed true AIE +0.0800, 95% CI (alpha = 0.05), nsim = 500
#>   power = 0.662 (MCSE 0.021)
#>   outcomes: detect 0.662 | inconclusive 0.338 | failed 0.000
```

``` stata
powerape power d, n(3000) claim(detect) nsim(500) seed(1)
```

``` stata-output
powerape -- detection claim (CI excludes 0, directional)
  probit, n = 3000, assumed true AIE + 0.0800,     95% CI (alpha = .05), nsim = 500
  fit: probit via engine(stata), effects(closed)
  power = 0.640 (MCSE 0.021)
  outcomes: detect 0.640 | inconclusive 0.360 | failed 0.000
```

Compare that with how easy an equally sized *main effect* would be:

- R
- Stata

``` r

d_main <- ape_dgp(model = "probit",
                  focal = pa_var("treat", "binary", p = 0.5),
                  covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
                  baseline = 0.30, signal = 0.10)
d_main <- set_ape(d_main, target = 0.08)
ape_power(d_main, n = 3000, claim = "detect", nsim = 500, seed = 2)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 3000, assumed true APE +0.0800, 95% CI (alpha = 0.05), nsim = 500
#>   power = 0.996 (MCSE 0.003)
#>   outcomes: detect 0.996 | inconclusive 0.004 | failed 0.000
```

``` stata
powerape dgp d_main, focal(treat, binary p(0.5))                          ///
    covariates((age, normal mean(45) sd(12))) baseline(0.30) signal(0.10)
powerape setape d_main, target(0.08)

powerape power d_main, n(3000) claim(detect) nsim(500) seed(2)
```

``` stata-output
powerape -- detection claim (CI excludes 0, directional)
  probit, n = 3000, assumed true APE + 0.0800,     95% CI (alpha = .05), nsim = 500
  fit: probit via engine(stata), effects(closed)
  power = 0.992 (MCSE 0.004)
  outcomes: detect 0.992 | inconclusive 0.008 | failed 0.000
```

Same 8-point effect, same n – a fraction of the power. The folk theorem
that interaction studies need several times the sample is here made
exact, in the units researchers actually interpret, via
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md):

- R
- Stata

``` r

ape_n(d, power = 0.80, claim = "detect", nsim = 500, seed = 3)
#> powerape required sample size -- detect claim
#>   n = 4490 for 80% target power (confirmed 0.821, MCSE 0.009)
#>   assumed true AIE +0.0800, sesoi -, 95% CI (alpha = 0.05), probit
#>   search: 1 step(s); confirmed in 2 round(s) at nsim = 2000.
```

``` stata
powerape n d, power(0.80) claim(detect) nsim(500) seed(3) engine(mata)
```

``` stata-output
powerape required sample size -- detect claim
  n = 4279 for 80% target power (confirmed 0.798, MCSE 0.009)
  assumed true AIE + 0.0800, sesoi -,     95% CI (alpha =    .05), probit
  fit: Mata engine (Fisher scoring)
  search: 1 step(s); confirmed in 1 round(s) at nsim = 2000.
```

All claims work for AIEs: power a *meaningful* interaction with
`claim = "minimum"` and an AIE SESOI, or show treatment-effect
homogeneity with `claim = "equivalence"`.

## Correspondence with ginteff (and Stata)

powerape’s binary AIE is the 0-to-1 counterfactual contrast. In
`ginteff` that corresponds to **factor** variables passed to `dydxs`
(the analog of Stata’s `i.` factor syntax): discrete change from the
base level. Numeric variables in `firstdiff` instead get “+1 unit from
the observed value” – for a 0/1 dummy that extrapolates treated units to
2, a different estimand. When validating a powerape design against
`ginteff` on pilot data, declare binary variables as factors.

In Stata, `effects(margins)` computes the AIE of every simulated study
with Radean’s (2023) `ginteff` command itself instead of the closed-form
delta method. On the same seed the two reach the same decision in every
study; `ginteff` is about ten times slower, so it serves as a check.

- R
- Stata

The R package computes the AIE with its internal reimplementation,
unit-tested against `ginteff` for R.

``` stata
powerape power d, n(3000) claim(detect) nsim(25) seed(1)

powerape power d, n(3000) claim(detect) nsim(25) seed(1) effects(margins)
```

``` stata-output
powerape -- detection claim (CI excludes 0, directional)
  probit, n = 3000, assumed true AIE + 0.0800,     95% CI (alpha = .05), nsim = 25
  fit: probit via engine(stata), effects(closed)
  power = 0.760 (MCSE 0.085)
  outcomes: detect 0.760 | inconclusive 0.240 | failed 0.000

powerape -- detection claim (CI excludes 0, directional)
  probit, n = 3000, assumed true AIE + 0.0800,     95% CI (alpha = .05), nsim = 25
  fit: probit via engine(stata), effects(margins)
  power = 0.760 (MCSE 0.085)
  outcomes: detect 0.760 | inconclusive 0.240 | failed 0.000
```

## References

- Ai, C., & Norton, E. C. (2003). Interaction terms in logit and probit
  models. *Economics Letters*, 80(1), 123-129.
- Radean, M. (2023). ginteff: A generalized command for computing
  interaction effects. *The Stata Journal*, 23(2).
