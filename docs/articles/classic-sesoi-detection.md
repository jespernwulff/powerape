# The classic SESOI power analysis, and what it cannot tell you

## The classic recipe

The most common SESOI power analysis powers a *null-hypothesis test* at
a true effect equal to the SESOI: “if the effect is exactly as small as
the smallest thing I care about, I still detect it with 80%
probability”. In `powerape` this is `claim = "detect"` with the target
set to the SESOI:

``` r

library(powerape)

d <- ape_dgp(model = "probit",
             focal = pa_var("treat", "binary", p = 0.5),
             baseline = 0.30)
d <- set_ape(d, target = 0.05)   # truth = SESOI = 5 pp

ape_power(d, n = 5000, claim = "detect", sesoi = 0.05,
          nsim = 600, seed = 1)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 5000, assumed true APE +0.0500, 95% CI, nsim = 600
#>   power = 0.962 (MCSE 0.008)
#>   outcomes: minimum 0.030 | detect-only 0.932 | inconclusive 0.008 | equivalence 0.030 | failed 0.000
```

Note what the outcome distribution shows: even with ample detection
power, essentially **no study can claim the effect is meaningfully
large** (`minimum` outcome near zero), because when the truth sits at
the SESOI the CI lower bound almost never clears it. Detection is all
the classic recipe buys.

## The boundary is a wall, not a target

Asking for minimum-effect power *at* the SESOI is statistically
degenerate, and `powerape` refuses with an explanation rather than
returning a useless number:

``` r

ape_power(d, n = 5000, claim = "minimum", sesoi = 0.05)
#> Error:
#> ! Assumed true APE (0.0500 in magnitude) does not exceed the SESOI (0.0500). On or below the null boundary of the minimum-effect test the claim succeeds with probability ~2.5% (the test size) at ANY sample size, so no n delivers meaningful power. Either assume a planning value above the SESOI, or run the classic SESOI analysis instead: claim = "detect" with the target set to the SESOI.
```

On the null boundary of the minimum-effect test, the claim succeeds with
probability equal to the test size (about 2.5% with a 95% CI) at *any*
sample size.

## The two coherent moves

Either stay with detection at truth = SESOI (this vignette), or assume a
planning value **above** the SESOI and power the minimum-effect claim.
The second buys a stronger conclusion and costs more observations:

``` r

## detection power for truth = SESOI = 0.05
ape_n(d, power = 0.80, claim = "detect", nsim = 500, seed = 2)
#> powerape required sample size -- detect claim
#>   n = 2742 for 80% target power (confirmed 0.799, MCSE 0.009)
#>   assumed true APE +0.0500, sesoi -, 95% CI, probit
#>   search: 1 step(s); confirmed in 1 round(s) at nsim = 2000.

## minimum-effect power for truth = 0.08 against SESOI = 0.05
d8 <- set_ape(d, target = 0.08)
ape_n(d8, power = 0.80, claim = "minimum", sesoi = 0.05,
      nsim = 500, seed = 3)
#> powerape required sample size -- minimum claim
#>   n = 8054 for 80% target power (confirmed 0.826, MCSE 0.008)
#>   assumed true APE +0.0800, sesoi 0.050, 95% CI, probit
#>   search: 1 step(s); confirmed in 2 round(s) at nsim = 2000.
```

The gap between those two n’s is the price of concluding “meaningfully
large” instead of “not zero”. Making that price visible – in APE units,
before data collection – is precisely what a SESOI-based design
discussion needs.

## Size audit

A useful sanity check of the whole simulation pipeline: at a true APE of
zero, directional detection should fire at half the nominal two-sided
rate.

``` r

d0 <- set_ape(d, target = 0)
ape_power(d0, n = 2000, claim = "detect", nsim = 600, seed = 4)
#> powerape -- detection claim (CI excludes 0, directional)
#>   probit, n = 2000, assumed true APE +0.0000, 95% CI, nsim = 600
#>   power = 0.022 (MCSE 0.006)
#>   outcomes: detect 0.022 | inconclusive 0.978 | failed 0.000
```
