# Power curve over a grid of sample sizes

Runs the power simulation at each element of `n` and collects the
results. Works for APE and AIE designs alike.

## Usage

``` r
ape_curve(
  dgp,
  n,
  claim = c("minimum", "detect", "equivalence"),
  sesoi = NULL,
  alpha = 0.05,
  conf = NULL,
  nsim = 1000,
  seed = NULL,
  se = c("model", "robust"),
  separation = c("fail", "keep")
)

# S3 method for class 'powerape_curve'
plot(x, target_power = NULL, ...)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
  or
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md).

- n:

  Integer vector (length \>= 2) of total sample sizes.

- claim:

  `"minimum"` (default), `"detect"`, or `"equivalence"`.

- sesoi:

  Smallest effect size of interest, in APE units. Required for
  `"minimum"` and `"equivalence"`; optional for `"detect"` (if supplied,
  the outcome table is still broken out against it).

- alpha:

  The claim's error rate (default 0.05): two-sided for `"detect"`,
  one-sided for `"minimum"`, TOST for `"equivalence"`. May also be a
  **sample-size rule**: a function of `n` returning the error rate to
  use at that `n`, for example `function(n) alphaN::alphaN(n, BF = 3)`,
  the Bayes-factor calibration of Wulff and Taylor (2024) that lowers
  alpha as `n` grows. The searching functions evaluate the rule at every
  candidate `n`, so
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  designs jointly over the pair (alpha(n), n); results store the rule
  and the level it produced.

- conf:

  Optional override: the interval level used for the claim itself
  (`1 - alpha` for detection, `1 - 2 alpha` otherwise). Supplying
  `conf = 0.95` for a minimum-effect or equivalence claim reproduces the
  pre-1.8.0 behavior (one-sided error rate 2.5%).

- nsim:

  Replications per grid point.

- seed:

  Optional; grid point i uses `seed + i - 1`.

- se:

  Standard errors for the exogenous cross-sectional routes: `"model"`
  (default, expected-information ML) or `"robust"`
  (heteroskedasticity-robust HC0 sandwich, as in the sandwich package).
  Panel designs always use unit-clustered SEs and IV designs the stacked
  method-of-moments robust sandwich; `se` is ignored there.

- separation:

  How simulated studies with separation in an identifying cell are
  scored: `"fail"` (default; counted as failed, as the field's reference
  analysis refuses them) or `"keep"` (retained with their degenerate
  Wald intervals). See the section 'Separation and sparse cells'.

- x:

  A `powerape_curve` object.

- target_power:

  Optional horizontal reference line; when the curve crosses it, the
  crossing n is marked and annotated.

- ...:

  Passed to the base plot call.

## Value

A `powerape_curve` object with a `results` data frame (`n`, `alpha`,
`power`, `mcse`, `failed`, `separated`); `alpha` varies along the grid
when it is a sample-size rule, and `separated` is the share of simulated
studies with separation (counted within `failed` under the default
`separation = "fail"`).

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
cv <- ape_curve(d, n = c(300, 600, 900), claim = "detect", nsim = 300, seed = 1)
plot(cv, target_power = 0.8)

# }
```
