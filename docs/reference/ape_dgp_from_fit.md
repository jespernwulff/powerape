# Specify a DGP from a fitted pilot model

The third covariate route (DESIGN.md section 3.1): everything except the
focal coefficient is lifted from a fitted pilot `glm` — the covariate
rows come from the pilot's model matrix (resampled with replacement, so
their joint distribution is preserved), the nuisance coefficients are
the pilot estimates, and the intercept is the pilot's unless `baseline`
overrides it. The focal coefficient itself is *not* taken from the
pilot: it is re-solved in APE units by
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md),
keeping the "power for the effect size you care about" logic intact.

## Usage

``` r
ape_dgp_from_fit(fit, focal, baseline = NULL)
```

## Arguments

- fit:

  A fitted `glm` with `family = binomial("probit")` or
  `binomial("logit")`.

- focal:

  Name of the focal regressor's model-matrix column: a numeric variable
  in the model (binary 0/1 or continuous), or the single dummy column of
  a two-level factor or logical (e.g. `"treatTRUE"`).

- baseline:

  Optional override, in (0, 1): recalibrate the intercept so the average
  `P(Y = 1)` with the focal at reference equals this value. Default
  `NULL` keeps the pilot-implied baseline (reported in the object).

## Value

An object of class `powerape_dgp` (route `"empirical"`, builder
`"from_fit"`). Pin the effect with
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md),
then use the usual power functions;
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
can vary `baseline`.

## Details

The pilot fit must be a binomial `glm` with a probit or logit link, no
offset, and a main-effects role for the focal variable: the focal must
enter the model through exactly one column of the model matrix, and its
variable may appear in no other term – no interactions with it, no
polynomial or other transformed copies of it (`I(dose^2)` next to
`dose`, `log(dose)`, splines), because those columns would be held at
their pilot values while the focal changes, so the pinned effect would
not be the APE that `margins` or marginaleffects report for the same
model. Such pilots are refused with an error (use
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md) +
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
for interaction designs). Nuisance-side interactions and transformed
covariates are fine: they are just columns of the model matrix. A
precomputed column that is a function of the focal (a data column
`dose_sq`) cannot be detected; leave such columns out.

## Examples

``` r
set.seed(1)
pilot <- data.frame(treat = rbinom(600, 1, 0.5),
                    age = rnorm(600, 45, 12))
pilot$y <- rbinom(600, 1, pnorm(-0.6 + 0.3 * pilot$treat + 0.01 * pilot$age))
fit <- glm(y ~ treat + age, binomial("probit"), pilot)
d <- ape_dgp_from_fit(fit, focal = "treat")
d <- set_ape(d, target = 0.10)
d
#> powerape DGP -- probit (empirical covariates)
#>   focal: treat (binary, p = 0.48) [resampled jointly with covariates]
#>   covariates: age 
#>   pilot rows: 600 (resampled with replacement)
#>   pilot model: y ~ treat + age (nuisance coefficients lifted; focal coefficient re-solved)
#>   baseline P(Y=1 | reference) = 0.404, nuisance signal (latent pseudo-R2) = 0.032
#>   true APE pinned at +0.1000 (beta_focal = 0.2571, implied P(Y=1 | focal=1) = 0.504)
```
