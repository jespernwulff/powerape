# Specify a DGP with covariates resampled from pilot data

Instead of parametric marginals and a copula, the covariate distribution
is the empirical joint distribution of a pilot (or register) dataset:
simulation resamples its rows with replacement, and all calibration
integrals are exact averages over the rows. This preserves the covariate
dependence structure for free and is the recommended route when data
exist (DESIGN.md section 3.1).

## Usage

``` r
ape_dgp_empirical(
  model = c("probit", "logit"),
  data,
  focal,
  baseline,
  signal = 0,
  gamma = NULL,
  seed_int = 20260812L
)
```

## Arguments

- model:

  `"probit"` (default) or `"logit"`.

- data:

  Data frame or matrix of numeric columns, no missing values.
  Dummy-encode factors first.

- focal:

  Column name in `data`, or a
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md).

- baseline:

  Baseline outcome rate in (0, 1): the average `P(Y = 1)` with the focal
  at its reference value (0 for binary, the mean for continuous),
  covariates at their empirical distribution.

- signal:

  Nuisance-signal strength (latent pseudo-R-squared of the nuisance
  index), as in
  [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md);
  equal weights on the standardized scale, scaled on the empirical rows.
  Ignored when `gamma` is supplied.

- gamma:

  Optional explicit nuisance coefficients, one per covariate column.

- seed_int:

  Seed for the deterministic draw used only when a continuous `pa_var`
  focal is combined with empirical covariates.

## Value

An object of class `powerape_dgp` (route `"empirical"`).

## Details

The focal variable comes in two flavors:

- a **column name** in `data`: the focal variable is resampled jointly
  with the covariates, preserving its empirical dependence with them
  (observational mental model). Type is inferred (binary if the column
  is 0/1, continuous otherwise; a continuous focal's reference value is
  its mean).

- a
  **[`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)**:
  the focal variable is drawn independently of the resampled rows
  (randomized-treatment mental model); all of `data` is treated as
  covariates.

## Examples

``` r
pilot <- data.frame(treat = rbinom(600, 1, 0.5),
                    age = rnorm(600, 45, 12))
d <- ape_dgp_empirical(data = pilot, focal = "treat",
                       baseline = 0.30, signal = 0.15)
d
#> powerape DGP -- probit (empirical covariates)
#>   focal: treat (binary, p = 0.53) [resampled jointly with covariates]
#>   covariates: age 
#>   pilot rows: 600 (resampled with replacement)
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.150
#>   target effect not set yet - call set_ape() or set_aie().
```
