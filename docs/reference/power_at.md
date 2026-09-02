# Power at other error rates from one stored simulation

[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
keeps the estimate and standard error of every replication
(`keep_draws = TRUE`, the default), so the claim can be re-evaluated at
any error rate without simulating again: power at `alpha` is a
re-threshold of the same draws, exact on those draws and monotone in
`alpha`. This is the curve that
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
optimizes.

## Usage

``` r
power_at(x, alpha)
```

## Arguments

- x:

  A `powerape_power` object from
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
  with stored draws.

- alpha:

  Numeric vector of error rates in (0, 0.5), each read in the claim's
  conventional form (two-sided for detection, one-sided for the
  minimum-effect claim, TOST for equivalence; see
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)).

## Value

A data frame with one row per `alpha`: `alpha`, the interval level
`conf` it implies for the claim, `power`, its Monte Carlo standard error
`mcse`, and the outcome distribution.

## Examples

``` r
# \donttest{
d <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                     baseline = 0.30), 0.10)
pw <- ape_power(d, n = 800, claim = "minimum", sesoi = 0.03,
                nsim = 500, seed = 1)
power_at(pw, alpha = c(0.01, 0.025, 0.05, 0.10))
#>   alpha conf power       mcse minimum detect_only inconclusive equivalence
#> 1 0.010 0.98 0.416 0.02204287   0.416       0.252        0.332           0
#> 2 0.025 0.95 0.554 0.02222989   0.554       0.230        0.216           0
#> 3 0.050 0.90 0.678 0.02089574   0.678       0.176        0.146           0
#> 4 0.100 0.80 0.814 0.01740138   0.814       0.102        0.084           0
#>   failed
#> 1      0
#> 2      0
#> 3      0
#> 4      0
# }
```
