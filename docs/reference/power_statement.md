# Render a power analysis as a citable methods paragraph

Turns a
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
or
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
result into a self-contained methods paragraph stating the estimand, the
full data-generating assumptions, the claim and CI convention, and the
Monte Carlo precision – the transparency artifact for grant applications
and preregistrations.

## Usage

``` r
power_statement(x)
```

## Arguments

- x:

  A `powerape_power` or `powerape_n` object.

## Value

A character string of class `powerape_statement` (printed wrapped).

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
pw <- ape_power(d, n = 800, claim = "minimum", sesoi = 0.03,
                nsim = 500, seed = 1)
power_statement(pw)
#> We conducted a simulation-based power analysis for the average partial effect
#> (APE) of treat using the powerape package (version 1.3.0), following the
#> confidence-interval approach of Riesthuis (2024). The assumed data-generating
#> process was a probit model with focal variable treat (binary, prevalence
#> 0.50); no additional covariates; baseline outcome rate 0.300 with the focal
#> at reference; nuisance covariates contribute a latent pseudo-R-squared of
#> 0.00. The assumed true APE was 0.100 (10.0 percentage points). At a sample
#> size of n = 800, simulated power for the minimum-effect claim (the 95%
#> confidence interval's lower bound exceeding the smallest effect size of
#> interest, 0.030) was 0.554 (Monte Carlo SE 0.022; 500 replications). Across
#> replications, the probability of concluding a meaningful effect was 0.554, of
#> detection without meaningfulness 0.300, of an inconclusive result 0.146, and
#> of equivalence 0.000. Estimation used maximum likelihood with delta-method
#> Wald confidence intervals; replications that failed to converge (0.0%)
#> counted against the claim.
# }
```
