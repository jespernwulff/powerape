# Render a power analysis as a citable methods paragraph

Turns a
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
or
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
result into a self-contained methods paragraph stating the estimand, the
full data-generating assumptions, the claim and CI convention (including
a sample-size rule for `alpha` when one was used), and the Monte Carlo
precision – the transparency artifact for grant applications and
preregistrations. For an
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
result the paragraph states the justification of the error rate instead:
the objective and its weights, the chosen level and its flat range, the
cap, and the realized-size verification when it was run.

## Usage

``` r
power_statement(x)
```

## Arguments

- x:

  A `powerape_power`, `powerape_n`, `powerape_mde`, or `powerape_alpha`
  object.

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
#> (APE) of treat using the powerape package (version 1.11.0), following the
#> confidence-interval approach of Riesthuis (2024). The assumed data-generating
#> process was a probit model with focal variable treat (binary, prevalence
#> 0.50); no additional covariates; baseline outcome rate 0.300 with the focal
#> at reference; nuisance covariates contribute a latent pseudo-R-squared of
#> 0.00. The assumed true APE was 0.100 (10.0 percentage points). At a sample
#> size of n = 800, simulated power for the minimum-effect claim (the 90%
#> confidence interval's lower bound exceeding the smallest effect size of
#> interest, 0.030; a one-sided test at alpha = 5%) was 0.678 (Monte Carlo SE
#> 0.021; 500 replications). Across replications, the probability of concluding
#> a meaningful effect was 0.678, of detection without meaningfulness 0.176, of
#> an inconclusive result 0.146, and of equivalence 0.000. Estimation used
#> maximum likelihood with delta-method Wald confidence intervals; replications
#> without a maximum-likelihood estimate, from a failed fit or from separation
#> (a focal cell with no events or no non-events), counted against the claim
#> (0.0% of replications).
# }
```
