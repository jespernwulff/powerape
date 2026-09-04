# Closed-form justified alpha at the detectability frontier

The normal-approximation companion of
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
for pricing minimum detectable effects: the error rate that minimizes
the weighted combined Type I and Type II error rate (Mudge et al., 2012;
Maier & Lakens, 2022) for the effect that is *just detectable* at the
target power. Under the normal approximation, the minimize-mode
first-order condition for an effect-to-standard-error ratio `r` is
`m * cost * phi(z) = prior * phi(z - r)`, with `z` the claim's critical
value and `m = 2` for the two-sided detection claim (`m = 1` for the
one-sided minimum-effect claim). At the detectability frontier the ratio
`r` equals `z + qnorm(power)`, so the condition collapses to
`phi(z) = prior * phi(qnorm(power)) / (m * cost)` – one equation
involving neither the sample size nor the standard error. The frontier
level is therefore the **same at every n**, which makes it the natural
companion of an MDE grid: solve the level once, pass it as `alpha` to
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
(or
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)),
and the whole grid is priced under one justified error rate.

## Usage

``` r
alpha_frontier(
  cost = 1,
  power = 0.8,
  prior = 1,
  claim = c("detect", "minimum")
)
```

## Arguments

- cost:

  Relative cost of a Type I error against a Type II error (default 1).

- power:

  Target power defining the frontier effect (default 0.80).

- prior:

  Prior odds of the just-detectable effect against the claim boundary
  (default 1).

- claim:

  `"detect"` (default; two-sided) or `"minimum"` (one-sided). For
  equivalence run
  [`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
  on stored draws.

## Value

The frontier error rate (a numeric scalar), suitable as the `alpha`
argument of the power functions when it is below 0.5.

## Details

Two anchors worth knowing. At Cohen's 4:1 weighting and 80% power the
two-sided frontier level is 0.0274. At *equal* costs the one-sided
frontier level is exactly `1 - power` (the frontier balances
`alpha = beta` there).

**Relation to
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md).**
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
is exact for a specific design: it optimizes over the stored draws of
one simulation, applies the conventional `cap = 0.05`, and can verify
the realized size. `alpha_frontier()` is its closed-form, scale-free
limit at the just-detectable effect; it applies no cap (levels above .05
need a justification of their own; Maier & Lakens, 2022) and warns when
the optimum does not lie below 0.5, where it could not be used as an
`alpha` anyway.

## References

Maier, M., & Lakens, D. (2022). Justify your alpha: A primer on two
practical approaches. *Advances in Methods and Practices in
Psychological Science, 5*(2). Mudge, J. F., Baker, L. F., Edge, C. B., &
Houlahan, J. E. (2012). Setting an optimal alpha that minimizes errors
in null hypothesis significance tests. *PLOS ONE, 7*(2).

## Examples

``` r
alpha_frontier(cost = 4, power = 0.80)            # 0.0274
#> [1] 0.02737173
alpha_frontier(cost = 1, power = 0.80, claim = "minimum")  # exactly 0.20
#> [1] 0.2
```
