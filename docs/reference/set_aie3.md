# Pin the three-way interaction effect of a design with two moderators

For a design built with `ape_dgp(moderator = , moderator2 = )`, solves
all seven effect coefficients from targets in probability points. The
pinned estimand is the **three-way average interaction effect** (as in
the `ginteff` package): the two-way average interaction effect of the
focal variable and the moderator when the second moderator is set to 1
for everyone, minus the same when it is set to 0, each averaged over the
covariates. With three binary variables this is the triple difference of
the eight counterfactual cell rates. The effect is symmetric in its
three variables.

## Usage

``` r
set_aie3(dgp, target, main, two_way)
```

## Arguments

- dgp:

  A `powerape_dgp` from
  [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
  with a `moderator` and a `moderator2`.

- target:

  Assumed true three-way average interaction effect (the planning
  value), in probability points.

- main:

  Numeric vector of three: the main-effect anchors (focal, moderator,
  second moderator).

- two_way:

  Numeric vector of three: the two-way anchors (focal x moderator, focal
  x second moderator, moderator x second moderator).

## Value

The `dgp` with all effect coefficients and `target_est` filled in;
`estimand` becomes `"aie3"`.

## Details

All anchors are explicit; there are no silent defaults:

- `main`: the average partial effect of each variable with the other two
  held at their reference values (0 for binary, the mean for
  continuous), in the order focal, moderator, second moderator.

- `two_way`: the two-way average interaction effect of each pair with
  the third variable held **at its reference value** (not averaged over
  it), in the order focal x moderator, focal x second moderator,
  moderator x second moderator.

- `target`: the three-way effect. The focal-by-moderator effect is
  `two_way[1]` when the second moderator is 0 and `two_way[1] + target`
  when it is 1.

A two-way anchor of 0 says that the pair does not interact in
probability points when the third variable is at its reference value. It
does not say that the product's coefficient is 0: in a probit or logit,
main effects alone produce interaction effects on the probability scale,
of every order. Errors report the attainable range when a target is
infeasible.

Three-way designs run through
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
and
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
(which searches the three-way effect and holds the six other anchors),
under the index analysis or, with
[`set_analysis()`](https://jespernwulff.github.io/powerape/reference/set_analysis.md),
a linear probability analysis whose estimate is the coefficient of the
triple product. A three-way effect is expensive: with balanced binary
variables and the same effect in probability points, each interaction
order needs about four times the observations of the order below.

## Computing time

A design with three binary variables and no covariates is pinned at
once. With covariates the anchors are integrals over the covariate
distribution, and at the default integration size of
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
(`n_int = 4e5`) pinning takes about 15 seconds with three binary
variables, about 40 seconds when the focal variable or the moderator is
continuous, and about a minute when both are (effects of continuous
variables level off as their coefficient grows, so their coefficients
are found by a grid scan).
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
solves only the triple product's coefficient again at each candidate
effect (about 10 seconds with a continuous variable). A smaller `n_int`
speeds up exploration.

## Examples

``` r
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
             moderator = pa_var("female", "binary", p = 0.5),
             moderator2 = pa_var("public", "binary", p = 0.5),
             baseline = 0.30)
d <- set_aie3(d, target = 0.10, main = c(0.08, 0.04, 0.02),
              two_way = c(0.03, 0, 0))
d
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50)
#>   moderator: female (binary, p = 0.50)
#>   second moderator: public (binary, p = 0.50) (all pairwise products and the triple product included)
#>   covariates: none
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.000
#>   true three-way AIE (treat x female x public) pinned at +0.1000 (beta = 0.2526)
#>   main APEs at reference: treat +0.080, female +0.040, public +0.020
#>   two-way AIEs with the third variable at reference: treat x female +0.030, treat x public +0.000, female x public +0.000
#>   cell rates P(Y=1 | treat, female, public): 000 0.300 | 100 0.380 | 010 0.340 | 110 0.450 | 001 0.320 | 101 0.400 | 011 0.360 | 111 0.570
true_aie3(d)
#> [1] 0.1
```
