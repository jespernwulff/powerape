# Pin the true APE of a DGP (inversion step)

Solves for the focal coefficient such that the data-generating process
in `dgp` has average partial effect exactly `target` (in probability
points): the average discrete change for a binary focal variable, the
average derivative for a continuous one. Root-finding runs on the
deterministic integration set. Errors if the target is infeasible – for
a binary focal the APE must lie inside `(-p0, 1 - p0)` given baseline
rate p0; for a continuous focal the APE plateaus in the coefficient, and
the attainable range is reported. For DGPs with a moderator, use
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
instead.

## Usage

``` r
set_ape(dgp, target)
```

## Arguments

- dgp:

  A `powerape_dgp` from
  [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
  or
  [`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md),
  without a moderator.

- target:

  Assumed true APE (the planning value), e.g. `0.05` for five
  probability points. Always explicit – there is no default truth.

## Value

The `dgp` with the focal coefficient and `target_est` filled in (and,
for a binary focal, the implied `P(Y = 1 | focal = 1)`).

## Examples

``` r
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
d
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50)
#>   covariates: none
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.000
#>   true APE pinned at +0.1000 (beta_focal = 0.2711, implied P(Y=1 | focal=1) = 0.400)
```
