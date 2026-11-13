# Pin the true effect on the treated of a difference-in-differences design

Solves for the treatment coefficient such that the design in `dgp` has
an average treatment effect on the treated of exactly `target`
probability points in the post period: the treated group's post-period
rate with treatment minus the rate it would have had without, the latter
fixed by `baseline`, `gap`, `trend`, and the `trends` scale of
[`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md).

## Usage

``` r
set_att(dgp, target)
```

## Arguments

- dgp:

  A `powerape_dgp` from
  [`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md).

- target:

  Assumed true ATT (the planning value), e.g. `0.05` for five
  probability points; for a design with a subgroup, the two subgroups'
  ATTs, e.g. `c(0.10, 0.04)`.

## Value

The `dgp` with the treatment coefficient, `target_est`, and the treated
group's post-period rate filled in.

## Details

For a design with a `subgroup`, `target` holds the two subgroups'
effects, reference subgroup first, and the pinned estimand is their
difference (second minus first).

## Examples

``` r
d <- ape_dgp_did(baseline = 0.30, gap = 0.10, trend = 0.05)
d <- set_att(d, target = 0.08)
true_att(d)
#> [1] 0.08
```
