# Recompute the true effect on the treated implied by a pinned design

Recompute the true effect on the treated implied by a pinned design

## Usage

``` r
true_att(dgp)
```

## Arguments

- dgp:

  A `powerape_dgp` from
  [`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md)
  after
  [`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md).

## Value

The true ATT (numeric scalar); equals the
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md)
target up to root-finding tolerance. For a design with a subgroup: the
two subgroups' ATTs and their difference.
