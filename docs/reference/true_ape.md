# Recompute the true APE implied by a pinned DGP

Recompute the true APE implied by a pinned DGP

## Usage

``` r
true_ape(dgp)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md).

## Value

The true APE (numeric scalar), evaluated on the deterministic
integration set; equals the
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
target up to root-finding tolerance.
