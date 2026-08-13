# Recompute the true AIE implied by a pinned DGP

Recompute the true AIE implied by a pinned DGP

## Usage

``` r
true_aie(dgp)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md).

## Value

The true AIE (numeric scalar), evaluated on the deterministic
integration set; equals the
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
target up to root-finding tolerance.
