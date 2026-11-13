# Recompute the true three-way interaction effect implied by a pinned DGP

Recompute the true three-way interaction effect implied by a pinned DGP

## Usage

``` r
true_aie3(dgp)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_aie3()`](https://jespernwulff.github.io/powerape/reference/set_aie3.md).

## Value

The true three-way average interaction effect (numeric scalar),
evaluated on the deterministic integration set; equals the
[`set_aie3()`](https://jespernwulff.github.io/powerape/reference/set_aie3.md)
target up to root-finding tolerance.
