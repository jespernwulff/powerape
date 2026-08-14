# Declare a variable of the data-generating process

Describes the marginal distribution of one regressor (the focal variable
or a covariate) for
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md).
Dependence between variables is set via the `correlation` argument of
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
(a Gaussian copula).

## Usage

``` r
pa_var(
  name,
  type = c("binary", "normal"),
  p = NULL,
  mean = 0,
  sd = 1,
  icc = NULL
)
```

## Arguments

- name:

  Variable name (a string; used in printing and error messages).

- type:

  `"binary"` (Bernoulli) or `"normal"`.

- p:

  Success probability, strictly between 0 and 1 (binary only).

- mean, sd:

  Mean and standard deviation (normal only; `sd > 0`).

- icc:

  Panel designs only
  ([`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)):
  within-unit correlation of the variable over time, via a
  shared-component device on the latent normal scale
  (`u_it = sqrt(icc) v_i + sqrt(1-icc) e_it`) before the marginal
  transform. `icc = 0` redraws each period; `icc = 1` makes the variable
  time-constant (e.g., a unit-level treatment). For binary variables the
  realized outcome-scale correlation is attenuated relative to `icc`, as
  with the Gaussian-copula `correlation`. Ignored by the cross-sectional
  constructors.

## Value

An object of class `pa_var`.

## Examples

``` r
pa_var("treat", "binary", p = 0.5)
#> $name
#> [1] "treat"
#> 
#> $type
#> [1] "binary"
#> 
#> $p
#> [1] 0.5
#> 
#> $mean
#> [1] 0
#> 
#> $sd
#> [1] 1
#> 
#> $icc
#> NULL
#> 
#> attr(,"class")
#> [1] "pa_var"
pa_var("age", "normal", mean = 45, sd = 12)
#> $name
#> [1] "age"
#> 
#> $type
#> [1] "normal"
#> 
#> $p
#> NULL
#> 
#> $mean
#> [1] 45
#> 
#> $sd
#> [1] 12
#> 
#> $icc
#> NULL
#> 
#> attr(,"class")
#> [1] "pa_var"
```
