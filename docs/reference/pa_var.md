# Declare a variable of the data-generating process

Describes the marginal distribution of one regressor (the focal variable
or a covariate) for
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md).
Dependence between variables is set via the `correlation` argument of
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
(a Gaussian copula).

## Usage

``` r
pa_var(name, type = c("binary", "normal"), p = NULL, mean = 0, sd = 1)
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
#> attr(,"class")
#> [1] "pa_var"
```
