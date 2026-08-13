# Power curve over a grid of sample sizes

Runs the power simulation at each element of `n` and collects the
results. Works for APE and AIE designs alike.

## Usage

``` r
ape_curve(
  dgp,
  n,
  claim = c("minimum", "detect", "equivalence"),
  sesoi = NULL,
  conf = 0.95,
  nsim = 1000,
  seed = NULL
)

# S3 method for class 'powerape_curve'
plot(x, target_power = NULL, ...)
```

## Arguments

- dgp:

  A `powerape_dgp` after
  [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md)
  or
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md).

- n:

  Integer vector (length \>= 2) of total sample sizes.

- claim:

  `"minimum"` (default), `"detect"`, or `"equivalence"`.

- sesoi:

  Smallest effect size of interest, in APE units. Required for
  `"minimum"` and `"equivalence"`; optional for `"detect"` (if supplied,
  the outcome table is still broken out against it).

- conf:

  CI level (default 0.95).

- nsim:

  Replications per grid point.

- seed:

  Optional; grid point i uses `seed + i - 1`.

- x:

  A `powerape_curve` object.

- target_power:

  Optional horizontal reference line; when the curve crosses it, the
  crossing n is marked and annotated.

- ...:

  Passed to the base plot call.

## Value

A `powerape_curve` object with a `results` data frame (`n`, `power`,
`mcse`, `failed`).

## Examples

``` r
# \donttest{
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
d <- set_ape(d, target = 0.10)
cv <- ape_curve(d, n = c(300, 600, 900), claim = "detect", nsim = 300, seed = 1)
plot(cv, target_power = 0.8)

# }
```
