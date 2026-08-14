# Pin the true AIE of a moderated DGP (inversion step)

For a DGP with a moderator, solves all three effect coefficients from
targets in effect-size units: the two main effects from their
*conditional-at-reference APEs* (the APE of each variable with the other
held at its reference value – 0 for binary, the mean for continuous),
and the interaction coefficient from the target average interaction
effect (AIE), holding the mains fixed. The AIE is `ginteff`'s estimand:
the double difference for binary-by-binary, the difference of average
derivatives for mixed pairs, and the average cross-partial for
continuous-by-continuous.

## Usage

``` r
set_aie(dgp, target, main_focal, main_moderator)
```

## Arguments

- dgp:

  A `powerape_dgp` from
  [`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
  or
  [`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md),
  with a `moderator`.

- target:

  Assumed true AIE (the planning value), in APE units.

- main_focal, main_moderator:

  Conditional-at-reference APEs anchoring the two main effects.

## Value

The `dgp` with `beta_focal`, `beta_mod`, `beta_int`, and `target_est`
filled in; `estimand` becomes `"aie"`.

## Details

All three anchors are explicit – there are no silent defaults. Errors
report the attainable range when a target is infeasible (for
binary-by-binary, the four counterfactual cell rates must stay in (0,1);
otherwise the AIE plateaus in the interaction coefficient).

For panel DGPs
([`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)
with a `moderator`; binary x binary), all anchors are defined on the
average structural function: unit effects are integrated out and the
Mundlak means held fixed, so the panel AIE is the same population
quantity a cross-sectional AIE design targets.

## Examples

``` r
d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
             moderator = pa_var("female", "binary", p = 0.55),
             baseline = 0.30)
d <- set_aie(d, target = 0.04, main_focal = 0.10, main_moderator = 0.05)
d
#> powerape DGP -- probit
#>   focal: treat (binary, p = 0.50)
#>   moderator: female (binary, p = 0.55) (focal-by-moderator interaction included)
#>   covariates: none
#>   baseline P(Y=1 | reference) = 0.300, nuisance signal (latent pseudo-R2) = 0.000
#>   true AIE pinned at +0.0400 (beta_int = 0.0892)
#>   main APEs at reference: treat +0.100 (beta = 0.2711), female +0.050 (beta = 0.1391)
```
