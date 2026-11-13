# powerape: Power Analysis for Average Partial Effects

Simulation-based power analysis in which the effect size of interest is
the average partial effect (APE) or the average interaction effect (AIE)
from a nonlinear model such as probit or logit, including
smallest-effect-size-of-interest (SESOI) powering via
confidence-interval criteria for detection, minimum-effect, and
equivalence claims (Riesthuis, 2024) and robustness sweeps over
contextual assumptions (Hancock and Feng, 2026).

## Details

Designs are built by
[`ape_dgp()`](https://jespernwulff.github.io/powerape/reference/ape_dgp.md)
(cross-sections, with an optional moderator or two),
[`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md)
and
[`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md)
(pilot data),
[`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)
(panels),
[`ape_dgp_iv()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_iv.md)
(endogenous focal variables), and
[`ape_dgp_did()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did.md)
and
[`ape_dgp_did_pilot()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_did_pilot.md)
(difference-in-differences, with the effect on the treated in
probability points as the target). The effect is pinned with
[`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md),
[`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md),
[`set_aie3()`](https://jespernwulff.github.io/powerape/reference/set_aie3.md)
or
[`set_att()`](https://jespernwulff.github.io/powerape/reference/set_att.md),
and
[`set_analysis()`](https://jespernwulff.github.io/powerape/reference/set_analysis.md)
switches a design to a linear probability analysis for comparison.
[`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
[`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
[`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
and
[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
compute power, required sample sizes, detectable effects and
sensitivity;
[`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
justifies the error rate, and
[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
writes the methods paragraph.

## See also

Useful links:

- <https://jespernwulff.github.io/powerape/>

- <https://github.com/jespernwulff/powerape>

- Report bugs at <https://github.com/jespernwulff/powerape/issues>

## Author

**Maintainer**: Jesper Wulff <jwulff@econ.au.dk>
