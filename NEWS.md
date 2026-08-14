# powerape 1.2.1

* Fixed: for **logit** panel DGPs with a pure random-effects component
  (`rho > 0` and `cre_share < 1`), calibration and true values integrated
  the unit effect with the probit closed form (index divided by
  `sqrt(1 + var_a)`), which over-attenuates a logistic kernel. A requested
  baseline .30 / APE .10 logit world at `rho = .5, cre_share = 0` actually
  had baseline .269 and true APE .116, so simulated power for such designs
  was optimistic. The unit effect is now integrated by 20-node
  Gauss-Hermite quadrature for logit; probit keeps the exact closed form
  and is numerically unchanged. Realized-vs-requested calibration is
  regression-tested for both links.

# powerape 1.2.0

* Panel designs: `ape_dgp_panel()` specifies a correlated random effects
  (CRE) probit/logit world -- Mundlak heterogeneity on observed unit means,
  latent unit-effect share `rho` (xtprobit convention), `cre_share` for the
  heterogeneity-regressor correlation, and per-variable within-unit
  persistence via `pa_var(..., icc =)` (icc = 1 gives unit-level,
  time-constant variables). Estimation in the engine is pooled probit/logit
  with Mundlak means and unit-clustered standard errors; the target
  remains the ASF-based APE, which the pooled CRE estimator recovers even
  though coefficients are attenuated. Sample-size arguments count **units
  (clusters)**; each contributes `n_periods` observations.
* Clustered delta-method inference: score-based cluster-robust sandwich
  matching R's `sandwich::vcovCL` convention exactly (verified against
  marginaleffects with `vcov = ~id`); Stata's `vce(cluster)` differs only
  by its observed-Hessian bread (~0.3% in the validation example, with the
  attribution reproduced to 7 decimals). A warning is issued below 30
  clusters.
* Validation: exact-enumeration reduction (no-Mundlak case), the
  Donner-Klar cluster design-effect anchor, clustered-CI coverage, and the
  documented "Mundlak insurance premium" (at rho = 0 the CRE estimator is
  correctly less powerful than a cross-sectional analysis of the same
  information). Battery gains V8.
* Panel scope in this release: parametric route, balanced panels, no
  moderator (AIE) combination yet.

# powerape 1.1.0

* `ape_n()` gains a high-precision confirmation stage: by default
  (`confirm = TRUE`) the candidate n is re-measured at
  `nsim_confirm = 4 * nsim` replications and accepted only if the confirmed
  power is within `max(0.005, 1.5 * MCSE)` below the goal, pushing n upward
  otherwise (overshoot is never trimmed -- the conservative direction).
  The returned power and MCSE come from the confirmation run; the search
  history gains a `stage` column. `confirm = FALSE` restores the fast
  single-stage search.
* Validation battery added (`tests/testthat/test-exact-power.R` and the
  project-level battery): simulated power matches the *exact* finite-sample
  power of the engine's decision rule -- computed by triple binomial
  enumeration in the saturated case -- to within Monte Carlo error, for
  detection, minimum-effect, and equivalence claims. Cross-checks against
  `power.prop.test()`, Stata `power twoproportions`, `pwr` (arcsine), and
  Hsieh's logistic-regression formulas (`powerMediation`) agree.

# powerape 1.0.0

* Pilot-model route: `ape_dgp_from_fit()` lifts covariate rows and nuisance
  coefficients from a fitted binomial `glm`; the focal coefficient is
  re-solved in APE units by `set_ape()`. Optional `baseline` override.
* Four vignettes: the minimum-effect workflow, the classic SESOI detection
  analysis, equivalence testing, and average interaction effects.
* Validation: APE point estimates and delta-method SEs match Stata
  `margins, dydx()` to ~4e-7 on fixed data.
* pkgdown site and hex logo.

# powerape 0.3.0

* Average interaction effects: `set_aie()` pins a moderated DGP from three
  explicit anchors (conditional-at-reference main-effect APEs plus the
  target AIE); the internal AIE estimator covers all four variable-type
  pairs and reproduces `ginteff` exactly (estimates and SEs) on fixed data.
* `ape_robust()`: scenario sweeps over contextual assumptions with
  pin-the-APE / pin-the-coefficients modes, worst-case power, and the
  insurance-premium sample size n_max (Hancock & Feng, 2025).
* `power_statement()`: renders any result as a citable methods paragraph.

# powerape 0.2.0

* Continuous (normal) focal variables: average-derivative APE with
  plateau-aware inversion and feasibility reporting.
* Empirical-covariate route: `ape_dgp_empirical()` resamples pilot rows
  (focal jointly from a column, or independently as a `pa_var`).
* Validation: APE estimates and SEs match `marginaleffects` on fixed data;
  simulated power matches Stata `power twoproportions` analytics.

# powerape 0.1.0

* Initial engine: probit/logit, binary focal variable, Gaussian-copula
  covariates, APE-unit inversion with feasibility checks, CI-based claims
  (`detect`, `minimum`, `equivalence`; Riesthuis, 2024), `ape_power()`,
  `ape_curve()`, `ape_n()`, Monte Carlo SEs throughout.
