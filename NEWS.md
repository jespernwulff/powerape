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
