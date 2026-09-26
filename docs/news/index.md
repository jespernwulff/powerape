# Changelog

## powerape 1.11.0

Fixes from the September 2026 blindspot audit. Numbers from designs
without rare outcomes, small focal groups, or the affected options are
unchanged: every seeded result of the article’s worked designs
regenerates bit for bit, apart from
[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
(fixed below).

- **Separation is now detected, reported, and by default scored as a
  failed fit** (new argument `separation = c("fail", "keep")` in
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
  [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
  [`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md),
  and
  [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md);
  DESIGN.md section 16). When a simulated study has a cell that
  identifies the effect – a level of a binary focal variable, or for AIE
  designs a focal-by-moderator cell – with no events or no non-events,
  the maximum-likelihood effect does not exist. Stata’s probit/logit
  drop the perfect predictor and `margins` reports the effect as not
  estimable; R’s [`glm()`](https://rdrr.io/r/stats/glm.html) instead
  stops silently at a large finite coefficient whose delta-method
  standard error collapses to the other cell’s binomial SE. Through
  1.10.0 the engine kept such replications, contrary to its
  documentation, and they usually counted as detections: power was
  overstated, and required sample sizes understated, in rare-outcome
  designs with small focal groups (by up to about half in the audit’s
  unbalanced examples), with `failed 0.000` printed throughout. They now
  count as failed (`"fail"`, the default, matching the documentation and
  Stata’s default analysis); `separation = "keep"` restores the 1.10.0
  behavior (R’s [`glm()`](https://rdrr.io/r/stats/glm.html) +
  marginaleffects analysis, Stata’s `asis`). Either way results report
  the share of simulated studies with separation (`separated`) and the
  average smallest identifying-cell count of events or non-events
  (`min_cell`); print methods and
  [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  state both, and the functions warn when the share exceeds 1% or the
  smallest cell averages fewer than 10 events. A backup check flags
  quasi-separated fits the cell check cannot see (a focal-side
  coefficient with an absurd standardized standard error). The
  exact-enumeration certificate (`tests/testthat/helper-exact.R`,
  battery V1) now enumerates either convention and is exercised where
  the conventions differ.
- **Panel and IV routes: a continuous focal is now centered in the
  simulated outcome index**, as it already was in calibration,
  inversion, and truth evaluation. With a focal mean other than 0 the
  simulation used to run in a world with a different baseline and APE
  than the ones pinned and reported (e.g. an APE about 8% larger and a
  baseline of .365 instead of .300 at focal mean 2 on the panel route;
  for IV-AIE designs the moderator’s main effect was also displaced).
  Binary focals and mean-zero continuous focals are unaffected bit for
  bit; shifting the focal’s location now leaves every draw and estimate
  unchanged (unit-tested to 1e-6).
- **[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
  gains `direction`** (`"positive"` or `"negative"`, defaulting to the
  sign of the pinned effect): it used to search increases only and
  silently re-pinned a DGP pinned at a negative effect to a positive
  MDE. With a binary outcome the direction matters away from a 0.5
  baseline (a decrease is harder to detect than an increase of the same
  size above 0.5, easier below it); the returned `mde` is signed, and
  print and
  [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  name the direction. `ape_robust(mode = "mde")` passes `direction`
  through and takes the MDE largest in magnitude as its worst case.
- **[`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
  re-measures the standard error at its first proposal** before
  searching. The pilot SE was taken at a reference effect (half the
  distance from the baseline to the nearer bound), and at large `n` –
  where the MDE lies far below that reference – the first candidate
  overshot by 2-4% at a 0.30 baseline (more at rarer baselines) and was
  accepted within the search tolerance without ever being trimmed. MDE
  results change accordingly (the confirmation stage still pushes only
  outward).
- **[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)’s
  `n_max` is now the largest requirement across scenarios.** It used to
  be solved in the scenario with the lowest simulated power at the
  user’s `n`; when every scenario’s power saturated near 1, that tie
  went to the first grid row (often the most favorable scenario), so
  `n_max` could be understated by about half. Scenarios are now ranked
  by a criterion that does not saturate (the normal-approximation
  requirement implied by each scenario’s simulated standard error),
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  runs in the least favorable scenario and in any within 5% of it (up to
  three), and the result reports the scenario that sets `n_max`
  (`nmax_scenario`) and every candidate searched (`nmax_candidates`).
  When a scenario’s implied effect sits on or beyond the claim boundary
  (possible with `pin = "coefficients"`), `n_max` is withheld with a
  warning. Ties in the worst-power scenario are broken by the same
  criterion. The documented `nmax_power = 0.90` default is unchanged;
  the vignette now says so.
- **[`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md)’s
  focal guard is structural.** It checked the focal’s *variable* name
  while addressing its model-matrix *column*, so interactions of a
  logical or factor focal, transformed copies of the focal (`I(dose^2)`
  next to `dose`, `log(dose)`, `I(treat * female)`), multi-column terms,
  and offsets were accepted silently – and the pinned “APE” was then a
  partial effect holding the focal’s other columns fixed, not the APE
  that margins or marginaleffects report for the same model. All of
  these are now refused with an explanation; a logical or two-level
  factor focal entering as a single main-effect dummy is accepted.
- **[`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  describes the world it priced**: the latent (Gaussian-copula)
  correlation actually specified (or that the regressors are independent
  – it used to assert “Gaussian-copula dependence” even without one),
  the covariates’ marginals, each panel variable’s within-unit
  persistence, and the empirical route’s resampling mode, plus the
  separation convention.
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  now stores its standard-error type (`se`), so statements for robust-SE
  sample-size results no longer describe model-based standard errors;
  [`print()`](https://rdrr.io/r/base/print.html) shows it.

## powerape 1.10.0

- **[`alpha_frontier()`](https://jespernwulff.github.io/powerape/reference/alpha_frontier.md)**:
  the closed-form justified error rate at the detectability frontier –
  the level that minimizes the weighted combined Type I and Type II
  error rate (Mudge et al., 2012; Maier & Lakens, 2022) for the effect
  that is just detectable at the target power. Under the normal
  approximation the first-order condition collapses to
  `phi(z) = prior * phi(qnorm(power)) / (m * cost)` (`m = 2` two-sided
  detection, `m = 1` one-sided minimum-effect), so the level involves
  neither `n` nor the standard error: solve it once and price a whole
  MDE grid under one justified alpha
  (`ape_mde(..., alpha = alpha_frontier(cost = 4))`). Anchors: 0.0274 at
  Cohen’s 4:1 weighting and 80% power (two-sided); exactly `1 - power`
  at equal costs (one-sided), where the frontier balances alpha = beta.
  Errors when the weighting admits no interior optimum and warns when
  the level is not below 0.5. Battery check V22 verifies the closed form
  against
  [`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)’s
  simulation optimum on frontier-pinned draws.
- **[`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
  gains the `se` argument** of
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md):
  a design analyzed with heteroskedasticity-robust standard errors now
  prices its robustness sweep (and the `n_max` search and MDE mode)
  under the same standard errors instead of silently reverting to
  model-based ones. Panel and IV routes keep their fixed inference and
  warn, as in
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md).
- **Validation battery V21** covers the empirical-covariates and
  pilot-model routes, which previously had unit tests but no battery
  item: 95% CI coverage at the nominal rate, estimator consistency, and
  concordance with the matched parametric world (a large pilot drawn
  *from* that world must reproduce its power within Monte Carlo error)
  for both
  [`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md)
  and
  [`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md).
- Input-validation hardening: a `sesoi` supplied alongside
  `claim = "detect"` (it breaks out the outcome table) is now checked
  for being a single positive number instead of silently corrupting that
  table;
  [`pa_var()`](https://jespernwulff.github.io/powerape/reference/pa_var.md)
  warns when arguments irrelevant to the declared type are supplied (`p`
  for a normal variable, `mean`/`sd` for a binary one).
- [`plot.powerape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md)
  no longer forces isotonic smoothing when `alpha` is a sample-size rule
  – under an n-dependent level, power is legitimately non-monotone in
  `n`, and the smoother would misdraw it.
- [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  no longer prints an empty covariate list for an empirical-route DGP
  whose pilot data supply only the focal column.

## powerape 1.9.0

- **Justified error rates**
  ([`vignette("justified-alpha")`](https://jespernwulff.github.io/powerape/articles/justified-alpha.md);
  DESIGN.md section 15). Two routes to an alpha that is chosen rather
  than inherited, both opt-in; the default stays `alpha = 0.05` per
  claim.
  - **Sample-size rules.** `alpha` may be a function of `n` in every
    power function, for example
    `alpha = function(n) alphaN::alphaN(n, BF = 3)`, the Bayes-factor
    calibration of Wulff and Taylor (2024) that lowers alpha as `n`
    grows. The rule is evaluated at each candidate `n`, so
    [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
    designs jointly over the pair (alpha(n), n), and
    [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
    [`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md),
    and
    [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md)
    use the level that belongs to their `n`. Results store the rule and
    the level it produced (`alpha`, `alpha_rule`; curves and search
    histories gain an `alpha` column), print methods say “by sample-size
    rule”, and
    [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
    states it.
  - **Error-cost optimization.**
    [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md)
    now keeps the estimate and standard error of every replication
    (`keep_draws = TRUE`), so power at any alpha is a re-threshold of
    one simulation, exact on those draws.
    [`power_at()`](https://jespernwulff.github.io/powerape/reference/power_at.md)
    exposes that curve.
    [`ape_alpha()`](https://jespernwulff.github.io/powerape/reference/ape_alpha.md)
    minimizes or balances the weighted combined Type I and Type II error
    rate of Mudge et al. (2012) and Maier and Lakens (2022) for given
    error costs and prior odds; it reports the flat range of levels
    indistinguishable from the optimum at the run’s precision, confines
    the search to `cap = 0.05` (raising alpha needs a justification of
    its own) while stating when the unconstrained optimum lies above the
    cap, and with `size = TRUE` simulates the claim-boundary world to
    report the REALIZED error rate at the chosen level, which analytic
    compromise-power tools cannot. `print`, `plot`, and
    [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
    methods; a preregistration paragraph.
- Verification: new unit tests (`test-alpha.R`: constant rules reproduce
  scalar runs bit for bit,
  [`power_at()`](https://jespernwulff.github.io/powerape/reference/power_at.md)
  exactness, a closed-form twin of the optimizer on synthetic draws,
  Cohen’s 4:1 balance returning .05/.20 at 80% power, JustifyAlpha
  concordance on the same draws, the realized boundary size) and battery
  check V20.
- `alphaN` and `JustifyAlpha` join Suggests (examples and concordance
  tests only).

## powerape 1.8.0

- **Breaking: error-rate convention.** Every claim is now tested at
  `alpha` (new argument, default 0.05) in its conventional form.
  Detection stays the two-sided 5% test (95% interval, directional
  counting); the minimum-effect claim is a one-sided test at 5% and
  equivalence a TOST at 5%, both read off the **90%** interval – the
  convention of Lakens (2017), Lakens et al. (2018), TOSTER, and
  Riesthuis (2024) for equivalence. Before 1.8.0 all three claims read
  the 95% interval, which made the SESOI claims one-sided tests at 2.5%
  (a TOST at 2.5%): internally consistent but stricter than the field’s
  convention without a stated reason. Consequences: minimum-effect and
  equivalence power rise and required n falls (Design A’s minimum-claim
  requirement moves from ~2,600 to ~2,100); detection results are
  unchanged. `conf` remains as an explicit override of the claim’s own
  interval level (`conf = 0.95` on a SESOI claim reproduces the old
  numbers). The outcome distribution uses both intervals accordingly.
  [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  now names the claim’s error rate.
- Verification: the TOSTER concordance test and battery V7 now audit the
  DEFAULT convention (powerape’s default must match TOSTER at TOSTER’s
  default alpha with no remapping), and a new battery check V19 measures
  size at the claim boundaries (.05 for the SESOI claims, .025 for
  directional detection). Why the old battery did not flag the issue:
  V1/V7/V13 verified that the engine implements its stated rule exactly
  – V7 even matched TOSTER at a remapped alpha = .025 and guarded that
  .05 “did not match” – so a convention choice could not be caught by a
  convention-matched check. The new checks pin the default to the
  field’s default.

## powerape 1.7.0

- **Unbalanced panels** (Wooldridge, 2019, J. Econometrics).
  [`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)
  accepts a distribution over panel lengths: a vector `n_periods` with
  probabilities `p_periods` (equal weights by default), or a scalar
  `n_periods` with a per-wave `retention` rate for monotone attrition.
  Mundlak means are computed over each unit’s observed periods and the
  estimating model gains period-count cohort indicators (the Wooldridge
  2019 workhorse); selection is completely at random by construction,
  `rho` is the length-mixture average of the latent unit-effect share,
  and sample sizes still count units. Degenerate mixtures
  (`retention = 1`, single support) short-circuit to the balanced code
  path, bit-identically. `retention` is sweepable in
  [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md).
  Validated as battery V18: balanced-reduction identity; Stata
  `probit i.d ... i.T, vce(cluster)` + `margins` agreement to ~5e-7 on
  coefficients and all seven decimals on the APE, with the 0.12%
  clustered-SE gap reproduced to 4e-8 from the observed-information
  bread; an unbalanced `xtprobit` ASF anchor (\|diff\| .0012);
  consistency and 95% coverage; and content anchors (power monotone in
  retention; no unequal-cluster penalty for a mean-preserving length
  mixture at the tested settings).

## powerape 1.6.0

- **Breaking: panel designs are probit-only.**
  `ape_dgp_panel(model = "logit")` is now refused with a teaching error
  instead of building a world. Rationale: the exact ASF recovery that
  justifies the CRE panel route is probit-specific (a normal unit effect
  rescales the index by `sqrt(1 + var_a)`); a pooled CRE *logit*
  estimates only a quasi-ML approximation of the ASF (measured gap in
  the least favorable validated configuration: -0.0002, Monte Carlo SE
  0.001 – small, but an approximation, and `powerape` does not price
  designs with an approximate estimand). Calibrated to the same baseline
  and probability-point target, the probit panel world is practically
  indistinguishable on the probability scale. The 1.2.1 Gauss-Hermite
  calibration machinery remains internally (and tested) for the mixture
  integrals; no public route reaches its logit branch.

## powerape 1.5.1

- Verification expansion (battery V12-V17, adversarially refereed):
  exact enumeration of the *clustered* panel decision rule (within half
  a Monte Carlo SE of simulation at 24 and 40 clusters); a size audit
  across all routes (nominal at the standard/panel routes and at the
  minimum-claim boundary; the small-G cluster over-rejection and the IV
  route’s conservative null quantified and decomposed); a recursive
  bivariate-probit anchor separating the binary-endogenous control
  function’s plim gap (-0.011) from the truth; an asymptotic-SE identity
  for the CF route; xtprobit agreement (pure-RE and Mundlak-augmented)
  to ~0.001; and bootstrap-vs-delta SE ratios within 2% on all three
  routes. New stochastic consistency/coverage test for the IV-AIE.
- Documentation:
  [`ape_dgp_iv()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_iv.md)
  now records two properties of the control-function estimator found in
  validation – the conservative null of the binary-endogenous
  generalized residual, and the observed-pairs vs product-measure ASF
  wedge for continuous endogenous focal variables (the margins
  convention, replicated exactly; a product-measure option is on the
  roadmap).

## powerape 1.5.0

- Minimum detectable effect:
  [`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
  inverts the design question – the sample size is fixed (an archive of
  so many firms, a panel of so many units and waves) and the function
  finds the smallest APE or AIE the design reliably concludes at a
  target power. Detect and minimum-effect claims search over the assumed
  effect, re-running the APE/AIE inversion at every candidate; the
  equivalence claim holds the pinned truth and returns the tightest
  establishable margin instead. Every route works unchanged (parametric,
  empirical, pilot-model, panel with n in units, IV), and the answer is
  verified with the same conservative confirmation stage as
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  (pushed upward, never down, when the high-precision run falls short).
  Feasibility ceilings produce a teaching error reporting the largest
  attainable effect and the power available there.
- Robustness for the MDE: `ape_robust(mode = "mde")` sweeps the
  contextual assumptions and reports the minimum detectable effect per
  scenario; the worst case is the *largest* MDE – the smallest effect
  the design finds even under the least favorable assumptions.
- [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  renders
  [`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md)
  results (“the minimum detectable APE at n = … was …”), units-aware for
  panel designs.
- Validation: battery gains V11 – the MDE at the four-tool anchor (n
  = 712) returns the anchor’s target (~0.10), the exact-enumeration
  power of the decision rule at the returned MDE sits at the goal, and
  the analytic inverse from
  [`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html)
  agrees.
- The minimum-effect and panel vignettes gain fixed-n sections
  demonstrating
  [`ape_mde()`](https://jespernwulff.github.io/powerape/reference/ape_mde.md).

## powerape 1.4.0

- IV designs via control-function probit:
  [`ape_dgp_iv()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_iv.md)
  specifies a world with an **endogenous focal variable** – linear first
  stage for a continuous focal, probit first stage for a binary one,
  instruments with a friendly strength knob (`iv_strength` = the share
  of the focal’s (latent) variance the instruments explain), and
  `endogeneity` = the structural/first-stage error correlation.
  Estimation in the engine replicates **Stata 18.5’s `cfprobit`**
  exactly: two-step control function with stacked method-of-moments
  standard errors (no bootstrap; robust and cluster flavors verified
  against Stata to ~1e-7, including the convention details – observed
  score derivatives, generated regressors as fixed instruments in the
  first-stage Jacobian, no small-sample factors). The estimand stays the
  ASF-based APE/AIE, which `margins` after `cfprobit` targets with the
  control functions held fixed.
- IV interaction designs: with a binary exogenous `moderator`, the main
  equation gains the focal-by-moderator term, the control function, and
  the control-function-by-moderator interaction (Stata’s `interact()` +
  `mainonly()`), and the first stage includes the moderator **and the
  instrument-by-moderator interactions** – the general IV requirement
  for interaction models. Validated against Stata’s exact syntax.
- The estimator internals also cover `fprobit` and `poisson` first
  stages, matching `cfprobit`’s fit exactly. For fractional first stages
  we found Stata’s `margins` to be internally inconsistent with its own
  estimator (its `_remake_cfs.ado` rebuilds the control function with
  the binary-outcome branch, ignoring the fractional value); powerape
  uses the estimation-consistent generalized residuals, and the
  discrepancy is reproduced to 7 decimals in the validation battery.
- Robust standard errors for the standard routes:
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
  [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
  and
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  gain `se = c("model", "robust")` – heteroskedasticity-robust (HC0)
  sandwich matching `sandwich::vcovHC(type = "HC0")` exactly. Panel
  designs keep their unit-clustered SEs and IV designs their stacked
  robust sandwich (a warning explains this if `se` is passed there).
- Validation: battery gains V10 (the seven-case Stata equivalence table,
  the fprobit attribution, engine consistency and coverage, and the “IV
  price” anchor – at zero endogeneity the CF estimator’s SE exceeds the
  exogenous estimator’s by about sqrt(1/iv_strength), the 2SLS variance
  logic).

## powerape 1.3.0

- Panel AIE designs (binary x binary):
  [`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)
  gains a `moderator` argument, and
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
  pins the average interaction effect of a CRE panel world. All anchors
  – baseline, the two conditional-at-reference main-effect APEs, and the
  target AIE – are defined on the average structural function (unit
  effects integrated exactly, Mundlak means held fixed), so the panel
  AIE is the same population quantity a cross-sectional AIE design
  targets. The estimating model includes the interaction’s own Mundlak
  mean whenever the product term is time-varying (the mean of a product
  is not the product of means); inference is unit-clustered throughout,
  and
  [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md),
  [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
  and
  [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md)
  all work unchanged.
- Validation: the ASF double difference matches Stata
  (`probit i.d##i.m ... , vce(cluster id)` + `margins` + `lincom`) to
  all printed digits on a fixed panel, with the clustered-SE gap again
  fully attributed to Stata’s observed-Hessian bread (reproduced to 7
  decimals); ginteff reproduces the double difference and model-based SE
  on the identical fit to 1e-6; exact reduction of the inversion to the
  cross-sectional case at rho = 0; estimator consistency and nominal
  clustered-CI coverage of the true AIE. Battery gains V9.
- Continuous focal or moderator pairs in panel AIE designs are not yet
  supported (clean error); cross-sectional AIE keeps all four type
  pairs.
- New vignette `panel-designs`: the CRE workflow end to end – the panel
  knobs, units-not-observations accounting, the panel-vs-cross-section
  comparison, rho robustness, and the panel AIE. The README now covers
  the panel route.

## powerape 1.2.1

- Fixed: for **logit** panel DGPs with a pure random-effects component
  (`rho > 0` and `cre_share < 1`), calibration and true values
  integrated the unit effect with the probit closed form (index divided
  by `sqrt(1 + var_a)`), which over-attenuates a logistic kernel. A
  requested baseline .30 / APE .10 logit world at
  `rho = .5, cre_share = 0` actually had baseline .269 and true APE
  .116, so simulated power for such designs was optimistic. The unit
  effect is now integrated by 20-node Gauss-Hermite quadrature for
  logit; probit keeps the exact closed form and is numerically
  unchanged. Realized-vs-requested calibration is regression-tested for
  both links.

## powerape 1.2.0

- Panel designs:
  [`ape_dgp_panel()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_panel.md)
  specifies a correlated random effects (CRE) probit/logit world –
  Mundlak heterogeneity on observed unit means, latent unit-effect share
  `rho` (xtprobit convention), `cre_share` for the
  heterogeneity-regressor correlation, and per-variable within-unit
  persistence via `pa_var(..., icc =)` (icc = 1 gives unit-level,
  time-constant variables). Estimation in the engine is pooled
  probit/logit with Mundlak means and unit-clustered standard errors;
  the target remains the ASF-based APE, which the pooled CRE estimator
  recovers even though coefficients are attenuated. Sample-size
  arguments count **units (clusters)**; each contributes `n_periods`
  observations.
- Clustered delta-method inference: score-based cluster-robust sandwich
  matching R’s
  [`sandwich::vcovCL`](https://sandwich.R-Forge.R-project.org/reference/vcovCL.html)
  convention exactly (verified against marginaleffects with
  `vcov = ~id`); Stata’s `vce(cluster)` differs only by its
  observed-Hessian bread (~0.3% in the validation example, with the
  attribution reproduced to 7 decimals). A warning is issued below 30
  clusters.
- Validation: exact-enumeration reduction (no-Mundlak case), the
  Donner-Klar cluster design-effect anchor, clustered-CI coverage, and
  the documented “Mundlak insurance premium” (at rho = 0 the CRE
  estimator is correctly less powerful than a cross-sectional analysis
  of the same information). Battery gains V8.
- Panel scope in this release: parametric route, balanced panels, no
  moderator (AIE) combination yet.

## powerape 1.1.0

- [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md)
  gains a high-precision confirmation stage: by default
  (`confirm = TRUE`) the candidate n is re-measured at
  `nsim_confirm = 4 * nsim` replications and accepted only if the
  confirmed power is within `max(0.005, 1.5 * MCSE)` below the goal,
  pushing n upward otherwise (overshoot is never trimmed – the
  conservative direction). The returned power and MCSE come from the
  confirmation run; the search history gains a `stage` column.
  `confirm = FALSE` restores the fast single-stage search.
- Validation battery added (`tests/testthat/test-exact-power.R` and the
  project-level battery): simulated power matches the *exact*
  finite-sample power of the engine’s decision rule – computed by triple
  binomial enumeration in the saturated case – to within Monte Carlo
  error, for detection, minimum-effect, and equivalence claims.
  Cross-checks against
  [`power.prop.test()`](https://rdrr.io/r/stats/power.prop.test.html),
  Stata `power twoproportions`, `pwr` (arcsine), and Hsieh’s
  logistic-regression formulas (`powerMediation`) agree.

## powerape 1.0.0

- Pilot-model route:
  [`ape_dgp_from_fit()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_from_fit.md)
  lifts covariate rows and nuisance coefficients from a fitted binomial
  `glm`; the focal coefficient is re-solved in APE units by
  [`set_ape()`](https://jespernwulff.github.io/powerape/reference/set_ape.md).
  Optional `baseline` override.
- Four vignettes: the minimum-effect workflow, the classic SESOI
  detection analysis, equivalence testing, and average interaction
  effects.
- Validation: APE point estimates and delta-method SEs match Stata
  `margins, dydx()` to ~4e-7 on fixed data.
- pkgdown site and hex logo.

## powerape 0.3.0

- Average interaction effects:
  [`set_aie()`](https://jespernwulff.github.io/powerape/reference/set_aie.md)
  pins a moderated DGP from three explicit anchors
  (conditional-at-reference main-effect APEs plus the target AIE); the
  internal AIE estimator covers all four variable-type pairs and
  reproduces `ginteff` exactly (estimates and SEs) on fixed data.
- [`ape_robust()`](https://jespernwulff.github.io/powerape/reference/ape_robust.md):
  scenario sweeps over contextual assumptions with pin-the-APE /
  pin-the-coefficients modes, worst-case power, and the
  insurance-premium sample size n_max (Hancock & Feng, 2025).
- [`power_statement()`](https://jespernwulff.github.io/powerape/reference/power_statement.md):
  renders any result as a citable methods paragraph.

## powerape 0.2.0

- Continuous (normal) focal variables: average-derivative APE with
  plateau-aware inversion and feasibility reporting.
- Empirical-covariate route:
  [`ape_dgp_empirical()`](https://jespernwulff.github.io/powerape/reference/ape_dgp_empirical.md)
  resamples pilot rows (focal jointly from a column, or independently as
  a `pa_var`).
- Validation: APE estimates and SEs match `marginaleffects` on fixed
  data; simulated power matches Stata `power twoproportions` analytics.

## powerape 0.1.0

- Initial engine: probit/logit, binary focal variable, Gaussian-copula
  covariates, APE-unit inversion with feasibility checks, CI-based
  claims (`detect`, `minimum`, `equivalence`; Riesthuis, 2024),
  [`ape_power()`](https://jespernwulff.github.io/powerape/reference/ape_power.md),
  [`ape_curve()`](https://jespernwulff.github.io/powerape/reference/ape_curve.md),
  [`ape_n()`](https://jespernwulff.github.io/powerape/reference/ape_n.md),
  Monte Carlo SEs throughout.
