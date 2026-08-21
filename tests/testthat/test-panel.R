# Panel (CRE probit) designs: construction, inversion, engine, anchors.

panel_base <- function(rho = 0.3, cre = 0.5, T_p = 4, target = NULL) {
  d <- ape_dgp_panel(
    focal = pa_var("treat", "binary", p = 0.5, icc = 1),
    covariates = list(pa_var("size", "normal", icc = 0.6)),
    n_periods = T_p, rho = rho, cre_share = cre,
    baseline = 0.30, signal = 0.10, n_int = 4e4
  )
  if (!is.null(target)) d <- set_ape(d, target)
  d
}

test_that("panel construction calibrates baseline, rho, and icc as stated", {
  d <- panel_base()
  ## baseline reproduced exactly on the integration draw (a integrated in closed form)
  id <- powerape:::panel_integration(d)
  expect_equal(mean(pnorm((d$beta0 + id$q) / d$scale_a)), 0.30, tolerance = 1e-6)
  ## Var(c) decomposition: mundlak part + var_a = rho/(1-rho)
  var_c <- d$rho / (1 - d$rho)
  expect_equal(d$var_a, (1 - d$cre_share) * var_c, tolerance = 1e-12)
  ## realized within-unit correlation of the normal covariate ~ its icc
  vars <- c(list(d$focal), d$covariates)
  rd <- powerape:::panel_regressor_draw(vars, d$iccs, d$R, 20000, 2, seed = 5)
  z <- matrix(rd$x[[2]], nrow = 2)   # 2 periods per unit, unit-major
  expect_equal(cor(z[1, ], z[2, ]), 0.6, tolerance = 0.02)
  ## time-constant focal really is time-constant
  tr <- matrix(rd$x[[1]], nrow = 2)
  expect_true(all(tr[1, ] == tr[2, ]))
})

test_that("panel guards fire", {
  expect_error(ape_dgp_panel(focal = pa_var("t", "binary", p = .5),
                             n_periods = 1, baseline = .3),
               "at least two periods")
  expect_error(ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                             n_periods = 4, rho = .3, cre_share = .5,
                             baseline = .3),
               "time-varying")
  d <- panel_base()
  expect_error(powerape:::draw_x_panel(d, 100), "set_ape")
})

test_that("panel inversion hits the target and the estimator recovers it", {
  d <- panel_base(target = 0.08)
  expect_equal(true_ape(d), 0.08, tolerance = 1e-7)
  ## consistency: mean of CRE-probit APE estimates equals the target
  sims <- powerape:::sim_ci(d, 400, 300, 0.95, seed = 21)
  expect_gt(mean(sims$ok), 0.99)
  est <- (sims$l + sims$u)[sims$ok] / 2
  expect_equal(mean(est), 0.08, tolerance = 0.01)
  ## clustered-SE calibration: mean SE tracks the sd of estimates
  ratio <- mean(sims$se[sims$ok]) / sd(est)
  expect_lt(abs(ratio - 1), 0.25)
})

test_that("clustered CIs cover the true APE at nominal rate", {
  d <- panel_base(target = 0.08)
  s <- powerape:::sim_ci(d, 250, 1500, 0.95, seed = 31)
  cov <- mean((s$l[s$ok] <= 0.08) & (s$u[s$ok] >= 0.08))
  expect_lt(abs(cov - 0.95), 0.02)
})

test_that("panel with no heterogeneity and no Mundlak terms reduces to the exact saturated case", {
  ## icc = 1 focal, no covariates, rho = 0: no Mundlak columns and iid rows,
  ## so the exact enumeration applies at n = N*T (up to the coarser
  ## unit-level assignment of treatment).
  T_p <- 3
  N <- 250
  dp <- ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                      n_periods = T_p, rho = 0, cre_share = 0,
                      baseline = 0.30, n_int = 4e4)
  dp <- set_ape(dp, 0.10)
  pw <- ape_power(dp, n = N, claim = "detect", nsim = 2000, seed = 41)
  ex <- exact_power_sat(N * T_p, 0.5, 0.30, 0.40, 0.95, "detect")
  expect_lt(abs(pw$power - ex), 4 * sqrt(ex * (1 - ex) / 2000) + 0.01)
})

test_that("the Mundlak specification's efficiency cost at rho = 0 is priced, not hidden", {
  ## With time-varying regressors the CRE estimator holds unit means fixed,
  ## spending ~1/T of the identifying variation. At rho = 0 that is pure
  ## insurance premium: the panel design must show LESS power than the same
  ## information analyzed cross-sectionally -- by design, not by accident.
  dp <- ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5),
                      covariates = list(pa_var("z", "normal")),
                      n_periods = 3, rho = 0, cre_share = 0,
                      baseline = 0.30, signal = 0.10, n_int = 4e4)
  dp <- set_ape(dp, 0.10)
  dc <- set_ape(ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                        covariates = list(pa_var("z", "normal")),
                        baseline = 0.30, signal = 0.10), 0.10)
  pp <- ape_power(dp, n = 250, claim = "detect", nsim = 600, seed = 41)
  pc <- ape_power(dc, n = 750, claim = "detect", nsim = 600, seed = 42)
  expect_lt(pp$power, pc$power)
  expect_lt(pc$power - pp$power, 0.30)
})

test_that("unit-level treatment matches the cluster design-effect anchor", {
  ## Pure RE, unit-level binary treatment, no covariates: classic
  ## cluster-randomized design. Simulated power should track
  ## power.prop.test at effective n = N*T / DEFF, DEFF = 1 + (T-1) ICC_y.
  T_p <- 4
  d <- ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                     n_periods = T_p, rho = 0.3, cre_share = 0,
                     baseline = 0.30, signal = 0, n_int = 4e4)
  d <- set_ape(d, 0.10)

  ## outcome-scale ICC by numerical integration of the RE probit
  icc_y <- function(m, va) {
    s <- sqrt(1 + va)
    p <- pnorm(m / s)
    p11 <- integrate(function(a) pnorm(m + a)^2 * dnorm(a, 0, sqrt(va)),
                     -Inf, Inf, rel.tol = 1e-10)$value
    (p11 - p^2) / (p * (1 - p))
  }
  i0 <- icc_y(d$beta0, d$var_a)
  i1 <- icc_y(d$beta0 + d$beta_focal, d$var_a)
  deff <- 1 + (T_p - 1) * mean(c(i0, i1))

  N <- 280
  eff_group <- (N / 2) * T_p / deff
  p_ref <- power.prop.test(n = eff_group, p1 = 0.30, p2 = 0.40)$power
  pw <- ape_power(d, n = N, claim = "detect", nsim = 800, seed = 51)
  expect_lt(abs(pw$power - p_ref), 0.05)
})

test_that("panel APE and clustered SE match marginaleffects + sandwich on fixed data", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  d <- panel_base(target = 0.08)
  set.seed(77)
  xx <- powerape:::draw_x_panel(d, 400)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  dat <- data.frame(y = y, treat = xx$X[, 2], size = xx$X[, 3],
                    sizebar = xx$X[, 4], id = xx$id)
  fit <- glm(y ~ treat + size + sizebar, binomial("probit"), dat)

  X <- cbind(1, dat$treat, dat$size, dat$sizebar)
  ft <- powerape:::fit_index_model(X, y, "probit")
  est <- powerape:::ape_est(ft$fit$coefficients, X, 2L, "probit", "binary")
  V <- powerape:::vcov_cluster(ft$fit, X, dat$id, "probit")
  se <- sqrt(drop(t(est$jac) %*% V %*% est$jac))

  mfx <- marginaleffects::avg_comparisons(fit, variables = list(treat = c(0, 1)),
                                          vcov = ~id)
  expect_equal(est$ape, mfx$estimate, tolerance = 1e-6)
  expect_equal(se, mfx$std.error, tolerance = 1e-3)
})

test_that("panel route is probit-only: logit refused with a teaching error", {
  ## The panel logit option was removed in 1.6.0 (pooled CRE logit is only a
  ## QMLE approximation of the ASF; the exact route is probit-specific). The
  ## refusal must teach, not just fail.
  expect_error(
    ape_dgp_panel("logit",
                  focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                  covariates = list(pa_var("size", "normal", icc = 0.6)),
                  n_periods = 4, rho = 0.5, cre_share = 0,
                  baseline = 0.30, signal = 0.10, n_int = 4e4, seed_int = 11L),
    "probit-only"
  )
  expect_error(
    ape_dgp_panel("logit",
                  focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                  n_periods = 4, rho = 0.5, baseline = 0.30),
    "quasi-ML approximation"
  )
  ## anything else non-probit fails plainly
  expect_error(
    ape_dgp_panel("cloglog",
                  focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                  n_periods = 4, rho = 0.5, baseline = 0.30),
    "must be \"probit\""
  )
})

test_that("gauss_hermite integrates normal moments exactly", {
  gh <- powerape:::gauss_hermite(20L)
  ## E[t^2] under exp(-t^2)/sqrt(pi) is 1/2; E[t^4] is 3/4
  expect_equal(sum(gh$weights * gh$nodes^2) / sqrt(pi), 0.5, tolerance = 1e-12)
  expect_equal(sum(gh$weights * gh$nodes^4) / sqrt(pi), 0.75, tolerance = 1e-12)
})

## ---- panel AIE (binary x binary) --------------------------------------------

panel_aie_base <- function(rho = 0.3, cre = 0.5, T_p = 4, pins = TRUE) {
  d <- ape_dgp_panel(
    focal = pa_var("treat", "binary", p = 0.5, icc = 1),
    moderator = pa_var("public", "binary", p = 0.4, icc = 0.7),
    covariates = list(pa_var("size", "normal", icc = 0.6)),
    n_periods = T_p, rho = rho, cre_share = cre,
    baseline = 0.30, signal = 0.10, n_int = 4e4, seed_int = 3L
  )
  if (pins) d <- set_aie(d, 0.06, main_focal = 0.10, main_moderator = 0.05)
  d
}

test_that("panel AIE construction, anchors, and guards", {
  d <- panel_aie_base()
  expect_equal(true_aie(d), 0.06, tolerance = 1e-7)
  ## the two main anchors are reproduced on the a-integrated ASF
  id <- powerape:::panel_integration(d)
  P <- d$mix_P
  p00 <- mean(P(d$beta0 + id$q))
  expect_equal(mean(P(d$beta0 + d$beta_focal + id$q)) - p00, 0.10,
               tolerance = 1e-7)
  expect_equal(mean(P(d$beta0 + d$beta_mod + id$q)) - p00, 0.05,
               tolerance = 1e-7)
  ## X layout: 1, d, m, dm, size, means(public, size), dm-bar
  set.seed(1)
  xx <- powerape:::draw_x_panel(d, 150)
  expect_identical(ncol(xx$X), 8L)
  expect_identical(length(xx$start), 8L)
  expect_true(all(xx$X[, 4] == xx$X[, 2] * xx$X[, 3]))
  ## unit-level pair: no Mundlak means, no interaction mean
  d2 <- ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                      moderator = pa_var("m", "binary", p = .4, icc = 1),
                      n_periods = 3, rho = .2, baseline = .3,
                      n_int = 2e4, seed_int = 4L)
  d2 <- set_aie(d2, 0.05, main_focal = 0.08, main_moderator = 0.04)
  set.seed(2)
  x2 <- powerape:::draw_x_panel(d2, 100)
  expect_identical(ncol(x2$X), 4L)
  ## guards
  expect_error(ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                             moderator = pa_var("w", "normal"),
                             n_periods = 3, baseline = .3),
               "binary focal x binary moderator")
  expect_error(set_ape(panel_aie_base(pins = FALSE), 0.1), "set_aie")
  expect_error(set_aie(panel_base(), 0.05, 0.1, 0.05), "no moderator")
  expect_error(set_aie(panel_aie_base(pins = FALSE), 0.7,
                       main_focal = 0.10, main_moderator = 0.05),
               "counterfactual cell")
})

test_that("panel AIE and model-based SE match ginteff on a fixed panel", {
  skip_if_not_installed("ginteff")
  d <- panel_aie_base()
  set.seed(42)
  xx <- powerape:::draw_x_panel(d, 400)
  y <- rbinom(nrow(xx$X), 1, xx$pr)
  X <- xx$X
  ## one fitted model, two estimators: ginteff on the glm (binaries as
  ## factors in dydxs, Mundlak means as plain columns = the ASF
  ## convention), our aie_est on the same fit with the design reordered
  ## into the engine's (1, d, m, dm, ...) layout -- isolates the estimand
  ## correspondence from fit-convergence noise
  df <- data.frame(y = y, d = factor(X[, 2]), m = factor(X[, 3]),
                   size = X[, 5], mbar = X[, 6], sbar = X[, 7],
                   dmbar = X[, 8])
  fit <- glm(y ~ d * m + size + mbar + sbar + dmbar, binomial("probit"),
             data = df)
  MM <- model.matrix(fit)
  ord <- c(1L, 2L, 3L, 8L, 4L, 5L, 6L, 7L)   # d:m sits last in the glm matrix
  e <- powerape:::aie_est(coef(fit)[ord], MM[, ord], "probit",
                          "binary", "binary")
  se_model <- sqrt(drop(t(e$jac) %*% vcov(fit)[ord, ord] %*% e$jac))
  g <- ginteff::ginteff(fit, dydxs = c("d", "m"))
  expect_equal(e$ape, unname(g$aie), tolerance = 1e-6)
  expect_equal(se_model, unname(g$se), tolerance = 1e-6)
})

test_that("panel AIE reduces to the cross-sectional AIE at rho = 0 with unit-level pairs", {
  dp <- ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                      moderator = pa_var("m", "binary", p = .4, icc = 1),
                      n_periods = 3, rho = 0, baseline = .30,
                      n_int = 4e4, seed_int = 6L)
  dp <- set_aie(dp, 0.06, main_focal = 0.10, main_moderator = 0.05)
  dc <- ape_dgp(focal = pa_var("t", "binary", p = .5),
                moderator = pa_var("m", "binary", p = .4),
                baseline = .30)
  dc <- set_aie(dc, 0.06, main_focal = 0.10, main_moderator = 0.05)
  ## with q = 0 both inversions solve the same scalar equations
  expect_equal(dp$beta_focal, dc$beta_focal, tolerance = 1e-6)
  expect_equal(dp$beta_mod, dc$beta_mod, tolerance = 1e-6)
  expect_equal(dp$beta_int, dc$beta_int, tolerance = 1e-6)
  expect_equal(true_aie(dp), true_aie(dc), tolerance = 1e-7)
  ## power: 300 units x 3 periods vs 900 independent observations (cluster
  ## SEs and design granularity differ, so equality is approximate)
  pp <- ape_power(dp, n = 300, claim = "detect", nsim = 800, seed = 7)
  pc <- ape_power(dc, n = 900, claim = "detect", nsim = 800, seed = 8)
  expect_lt(abs(pp$power - pc$power),
            4 * sqrt(pp$mcse^2 + pc$mcse^2) + 0.01)
})

test_that("panel AIE estimator recovers the target with calibrated clustered SEs", {
  d <- panel_aie_base()
  sims <- powerape:::sim_ci(d, 400, 300, 0.95, seed = 51)
  expect_gt(mean(sims$ok), 0.99)
  est <- (sims$l + sims$u)[sims$ok] / 2
  expect_equal(mean(est), 0.06, tolerance = 0.012)
  ratio <- mean(sims$se[sims$ok]) / sd(est)
  expect_lt(abs(ratio - 1), 0.25)
})

test_that("clustered CIs cover the true panel AIE at nominal rate", {
  d <- panel_aie_base()
  s <- powerape:::sim_ci(d, 300, 1200, 0.95, seed = 61)
  expect_gt(mean(s$ok), 0.99)
  cov <- mean((s$l[s$ok] <= 0.06) & (s$u[s$ok] >= 0.06))
  expect_gt(cov, 0.925)
  expect_lt(cov, 0.975)
})

test_that("panel AIE statement and print compose", {
  d <- panel_aie_base()
  pw <- ape_power(d, n = 200, claim = "detect", nsim = 60, seed = 71)
  txt <- unclass(power_statement(pw))
  expect_match(txt, "average interaction effect")
  expect_match(txt, "correlated random effects")
  expect_match(txt, "moderator public")
  expect_match(txt, "units observed over")
  expect_output(print(d), "moderator: public")
})
