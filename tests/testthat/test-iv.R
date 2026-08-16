# IV route: construction, inversion, engine, reduction, coverage.

iv_base <- function(rho = 0.4, ivs = 0.2, pins = TRUE) {
  d <- ape_dgp_iv(
    focal = pa_var("training", "binary", p = 0.4),
    covariates = list(pa_var("size", "normal")),
    instruments = pa_var("subsidy", "binary", p = 0.5),
    endogeneity = rho, iv_strength = ivs,
    baseline = 0.30, signal = 0.10, n_int = 4e4, seed_int = 2L
  )
  if (pins) d <- set_ape(d, 0.10)
  d
}

test_that("IV construction calibrates baseline, prevalence, and strength", {
  d <- iv_base(pins = FALSE)
  id <- powerape:::iv_integration(d)
  expect_equal(mean(pnorm(d$beta0 + id$idxz)), 0.30, tolerance = 1e-6)
  ## focal prevalence reproduced on a fresh draw
  rd <- powerape:::iv_draw(d, 60000, seed = 9L)
  expect_lt(abs(mean(rd$d) - 0.40), 0.01)
  ## latent instrument share ~ iv_strength
  qz <- drop(rd$Zi %*% d$pi_z)
  expect_equal(var(qz) / (var(qz) + 1), 0.2, tolerance = 0.01)
  ## continuous focal: variance shares
  dc <- ape_dgp_iv(focal = pa_var("rd", "normal", mean = 1, sd = 2),
                   instruments = pa_var("z1", "normal"),
                   endogeneity = 0.3, iv_strength = 0.3,
                   baseline = 0.25, n_int = 4e4, seed_int = 3L)
  rc <- powerape:::iv_draw(dc, 60000, seed = 10L)
  expect_equal(mean(rc$d), 1, tolerance = 0.03)
  expect_equal(sd(rc$d), 2, tolerance = 0.03)
  expect_equal(var(drop(rc$Zi %*% dc$pi_z)) / var(rc$d), 0.3, tolerance = 0.01)
  ## guards
  expect_error(ape_dgp_iv(focal = pa_var("t", "binary", p = .4),
                          instruments = pa_var("z", "normal"),
                          endogeneity = 1.2, iv_strength = .2, baseline = .3),
               "endogeneity")
  expect_error(ape_dgp_iv(focal = pa_var("t", "binary", p = .4),
                          moderator = pa_var("w", "normal"),
                          instruments = pa_var("z", "normal"),
                          endogeneity = .3, iv_strength = .2, baseline = .3),
               "binary moderator")
})

test_that("IV inversion pins the ASF APE and AIE exactly", {
  d <- iv_base()
  expect_equal(true_ape(d), 0.10, tolerance = 1e-7)
  di <- ape_dgp_iv(
    focal = pa_var("training", "binary", p = 0.4),
    moderator = pa_var("public", "binary", p = 0.5),
    covariates = list(pa_var("size", "normal")),
    instruments = pa_var("subsidy", "binary", p = 0.5),
    endogeneity = 0.4, iv_strength = 0.25,
    baseline = 0.30, signal = 0.10, n_int = 4e4, seed_int = 4L
  )
  di <- set_aie(di, 0.06, main_focal = 0.10, main_moderator = 0.05)
  expect_equal(true_aie(di), 0.06, tolerance = 1e-7)
})

test_that("CF estimator recovers the true APE (consistency + SE calibration)", {
  ## continuous focal, linear first stage: the exactly identified case.
  ## Two-step CF carries a small finite-sample attenuation (~0.007 at
  ## n = 900 here), so consistency is asserted at a larger n.
  dc <- ape_dgp_iv(focal = pa_var("rd", "normal", mean = 0, sd = 1),
                   instruments = pa_var("z1", "normal"),
                   endogeneity = 0.5, iv_strength = 0.3,
                   baseline = 0.30, n_int = 4e4, seed_int = 5L)
  dc <- set_ape(dc, 0.08)
  s <- powerape:::sim_ci(dc, 2500, 250, 0.95, seed = 41)
  expect_gt(mean(s$ok), 0.99)
  est <- (s$l + s$u)[s$ok] / 2
  expect_lt(abs(mean(est) - 0.08), 0.008)
  ratio <- mean(s$se[s$ok]) / sd(est)
  expect_lt(abs(ratio - 1), 0.25)
})

test_that("binary-endogenous CF plim gap is small (documented approximation)", {
  d <- iv_base()
  s <- powerape:::sim_ci(d, 1500, 300, 0.95, seed = 51)
  expect_gt(mean(s$ok), 0.98)
  est <- (s$l + s$u)[s$ok] / 2
  ## the generalized-residual CF for a binary endogenous variable is the
  ## standard Wooldridge approximation; measured gap ~ -0.010 here
  ## (about 10% of the effect at endogeneity .4, strength .2)
  expect_lt(abs(mean(est) - 0.10), 0.02)
})

test_that("the IV price: CF uses only instrument variation (SE-ratio anchor)", {
  ## At endogeneity = 0 the CF route still identifies the effect from the
  ## instrumented variation alone, so its SE is inflated by about
  ## sqrt(1 / iv_strength) relative to the exogenous estimator -- the
  ## familiar 2SLS variance logic. Anchor the mechanism quantitatively.
  ivs <- 0.3
  dc0 <- ape_dgp_iv(focal = pa_var("rd", "normal", mean = 0, sd = 1),
                    instruments = pa_var("z1", "normal"),
                    endogeneity = 0, iv_strength = ivs,
                    baseline = 0.30, n_int = 4e4, seed_int = 7L)
  dc0 <- set_ape(dc0, 0.08)
  ds0 <- ape_dgp(focal = pa_var("rd", "normal", mean = 0, sd = 1),
                 baseline = 0.30)
  ds0 <- set_ape(ds0, 0.08)
  expect_equal(true_ape(dc0), true_ape(ds0), tolerance = 2e-3)
  s_iv <- powerape:::sim_ci(dc0, 1500, 250, 0.95, seed = 61)
  s_ex <- powerape:::sim_ci(ds0, 1500, 250, 0.95, seed = 62)
  ## consistency also holds at endogeneity 0
  expect_lt(abs(mean((s_iv$l + s_iv$u)[s_iv$ok] / 2) - 0.08), 0.01)
  ratio <- mean(s_iv$se[s_iv$ok]) / mean(s_ex$se[s_ex$ok])
  expect_gt(ratio, 0.80 * sqrt(1 / ivs))
  expect_lt(ratio, 1.25 * sqrt(1 / ivs))
})

test_that("IV CIs cover the true APE at nominal rate", {
  dc <- ape_dgp_iv(focal = pa_var("rd", "normal", mean = 0, sd = 1),
                   instruments = pa_var("z1", "normal"),
                   endogeneity = 0.4, iv_strength = 0.3,
                   baseline = 0.30, n_int = 4e4, seed_int = 6L)
  dc <- set_ape(dc, 0.08)
  s <- powerape:::sim_ci(dc, 800, 1000, 0.95, seed = 71)
  expect_gt(mean(s$ok), 0.98)
  cov <- mean((s$l[s$ok] <= 0.08) & (s$u[s$ok] >= 0.08))
  expect_gt(cov, 0.92)
  expect_lt(cov, 0.98)
})

test_that("IV-AIE estimator is consistent with calibrated SEs (stochastic)", {
  di <- ape_dgp_iv(
    focal = pa_var("t", "binary", p = .4),
    moderator = pa_var("m", "binary", p = .5),
    covariates = list(pa_var("x", "normal")),
    instruments = pa_var("z", "binary", p = .5),
    endogeneity = .4, iv_strength = .25,
    baseline = .30, signal = .10, n_int = 4e4, seed_int = 4L
  )
  di <- set_aie(di, .06, main_focal = .10, main_moderator = .05)
  s <- powerape:::sim_ci(di, 2000, 300, .95, seed = 172)
  expect_gt(mean(s$ok), 0.98)
  est <- (s$l + s$u)[s$ok] / 2
  expect_lt(abs(mean(est) - 0.06), 0.015)
  expect_lt(abs(mean(s$se[s$ok]) / sd(est) - 1), 0.25)
})

test_that("se option: robust works on standard routes, guarded elsewhere", {
  ds <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  ds <- set_ape(ds, 0.10)
  pr <- ape_power(ds, n = 400, claim = "detect", nsim = 100, seed = 81,
                  se = "robust")
  expect_s3_class(pr, "powerape_power")
  expect_identical(pr$se, "robust")
  expect_match(unclass(power_statement(pr)), "heteroskedasticity-robust")
  d <- iv_base()
  expect_warning(ape_power(d, n = 400, claim = "detect", nsim = 30, seed = 82,
                           se = "robust"),
                 "fixed by the route")
})

test_that("IV print and statement compose", {
  d <- iv_base()
  expect_output(print(d), "endogenous focal")
  pw <- ape_power(d, n = 600, claim = "detect", nsim = 60, seed = 91)
  txt <- unclass(power_statement(pw))
  expect_match(txt, "control-function probit")
  expect_match(txt, "instrumented by subsidy")
  expect_match(txt, "stacked method-of-moments")
})
