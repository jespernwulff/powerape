# AIE: inversion, truth evaluation, engine.

test_that("binary x binary AIE inversion matches the closed form (no covariates)", {
  d <- ape_dgp("probit",
               focal = pa_var("treat", "binary", p = 0.5),
               moderator = pa_var("female", "binary", p = 0.55),
               baseline = 0.30)
  d <- set_aie(d, target = 0.04, main_focal = 0.10, main_moderator = 0.05)
  ## implied counterfactual cell rates: p00 = .30, p10 = .40, p01 = .35, p11 = .49
  expect_equal(d$beta_focal, qnorm(0.40) - qnorm(0.30), tolerance = 1e-6)
  expect_equal(d$beta_mod, qnorm(0.35) - qnorm(0.30), tolerance = 1e-6)
  expect_equal(d$beta_int,
               qnorm(0.49) - qnorm(0.40) - qnorm(0.35) + qnorm(0.30),
               tolerance = 1e-5)
  expect_equal(true_aie(d), 0.04, tolerance = 1e-7)
  expect_identical(d$estimand, "aie")
})

test_that("AIE guards and feasibility errors fire", {
  d <- ape_dgp("probit",
               focal = pa_var("treat", "binary", p = 0.5),
               moderator = pa_var("female", "binary", p = 0.55),
               baseline = 0.30)
  expect_error(set_ape(d, 0.05), "set_aie")
  expect_error(set_aie(d, 0.04), "main-effect anchors")
  expect_error(set_aie(d, 0.60, main_focal = 0.10, main_moderator = 0.05),
               "infeasible")
  expect_error(set_aie(d, 0.04, main_focal = 0.80, main_moderator = 0.05),
               "infeasible")

  d_plain <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                     baseline = 0.30)
  expect_error(set_aie(d_plain, 0.04, main_focal = 0.1, main_moderator = 0.1),
               "no moderator")
})

test_that("continuous-focal x binary-moderator AIE pins and verifies on a fresh draw", {
  d <- ape_dgp("probit",
               focal = pa_var("x", "normal", mean = 0, sd = 1),
               moderator = pa_var("m", "binary", p = 0.45),
               covariates = list(pa_var("z", "normal")),
               correlation = 0.2, baseline = 0.35, signal = 0.15)
  d <- set_aie(d, target = 0.03, main_focal = 0.06, main_moderator = 0.05)
  expect_equal(true_aie(d), 0.03, tolerance = 1e-7)

  set.seed(77)
  xx <- powerape:::draw_x(d, 4e5)
  bt <- powerape:::beta_true(d)
  est <- powerape:::aie_est(bt, xx$X, "probit", "normal", "binary")
  expect_equal(est$ape, 0.03, tolerance = 4e-3)
})

test_that("continuous x continuous AIE pins and verifies on a fresh draw", {
  d <- ape_dgp("probit",
               focal = pa_var("x1", "normal", mean = 0, sd = 1),
               moderator = pa_var("x2", "normal", mean = 1, sd = 2),
               baseline = 0.40)
  d <- set_aie(d, target = 0.01, main_focal = 0.05, main_moderator = 0.03)
  expect_equal(true_aie(d), 0.01, tolerance = 1e-7)
  expect_error(set_aie(d, 2.0, main_focal = 0.05, main_moderator = 0.03),
               "infeasible")

  set.seed(78)
  xx <- powerape:::draw_x(d, 4e5)
  bt <- powerape:::beta_true(d)
  est <- powerape:::aie_est(bt, xx$X, "probit", "normal", "normal")
  expect_equal(est$ape, 0.01, tolerance = 4e-3)
})

test_that("AIE engine runs, SEs are calibrated, and the claims pipeline works", {
  d <- ape_dgp("probit",
               focal = pa_var("treat", "binary", p = 0.5),
               moderator = pa_var("female", "binary", p = 0.55),
               baseline = 0.30)
  d <- set_aie(d, target = 0.08, main_focal = 0.10, main_moderator = 0.05)

  sims <- powerape:::sim_ci(d, n = 3000, nsim = 300, conf = 0.95, seed = 55)
  expect_gt(mean(sims$ok), 0.99)
  est <- (sims$l + sims$u)[sims$ok] / 2
  ratio <- mean(sims$se[sims$ok]) / sd(est)
  expect_lt(abs(ratio - 1), 0.25)

  pw <- ape_power(d, n = 3000, claim = "minimum", sesoi = 0.02,
                  nsim = 400, seed = 56)
  expect_identical(pw$estimand, "aie")
  expect_true(pw$power > 0 && pw$power <= 1)
  expect_equal(sum(pw$outcomes), 1, tolerance = 1e-12)
})
