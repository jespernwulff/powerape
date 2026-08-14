# Minimum detectable effect: inverse-consistency anchors, claims, routes,
# robustness mode, guards.

test_that("MDE inverts the two-proportion anchor (exact enumeration)", {
  d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  m <- ape_mde(d, n = 712, claim = "detect", nsim = 600, seed = 11)
  ## ape_n's four-tool anchor: n = 712 for target .10 at 80% -- so the MDE
  ## at 712 must come back near .10, and the EXACT power of the decision
  ## rule at the returned MDE must sit at the goal (conservative side ok)
  expect_lt(abs(m$mde - 0.10), 0.012)
  ex <- exact_power_sat(712, 0.5, 0.30, 0.30 + m$mde, 0.95, "detect")
  expect_gt(ex, 0.775)
  expect_lt(ex, 0.865)
  ## analytic inverse from power.prop.test agrees
  p2 <- power.prop.test(n = 356, p1 = 0.30, power = 0.80)$p2
  expect_lt(abs(m$mde - (p2 - 0.30)), 0.012)
})

test_that("minimum-claim MDE exceeds the SESOI and the detect MDE", {
  d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  md <- ape_mde(d, n = 712, claim = "detect", nsim = 400, seed = 21,
                confirm = FALSE)
  mm <- ape_mde(d, n = 712, claim = "minimum", sesoi = 0.05, nsim = 400,
                seed = 22, confirm = FALSE)
  expect_gt(mm$mde, 0.05)
  expect_gt(mm$mde, md$mde)
})

test_that("equivalence MDE returns the smallest establishable margin", {
  d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  d0 <- set_ape(d, 0)
  m <- ape_mde(d0, n = 1400, claim = "equivalence", nsim = 500, seed = 31)
  expect_true(m$is_margin)
  expect_gt(m$mde, 0)
  expect_gt(m$power, 0.76)
  ## the four-tool battery says bounds +/-.07 give ~63% at n = 1400, so the
  ## 80% margin must be wider than .07
  expect_gt(m$mde, 0.07)
})

test_that("panel MDE inverts the Design D anchor (units, not rows)", {
  dp <- ape_dgp_panel(focal = pa_var("adopt", "binary", p = 0.5, icc = 1),
                      covariates = list(pa_var("size", "normal", icc = 0.6)),
                      n_periods = 4, rho = 0.30, cre_share = 0.50,
                      baseline = 0.30, signal = 0.10, n_int = 4e4)
  ## 719 units reach ~80% minimum-claim power at target .10 (Section 5,
  ## Design D), so the MDE at 719 units should return ~.10
  m <- ape_mde(dp, n = 719, claim = "minimum", sesoi = 0.05, nsim = 500,
               seed = 41)
  expect_lt(abs(m$mde - 0.10), 0.015)
  expect_output(print(m), "719 units")
})

test_that("AIE MDE reuses the pinned anchors and confirms", {
  da <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                moderator = pa_var("female", "binary", p = 0.55),
                baseline = 0.30)
  da <- set_aie(da, 0.08, main_focal = 0.10, main_moderator = 0.05)
  m <- ape_mde(da, n = 1500, claim = "detect", nsim = 400, seed = 51)
  expect_identical(m$estimand, "aie")
  expect_gt(m$mde, 0)
  expect_gt(m$power, m$goal - 0.03)
  txt <- unclass(power_statement(m))
  expect_match(txt, "minimum detectable AIE")
  expect_match(txt, "average interaction effect")
})

test_that("ape_robust mode = 'mde' sweeps the MDE and flags the worst", {
  d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  rb <- ape_robust(d, n = 712, claim = "detect", mode = "mde", power = 0.80,
                   vary = list(baseline = c(0.20, 0.40)), grid_points = 3,
                   nsim = 300, seed = 61)
  expect_identical(nrow(rb$scenarios), 3L)
  expect_true(all(is.finite(rb$scenarios$mde)))
  expect_equal(rb$worst$mde, max(rb$scenarios$mde))
  expect_output(print(rb), "MDE range")
})

test_that("MDE guards fire", {
  da <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                moderator = pa_var("female", "binary", p = 0.55),
                baseline = 0.30)
  expect_error(ape_mde(da, n = 1000, claim = "detect", nsim = 100),
               "main-effect anchors")
  d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  expect_error(ape_mde(d, n = 1000, claim = "equivalence", nsim = 100),
               "pin the DGP first")
  dh <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.85)
  expect_error(ape_mde(dh, n = 40, claim = "detect", nsim = 200, seed = 71),
               "feasible range")
  expect_error(ape_mde(d, n = 1000, claim = "minimum", nsim = 100),
               "sesoi")
})
