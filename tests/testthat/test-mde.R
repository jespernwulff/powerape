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

test_that("IV-route MDE inverts the Design E anchor", {
  de <- ape_dgp_iv(
    focal       = pa_var("training", "binary", p = 0.4),
    covariates  = list(pa_var("size", "normal")),
    instruments = pa_var("subsidy", "binary", p = 0.5),
    endogeneity = 0.4, iv_strength = 0.2,
    baseline    = 0.30, signal = 0.10, n_int = 4e4, seed_int = 2L
  )
  ## Section 5 Design E: n = 6,433 reaches ~80% detect power at target .10,
  ## so the MDE at that n must come back near .10 (verified at high nsim
  ## with confirmation: 0.104)
  m <- ape_mde(de, n = 6433, claim = "detect", nsim = 300, seed = 81,
               confirm = FALSE)
  expect_lt(abs(m$mde - 0.10), 0.02)
})

test_that("ape_mde follows the pinned direction and obeys the outcome-flip identity", {
  ## y -> 1 - y maps the world (baseline b, gamma, APE -t) exactly onto
  ## (1 - b, -gamma, +t) (probit symmetry; rbinom's inversion couples the
  ## draws), so the decrease MDE of one world is minus the increase MDE of
  ## the other -- a relation 1.10.0 violated by always searching increases
  dA <- ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                covariates = list(pa_var("z", "normal")), baseline = 0.08,
                gamma = 0.3)
  dB <- ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                covariates = list(pa_var("z", "normal")), baseline = 0.92,
                gamma = -0.3)
  mA <- ape_mde(set_ape(dA, -0.03), n = 3000, claim = "detect", nsim = 200,
                seed = 5, confirm = FALSE)
  mB <- ape_mde(set_ape(dB, 0.03), n = 3000, claim = "detect", nsim = 200,
                seed = 5, confirm = FALSE)
  expect_identical(mA$direction, "negative")
  expect_identical(mB$direction, "positive")
  expect_lt(mA$mde, 0)
  expect_equal(mA$mde, -mB$mde, tolerance = 1e-8)
  expect_equal(mA$dgp$target_est, mA$mde)
  ## away from a .5 baseline the direction matters: at .08 a decrease is
  ## easier to detect than an increase (two-proportion analytic: .0256 vs
  ## .0300 at n = 3000)
  mInc <- ape_mde(dA, n = 3000, claim = "detect", nsim = 200, seed = 5,
                  confirm = FALSE)
  expect_identical(mInc$direction, "positive")
  expect_gt(mInc$mde, 1.08 * abs(mA$mde))
  expect_lt(abs(abs(mA$mde) - 0.0256), 0.0015)
  expect_output(print(mA), "decreases")
  expect_match(unclass(power_statement(mA)), "a decrease of")
  ## an explicit direction on an unpinned DGP
  mD <- ape_mde(dA, n = 3000, claim = "detect", nsim = 200, seed = 5,
                confirm = FALSE, direction = "negative")
  expect_lt(mD$mde, 0)
  expect_error(ape_mde(dA, n = 3000, claim = "detect", nsim = 50,
                       direction = "down"), "arg")
})

test_that("ape_mde re-measures the SE at its proposal: no first-step bias at large n", {
  ## at n = 10,000 the MDE (~.0175) is far below the reference effect of a
  ## .10 baseline (.05), where the APE's SE is ~6% larger; an SE taken there
  ## overshot the MDE by ~6% in 1.10.0, and a first candidate accepted
  ## within the search tolerance is never trimmed afterwards
  d <- ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5), baseline = 0.10)
  cf <- uniroot(function(dd) dd - (qnorm(0.975) + qnorm(0.8)) *
                  sqrt(0.1 * 0.9 / 5000 + (0.1 + dd) * (0.9 - dd) / 5000),
                c(0.001, 0.1))$root
  m <- ape_mde(d, n = 10000, claim = "detect", nsim = 200, seed = 7,
               confirm = FALSE)
  expect_lt(abs(m$mde / cf - 1), 0.025)
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
