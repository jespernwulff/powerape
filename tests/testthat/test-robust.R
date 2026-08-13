test_that("ape_robust sweeps scenarios and pins the APE by default", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  d <- set_ape(d, 0.10)
  rb <- ape_robust(d, n = 500, claim = "detect",
                   vary = list(baseline = c(0.20, 0.40)), grid_points = 3,
                   nsim = 300, seed = 5, nmax = FALSE)
  expect_s3_class(rb, "powerape_robust")
  expect_equal(nrow(rb$scenarios), 3L)
  expect_true(all(rb$scenarios$implied_effect == 0.10))
  expect_true(all(rb$scenarios$power >= 0 & rb$scenarios$power <= 1))
  expect_equal(rb$worst$power, min(rb$scenarios$power))
  expect_named(rb$marginals, "baseline")
  expect_output(print(rb), "robustness sweep")
})

test_that("pin = coefficients lets the implied effect drift", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  d <- set_ape(d, 0.10)
  rb <- ape_robust(d, n = 500, claim = "detect",
                   vary = list(baseline = c(0.20, 0.40)), grid_points = 3,
                   nsim = 300, seed = 5, pin = "coefficients", nmax = FALSE)
  expect_false(all(abs(rb$scenarios$implied_effect - 0.10) < 1e-6))
  ## the base scenario (baseline = .30) reproduces the base effect exactly
  mid <- which(abs(rb$scenarios$baseline - 0.30) < 1e-9)
  expect_equal(rb$scenarios$implied_effect[mid], 0.10, tolerance = 1e-6)
})

test_that("ape_robust validates vary inputs", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  d <- set_ape(d, 0.10)
  expect_error(ape_robust(d, 500, claim = "detect",
                          vary = list(nope = c(1, 2)), nsim = 100,
                          nmax = FALSE),
               "Cannot vary")
  expect_error(ape_robust(d, 500, claim = "detect", vary = list(c(0.2, 0.4)),
                          nsim = 100, nmax = FALSE),
               "named")
})

test_that("infeasible scenarios are reported, not fatal", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  d <- set_ape(d, 0.10)
  ## baseline = 0.95 makes a +0.10 APE infeasible (ceiling 0.05)
  rb <- ape_robust(d, n = 400, claim = "detect",
                   vary = list(baseline = c(0.30, 0.60, 0.95)),
                   nsim = 200, seed = 7, nmax = FALSE)
  expect_equal(sum(is.na(rb$scenarios$power)), 1L)
  expect_match(rb$scenarios$note[is.na(rb$scenarios$power)], "infeasible")
})

test_that("n_max runs in the worst scenario", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  d <- set_ape(d, 0.10)
  rb <- ape_robust(d, n = 400, claim = "detect",
                   vary = list(baseline = c(0.30, 0.45, 0.50)),
                   nsim = 400, seed = 9, nmax = TRUE, nmax_power = 0.80)
  expect_false(is.null(rb$nmax))
  expect_gt(rb$nmax$n, 100)
  expect_lt(abs(rb$nmax$power - 0.80), 0.06)
})
