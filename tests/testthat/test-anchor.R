# The primary correctness benchmark (DESIGN.md section 7.1): with a binary
# focal variable and no covariates, the APE is a difference in proportions,
# so simulated detection power must reproduce power.prop.test().

test_that("powerape reproduces power.prop.test in the two-proportion case", {
  ppt <- stats::power.prop.test(p1 = 0.30, p2 = 0.40, power = 0.80,
                                sig.level = 0.05)
  n_tot <- 2L * ceiling(ppt$n)

  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  d <- set_ape(d, target = 0.10)
  expect_equal(d$p1, 0.40, tolerance = 1e-8)

  pw <- ape_power(d, n = n_tot, claim = "detect", conf = 0.95,
                  nsim = 4000, seed = 42)
  expect_lt(abs(pw$power - 0.80), 0.035)
  expect_equal(pw$outcomes[["failed"]], 0)
})
