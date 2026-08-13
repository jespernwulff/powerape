make_d <- function(target) {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               baseline = 0.30)
  set_ape(d, target)
}

test_that("coherence guards fire", {
  d5 <- make_d(0.05)
  expect_error(ape_power(d5, n = 500, claim = "minimum", sesoi = 0.05),
               "does not exceed")
  expect_error(ape_power(d5, n = 500, claim = "minimum", sesoi = 0.10),
               "does not exceed")
  expect_error(ape_power(d5, n = 500, claim = "minimum"),
               "needs a single positive")
  expect_error(ape_power(d5, n = 500, claim = "equivalence", sesoi = 0.05),
               "Equivalence power requires")
  d_unset <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                     baseline = 0.30)
  expect_error(ape_power(d_unset, n = 500, claim = "detect"),
               "set_ape")
})

test_that("detection power at target 0 equals the directional test size", {
  d0 <- make_d(0)
  pw <- ape_power(d0, n = 600, claim = "detect", nsim = 3000, seed = 7)
  expect_lt(abs(pw$power - 0.025), 0.012)
})

test_that("outcome distribution is a proper partition", {
  d10 <- make_d(0.10)
  pw <- ape_power(d10, n = 400, claim = "minimum", sesoi = 0.03,
                  nsim = 500, seed = 3)
  expect_named(pw$outcomes,
               c("minimum", "detect_only", "inconclusive", "equivalence", "failed"))
  expect_equal(sum(pw$outcomes), 1, tolerance = 1e-12)
  expect_true(pw$power >= 0 && pw$power <= 1)
  expect_equal(pw$power, pw$outcomes[["minimum"]], tolerance = 1e-12)
})

test_that("equivalence power at truth 0 grows with n", {
  d0 <- make_d(0)
  p_small <- ape_power(d0, n = 300, claim = "equivalence", sesoi = 0.05,
                       nsim = 800, seed = 11)
  p_large <- ape_power(d0, n = 3000, claim = "equivalence", sesoi = 0.05,
                       nsim = 800, seed = 12)
  expect_gt(p_large$power, p_small$power)
  expect_gt(p_large$power, 0.5)
})

test_that("negative targets mirror cleanly", {
  dneg <- make_d(-0.10)
  pw <- ape_power(dneg, n = 800, claim = "minimum", sesoi = 0.03,
                  nsim = 500, seed = 5)
  expect_gt(pw$power, 0.2)
  expect_equal(sum(pw$outcomes), 1, tolerance = 1e-12)
})

test_that("ape_curve and ape_n run and land near the analytic anchor", {
  d10 <- make_d(0.10)
  cv <- ape_curve(d10, n = c(200, 500), claim = "detect", nsim = 300, seed = 21)
  expect_s3_class(cv, "powerape_curve")
  expect_equal(nrow(cv$results), 2L)

  an <- ape_n(d10, power = 0.80, claim = "detect", nsim = 800, seed = 31)
  ppt <- stats::power.prop.test(p1 = 0.30, p2 = 0.40, power = 0.80,
                                sig.level = 0.05)
  n_ref <- 2 * ceiling(ppt$n)
  expect_lt(abs(an$n - n_ref), 150)
  expect_lt(abs(an$power - 0.80), 0.05)
})
