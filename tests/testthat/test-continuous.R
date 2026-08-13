# Continuous (normal) focal variable: average-derivative APE.

test_that("continuous-focal probit inversion matches the closed form (no covariates)", {
  ## For u ~ N(0, s2): E[phi(a + b*u)] = phi(a / sqrt(1 + b^2 s2)) / sqrt(1 + b^2 s2)
  d <- ape_dgp("probit", focal = pa_var("x", "normal", mean = 2, sd = 1.5),
               baseline = 0.30)
  d <- set_ape(d, 0.08)
  b <- d$beta_focal
  a <- d$beta0
  s2 <- 1.5^2
  delta_closed <- b * dnorm(a / sqrt(1 + b^2 * s2)) / sqrt(1 + b^2 * s2)
  expect_equal(delta_closed, 0.08, tolerance = 2e-3)
  expect_equal(true_ape(d), 0.08, tolerance = 1e-7)
})

test_that("continuous-focal plateau makes large targets infeasible", {
  ## sup_b APE = phi(0)/sd = 0.266 here, so 0.30 must error
  d <- ape_dgp("probit", focal = pa_var("x", "normal", mean = 2, sd = 1.5),
               baseline = 0.30)
  expect_error(set_ape(d, 0.30), "infeasible")
})

test_that("negative continuous targets mirror", {
  d <- ape_dgp("probit", focal = pa_var("x", "normal", mean = 0, sd = 1),
               baseline = 0.40)
  d <- set_ape(d, -0.05)
  expect_lt(d$beta_focal, 0)
  expect_equal(true_ape(d), -0.05, tolerance = 1e-7)
})

test_that("continuous-focal inversion works with correlated covariates", {
  d <- ape_dgp("probit", focal = pa_var("x", "normal", mean = 0, sd = 1),
               covariates = list(pa_var("z1", "normal"),
                                 pa_var("z2", "binary", p = 0.4)),
               correlation = 0.3, baseline = 0.35, signal = 0.20)
  d <- set_ape(d, 0.06)
  expect_equal(true_ape(d), 0.06, tolerance = 1e-7)

  ## independent fresh-draw check of the average derivative
  set.seed(123)
  xx <- powerape:::draw_x(d, 3e5)
  bt <- powerape:::beta_true(d)
  a <- drop(xx$X %*% bt)
  ape_mc <- bt[2] * mean(dnorm(a))
  expect_equal(ape_mc, 0.06, tolerance = 4e-3)
})

test_that("continuous-focal engine runs and its SEs are calibrated", {
  d <- ape_dgp("probit", focal = pa_var("x", "normal", mean = 0, sd = 1),
               baseline = 0.30)
  d <- set_ape(d, 0.06)
  sims <- powerape:::sim_ci(d, n = 2000, nsim = 300, conf = 0.95, seed = 44)
  expect_gt(mean(sims$ok), 0.99)
  est <- (sims$l + sims$u)[sims$ok] / 2
  ratio <- mean(sims$se[sims$ok]) / sd(est)
  expect_lt(abs(ratio - 1), 0.25)

  pw <- ape_power(d, n = 1500, claim = "detect", nsim = 300, seed = 45)
  expect_gt(pw$power, 0.9)
})
