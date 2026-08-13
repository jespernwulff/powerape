test_that("no-covariate probit inversion matches the closed form", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  d <- set_ape(d, 0.10)
  expect_equal(d$beta_focal, qnorm(0.40) - qnorm(0.30), tolerance = 1e-6)
  expect_equal(d$p1, 0.40, tolerance = 1e-8)
  expect_equal(true_ape(d), 0.10, tolerance = 1e-8)
})

test_that("infeasible targets error informatively", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5), baseline = 0.70)
  expect_error(set_ape(d, 0.50), "infeasible")
  expect_error(set_ape(d, -0.75), "infeasible")
})

test_that("inversion with covariates hits the target and is monotone", {
  covs <- list(pa_var("age", "normal", mean = 45, sd = 12),
               pa_var("female", "binary", p = 0.55))
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               covariates = covs, correlation = 0.2, baseline = 0.30,
               signal = 0.15)
  d5 <- set_ape(d, 0.05)
  expect_equal(true_ape(d5), 0.05, tolerance = 1e-7)
  d10 <- set_ape(d, 0.10)
  expect_gt(d10$beta_focal, d5$beta_focal)

  ## independent check on a fresh draw (different RNG stream than seed_int)
  set.seed(99)
  xx <- powerape:::draw_x(d5, 4e5)
  bt <- powerape:::beta_true(d5)
  X1 <- xx$X; X1[, 2] <- 1
  X0 <- xx$X; X0[, 2] <- 0
  ape_mc <- mean(pnorm(drop(X1 %*% bt)) - pnorm(drop(X0 %*% bt)))
  expect_equal(ape_mc, 0.05, tolerance = 4e-3)
})

test_that("baseline and signal calibration hold with covariates", {
  covs <- list(pa_var("x1", "normal"), pa_var("x2", "binary", p = 0.3))
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               covariates = covs, baseline = 0.25, signal = 0.30)
  idxz <- powerape:::integration_draw(d)$idxz
  p0 <- mean(pnorm(d$beta0 + idxz))
  expect_equal(p0, 0.25, tolerance = 1e-6)
  expect_equal(d$signal_implied, 0.30, tolerance = 1e-6)
})

test_that("logit inversion works too", {
  d <- ape_dgp("logit", focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
  d <- set_ape(d, 0.10)
  expect_equal(true_ape(d), 0.10, tolerance = 1e-8)
})
