test_that("saturated-case APE equals the difference in proportions and its SE the unpooled formula", {
  set.seed(1)
  n <- 4000
  d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.30 + 0.10 * d)
  X <- cbind(1, d)

  ft <- powerape:::fit_index_model(X, y, "probit")
  expect_true(ft$ok)
  est <- powerape:::ape_binary_est(ft$fit$coefficients, X, 2L, "probit")

  p1 <- mean(y[d == 1])
  p0 <- mean(y[d == 0])
  expect_equal(est$ape, p1 - p0, tolerance = 1e-6)

  V <- powerape:::vcov_from_glmfit(ft$fit)
  se <- sqrt(drop(t(est$jac) %*% V %*% est$jac))
  n1 <- sum(d == 1)
  n0 <- n - n1
  se_ref <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  expect_equal(se, se_ref, tolerance = 1e-5)
})

test_that("the same equivalences hold for logit", {
  set.seed(2)
  n <- 4000
  d <- rbinom(n, 1, 0.5)
  y <- rbinom(n, 1, 0.25 + 0.15 * d)
  X <- cbind(1, d)

  ft <- powerape:::fit_index_model(X, y, "logit")
  expect_true(ft$ok)
  est <- powerape:::ape_binary_est(ft$fit$coefficients, X, 2L, "logit")

  p1 <- mean(y[d == 1])
  p0 <- mean(y[d == 0])
  expect_equal(est$ape, p1 - p0, tolerance = 1e-6)

  V <- powerape:::vcov_from_glmfit(ft$fit)
  se <- sqrt(drop(t(est$jac) %*% V %*% est$jac))
  n1 <- sum(d == 1)
  n0 <- n - n1
  se_ref <- sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  expect_equal(se, se_ref, tolerance = 1e-5)
})
