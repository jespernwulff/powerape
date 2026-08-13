# DESIGN.md sections 2.1 and 7.2: the internal AIE must reproduce ginteff
# exactly on fixed datasets. Binary variables enter ginteff as factors in
# `dydxs` (the analog of Stata's i. syntax): discrete change from the base
# level. Numeric `firstdiff` would instead mean "+1 unit from the observed
# value" -- a different estimand.

pa_se <- function(est, V) sqrt(drop(t(est$jac) %*% V %*% est$jac))

test_that("AIE matches ginteff exactly: binary x binary (probit)", {
  skip_if_not_installed("ginteff")
  set.seed(20)
  n <- 1200
  d <- rbinom(n, 1, 0.5); m <- rbinom(n, 1, 0.45); x <- rnorm(n)
  y <- rbinom(n, 1, pnorm(-0.4 + 0.35 * d + 0.25 * m + 0.3 * d * m + 0.4 * x))
  fit <- glm(y ~ d * m + x, binomial("probit"),
             data = data.frame(y, d = factor(d), m = factor(m), x))

  X <- cbind(1, d, m, d * m, x)
  ft <- powerape:::fit_index_model(X, y, "probit")
  e <- powerape:::aie_est(ft$fit$coefficients, X, "probit", "binary", "binary")
  se <- pa_se(e, powerape:::vcov_from_glmfit(ft$fit))

  g <- ginteff::ginteff(fit, dydxs = c("d", "m"))
  expect_equal(e$ape, unname(g$aie), tolerance = 1e-5)
  expect_equal(se, unname(g$se), tolerance = 1e-4)
})

test_that("AIE matches ginteff exactly: continuous x binary (probit)", {
  skip_if_not_installed("ginteff")
  set.seed(21)
  n <- 1200
  m <- rbinom(n, 1, 0.45); w <- rnorm(n); x <- rnorm(n)
  y <- rbinom(n, 1, pnorm(-0.4 + 0.3 * w + 0.25 * m + 0.2 * w * m + 0.4 * x))
  fit <- glm(y ~ w * m + x, binomial("probit"),
             data = data.frame(y, w, m = factor(m), x))

  X <- cbind(1, w, m, w * m, x)
  ft <- powerape:::fit_index_model(X, y, "probit")
  e <- powerape:::aie_est(ft$fit$coefficients, X, "probit", "normal", "binary")
  se <- pa_se(e, powerape:::vcov_from_glmfit(ft$fit))

  g <- ginteff::ginteff(fit, dydxs = c("w", "m"))
  expect_equal(e$ape, unname(g$aie), tolerance = 1e-5)
  expect_equal(se, unname(g$se), tolerance = 1e-4)
})

test_that("AIE matches ginteff exactly: binary x continuous (logit)", {
  skip_if_not_installed("ginteff")
  set.seed(23)
  n <- 1200
  d <- rbinom(n, 1, 0.5); w <- rnorm(n); x <- rnorm(n)
  y <- rbinom(n, 1, plogis(-0.4 + 0.5 * d + 0.3 * w + 0.4 * d * w + 0.5 * x))
  fit <- glm(y ~ d * w + x, binomial("logit"),
             data = data.frame(y, d = factor(d), w, x))

  X <- cbind(1, d, w, d * w, x)
  ft <- powerape:::fit_index_model(X, y, "logit")
  e <- powerape:::aie_est(ft$fit$coefficients, X, "logit", "binary", "normal")
  se <- pa_se(e, powerape:::vcov_from_glmfit(ft$fit))

  g <- ginteff::ginteff(fit, dydxs = c("d", "w"))
  expect_equal(e$ape, unname(g$aie), tolerance = 1e-5)
  expect_equal(se, unname(g$se), tolerance = 1e-4)
})

test_that("AIE matches ginteff exactly: continuous x continuous (probit)", {
  skip_if_not_installed("ginteff")
  set.seed(22)
  n <- 1200
  w <- rnorm(n); v <- rnorm(n)
  y <- rbinom(n, 1, pnorm(-0.3 + 0.3 * w + 0.2 * v + 0.15 * w * v))
  fit <- glm(y ~ w * v, binomial("probit"), data = data.frame(y, w, v))

  X <- cbind(1, w, v, w * v)
  ft <- powerape:::fit_index_model(X, y, "probit")
  e <- powerape:::aie_est(ft$fit$coefficients, X, "probit", "normal", "normal")
  se <- pa_se(e, powerape:::vcov_from_glmfit(ft$fit))

  g <- ginteff::ginteff(fit, dydxs = c("w", "v"))
  expect_equal(e$ape, unname(g$aie), tolerance = 1e-5)
  expect_equal(se, unname(g$se), tolerance = 1e-3)
})
