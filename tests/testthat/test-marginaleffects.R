# DESIGN.md section 7.2: APE point estimates and delta-method SEs must match
# marginaleffects on fixed datasets.

test_that("APE and SE match marginaleffects (probit, binary focal)", {
  skip_if_not_installed("marginaleffects")
  set.seed(10)
  n <- 1500
  dat <- data.frame(d = rbinom(n, 1, 0.5), x = rnorm(n))
  dat$y <- rbinom(n, 1, pnorm(-0.5 + 0.4 * dat$d + 0.6 * dat$x))
  fit <- glm(y ~ d + x, family = binomial("probit"), data = dat)
  X <- model.matrix(fit)

  est <- powerape:::ape_est(coef(fit), X, 2L, "probit", "binary")
  se <- sqrt(drop(t(est$jac) %*% vcov(fit) %*% est$jac))

  mfx <- marginaleffects::avg_comparisons(fit, variables = list(d = c(0, 1)))
  expect_equal(est$ape, mfx$estimate, tolerance = 1e-6)
  expect_equal(se, mfx$std.error, tolerance = 1e-4)
})

test_that("APE and SE match marginaleffects (probit and logit, continuous focal)", {
  skip_if_not_installed("marginaleffects")
  set.seed(11)
  n <- 1500
  dat <- data.frame(w = rnorm(n, 1, 2), x = rnorm(n))
  for (fam in c("probit", "logit")) {
    G <- if (fam == "probit") pnorm else plogis
    dat$y <- rbinom(n, 1, G(-0.4 + 0.3 * dat$w - 0.5 * dat$x))
    fit <- glm(y ~ w + x, family = binomial(fam), data = dat)
    X <- model.matrix(fit)

    est <- powerape:::ape_est(coef(fit), X, 2L, fam, "normal")
    se <- sqrt(drop(t(est$jac) %*% vcov(fit) %*% est$jac))

    mfx <- marginaleffects::avg_slopes(fit, variables = "w")
    expect_equal(est$ape, mfx$estimate, tolerance = 1e-4)
    expect_equal(se, mfx$std.error, tolerance = 1e-3)
  }
})
