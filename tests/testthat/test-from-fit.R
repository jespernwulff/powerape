make_fit <- function(n = 1500, seed = 42, link = "probit") {
  set.seed(seed)
  age <- rnorm(n, 45, 12)
  female <- rbinom(n, 1, 0.55)
  treat <- rbinom(n, 1, plogis(0.5 * scale(age)[, 1]))
  G <- if (link == "probit") pnorm else plogis
  y <- rbinom(n, 1, G(-0.8 + 0.35 * treat + 0.015 * age + 0.2 * female))
  glm(y ~ treat + age + female, family = binomial(link),
      data = data.frame(y, treat, age, female))
}

test_that("from-fit route lifts coefficients and calibrates exactly", {
  fit <- make_fit()
  d <- ape_dgp_from_fit(fit, focal = "treat")
  expect_identical(d$route, "empirical")
  expect_identical(d$builder, "from_fit")
  expect_equal(d$gamma, unname(coef(fit)[c("age", "female")]), tolerance = 1e-12)
  expect_equal(d$beta0, unname(coef(fit)[["(Intercept)"]]), tolerance = 1e-12)

  ## implied baseline = counterfactual mean with treat = 0 over pilot rows
  X <- model.matrix(fit)
  man <- mean(pnorm(coef(fit)[["(Intercept)"]] +
                      drop(X[, c("age", "female")] %*% coef(fit)[c("age", "female")])))
  expect_equal(d$baseline, man, tolerance = 1e-12)

  d <- set_ape(d, 0.05)
  expect_equal(true_ape(d), 0.05, tolerance = 1e-7)
  pw <- ape_power(d, n = 600, claim = "detect", nsim = 200, seed = 3)
  expect_true(pw$power > 0 && pw$power <= 1)
  expect_output(print(d), "pilot model")
})

test_that("baseline override recalibrates the intercept", {
  fit <- make_fit()
  d <- ape_dgp_from_fit(fit, focal = "treat", baseline = 0.25)
  id <- powerape:::integration_draw(d)
  expect_equal(mean(pnorm(d$beta0 + id$idxz)), 0.25, tolerance = 1e-6)
  expect_equal(d$baseline, 0.25)
})

test_that("continuous focal from a logit fit works", {
  fit <- make_fit(link = "logit")
  d <- ape_dgp_from_fit(fit, focal = "age")
  expect_identical(d$focal$type, "normal")
  expect_identical(d$model, "logit")
  d <- set_ape(d, 0.004)
  expect_equal(true_ape(d), 0.004, tolerance = 1e-8)
})

test_that("from-fit validation errors fire", {
  fit <- make_fit()
  expect_error(ape_dgp_from_fit(fit, focal = "nope"), "not a column")
  expect_error(ape_dgp_from_fit(lm(mpg ~ wt, mtcars), focal = "wt"), "binomial")

  set.seed(1)
  dd <- data.frame(y = rbinom(800, 1, 0.4), d = rbinom(800, 1, 0.5),
                   x = rnorm(800))
  fit2 <- glm(y ~ d * x, binomial("probit"), dd)
  expect_error(ape_dgp_from_fit(fit2, focal = "d"), "nteraction")

  fit3 <- glm(y ~ d + x, binomial("cauchit"), dd)
  expect_error(ape_dgp_from_fit(fit3, focal = "d"), "probit or logit")
})

test_that("ape_robust varies baseline for from-fit DGPs and blocks signal", {
  fit <- make_fit()
  d <- set_ape(ape_dgp_from_fit(fit, focal = "treat"), 0.05)
  rb <- ape_robust(d, n = 400, claim = "detect",
                   vary = list(baseline = c(0.20, 0.30, 0.40)),
                   nsim = 150, seed = 4, nmax = FALSE)
  expect_equal(nrow(rb$scenarios), 3L)
  expect_true(all(rb$scenarios$implied_effect == 0.05))
  expect_error(ape_robust(d, 400, claim = "detect",
                          vary = list(signal = c(0, 0.2)), nmax = FALSE),
               "Cannot vary")
})
