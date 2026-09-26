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

test_that("the focal guard refuses every other route by which the focal enters", {
  ## each of these pilots would hold a function of the focal fixed while the
  ## focal changes, so the pinned effect would not be the model's APE
  set.seed(3)
  pd <- data.frame(dose = rnorm(900), age = rnorm(900),
                   female = rbinom(900, 1, 0.5), treatL = runif(900) < 0.5,
                   expo = runif(900, 1, 3))
  pd$treat <- as.integer(pd$treatL)
  pd$y <- rbinom(900, 1, pnorm(-0.5 + 0.3 * pd$dose + 0.2 * pd$age +
                                 0.2 * pd$treat))
  f <- function(fml) glm(fml, binomial("probit"), pd)
  expect_error(ape_dgp_from_fit(f(y ~ dose + I(dose^2) + age), "dose"),
               "also enters")
  expect_error(ape_dgp_from_fit(f(y ~ expo + log(expo) + age), "expo"),
               "also enters")
  expect_error(ape_dgp_from_fit(f(y ~ treatL * female + age), "treatLTRUE"),
               "nteract")
  expect_error(ape_dgp_from_fit(f(y ~ factor(treat) * female + age),
                                "factor(treat)1"), "nteract")
  expect_error(ape_dgp_from_fit(f(y ~ treat + I(treat * female) + age), "treat"),
               "also enters")
  expect_error(ape_dgp_from_fit(f(y ~ poly(dose, 2) + age), "poly(dose, 2)1"),
               "single column")
  expect_error(ape_dgp_from_fit(f(y ~ treat + age + offset(0.1 * age)), "treat"),
               "offset")
  expect_error(ape_dgp_from_fit(glm(y ~ treat + age, binomial("probit"), pd,
                                    offset = 0.1 * age), "treat"),
               "offset")
  ## clean pilots still pass, including nuisance-side interactions and a
  ## logical focal entering as a single main-effect dummy
  pd$industry <- sample(c("a", "b", "c"), 900, TRUE)
  d1 <- ape_dgp_from_fit(f(y ~ treat + age + factor(industry) + age:female),
                         "treat")
  expect_identical(d1$focal$type, "binary")
  d2 <- ape_dgp_from_fit(f(y ~ treatL + age), "treatLTRUE")
  expect_identical(d2$focal$type, "binary")
  expect_equal(true_ape(set_ape(d2, 0.05)), 0.05, tolerance = 1e-7)
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
