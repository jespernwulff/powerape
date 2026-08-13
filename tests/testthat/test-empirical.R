make_pilot <- function(n = 2000, seed = 5) {
  set.seed(seed)
  age <- rnorm(n, 45, 12)
  female <- rbinom(n, 1, 0.55)
  treat <- rbinom(n, 1, plogis(0.4 * scale(age)[, 1]))
  data.frame(treat = treat, age = age, female = female)
}

test_that("empirical route: joint binary focal calibrates and inverts exactly on the rows", {
  pilot <- make_pilot()
  d <- ape_dgp_empirical("probit", data = pilot, focal = "treat",
                         baseline = 0.30, signal = 0.20)
  id <- powerape:::integration_draw(d)
  expect_equal(mean(pnorm(d$beta0 + id$idxz)), 0.30, tolerance = 1e-6)
  expect_equal(d$signal_implied, 0.20, tolerance = 1e-6)
  d <- set_ape(d, 0.05)
  expect_equal(true_ape(d), 0.05, tolerance = 1e-7)

  pw <- ape_power(d, n = 800, claim = "detect", nsim = 300, seed = 8)
  expect_true(pw$power > 0 && pw$power <= 1)
  expect_lt(pw$outcomes[["failed"]], 0.02)
})

test_that("empirical route: continuous focal from a data column", {
  pilot <- make_pilot()[, c("age", "female")]
  d <- ape_dgp_empirical("probit", data = pilot, focal = "age",
                         baseline = 0.30, signal = 0.10)
  expect_identical(d$focal$type, "normal")
  d <- set_ape(d, 0.005)
  expect_equal(true_ape(d), 0.005, tolerance = 1e-8)
  expect_output(print(d), "empirical")

  pw <- ape_power(d, n = 700, claim = "detect", nsim = 200, seed = 9)
  expect_true(pw$power > 0 && pw$power <= 1)
})

test_that("empirical route: independent pa_var focal over pilot covariates", {
  pilot <- make_pilot()[, c("age", "female")]
  d <- ape_dgp_empirical("probit", data = pilot,
                         focal = pa_var("treat", "binary", p = 0.5),
                         baseline = 0.30, signal = 0.20)
  d <- set_ape(d, 0.06)
  expect_equal(true_ape(d), 0.06, tolerance = 1e-7)
  pw <- ape_power(d, n = 600, claim = "detect", nsim = 200, seed = 10)
  expect_true(pw$power > 0 && pw$power <= 1)
})

test_that("empirical route input validation", {
  pilot <- make_pilot(600)
  expect_error(ape_dgp_empirical(data = pilot, focal = "nope", baseline = 0.3),
               "not found")
  p2 <- pilot
  p2$lab <- sample(letters, nrow(p2), TRUE)
  expect_error(ape_dgp_empirical(data = p2, focal = "treat", baseline = 0.3),
               "numeric")
  p3 <- pilot
  p3$konst <- 1
  expect_error(ape_dgp_empirical(data = p3, focal = "treat", baseline = 0.3,
                                 signal = 0.2),
               "Constant covariate")
  expect_warning(ape_dgp_empirical(data = pilot[1:100, ], focal = "treat",
                                   baseline = 0.3),
                 "500 pilot rows")
  p4 <- pilot
  p4$age[5] <- NA
  expect_error(ape_dgp_empirical(data = p4, focal = "treat", baseline = 0.3),
               "complete cases")
})
