# Unbalanced panels (Wooldridge, 2019): length mixtures, monotone
# attrition, observed-period Mundlak means, T-cohort dummies. DESIGN.md
# section 13.

unbal_base <- function(..., target = 0.10) {
  set_ape(ape_dgp_panel(
    focal = pa_var("treat", "binary", p = 0.5, icc = 1),
    covariates = list(pa_var("size", "normal", icc = 0.6)),
    rho = 0.3, cre_share = 0.5, baseline = 0.30, signal = 0.10,
    n_int = 4e4, seed_int = 5L, ...), target)
}

test_that("degenerate mixtures short-circuit to the balanced path exactly", {
  db <- unbal_base(n_periods = 4)
  dr <- unbal_base(n_periods = 4, retention = 1)
  dp <- unbal_base(n_periods = c(4, 6), p_periods = c(1, 0))
  expect_identical(db$beta0, dr$beta0)
  expect_identical(db$beta_focal, dr$beta_focal)
  expect_identical(db$beta0, dp$beta0)
  expect_false(isTRUE(dr$unbalanced))
  expect_false(isTRUE(dp$unbalanced))
  pb <- ape_power(db, 200, claim = "detect", nsim = 60, seed = 11)
  pr <- ape_power(dr, 200, claim = "detect", nsim = 60, seed = 11)
  expect_identical(pb$power, pr$power)
})

test_that("unbalanced calibration is exact and E[T] is the mixture mean", {
  du <- unbal_base(n_periods = 6, retention = 0.8)
  expect_true(du$unbalanced)
  ## truncated-geometric attrition: P(T=t) = .8^(t-1) * .2 below the max,
  ## .8^5 at the max
  expect_equal(du$expected_T, sum(1:6 * c(0.2 * 0.8^(0:4), 0.8^5)),
               tolerance = 1e-12)
  expect_equal(true_ape(du), 0.10, tolerance = 1e-7)
  dm <- unbal_base(n_periods = c(2, 4, 6), p_periods = c(0.25, 0.5, 0.25))
  expect_equal(dm$expected_T, 4)
  expect_equal(true_ape(dm), 0.10, tolerance = 1e-7)
})

test_that("draws respect the length mixture and add cohort dummies", {
  du <- unbal_base(n_periods = c(2, 5), p_periods = c(0.5, 0.5))
  set.seed(42)
  xx <- powerape:::draw_x_panel(du, 400)
  T_i <- as.integer(table(xx$id))
  expect_true(all(T_i %in% c(2L, 5L)))
  expect_gt(sum(T_i == 2L), 100)
  expect_gt(sum(T_i == 5L), 100)
  ## X: intercept, d, size, sizebar, one cohort dummy (levels {2,5})
  expect_identical(ncol(xx$X), 5L)
  expect_identical(length(xx$start), 5L)
  ## the dummy column marks T_i = 5 rows (reference = first level, 2)
  dumcol <- xx$X[, 5L]
  expect_setequal(unique(dumcol), c(0, 1))
  expect_identical(unname(rowsum(dumcol, xx$id)[, 1] > 0),
                   unname(T_i == 5L))
  ## Mundlak mean of `size` is the observed-period mean, unit by unit
  size_col <- xx$X[, 3L]
  sbar_col <- xx$X[, 4L]
  hand <- (rowsum(size_col, xx$id)[, 1] / T_i)[xx$id]
  expect_equal(sbar_col, hand, tolerance = 1e-12)
})

test_that("realized outcome moments match the requested unbalanced world", {
  du <- unbal_base(n_periods = 6, retention = 0.8)
  set.seed(99)
  xx <- powerape:::draw_x_panel(du, 40000)
  arm <- xx$X[, 2L]
  p0 <- mean(xx$pr[arm == 0]); p1 <- mean(xx$pr[arm == 1])
  expect_lt(abs(p0 - 0.30), 0.012)
  expect_lt(abs((p1 - p0) - 0.10), 0.012)
})

test_that("clustered vcov on an unbalanced draw matches sandwich exactly", {
  skip_if_not_installed("sandwich")
  du <- unbal_base(n_periods = c(2, 4), p_periods = c(0.4, 0.6))
  set.seed(13)
  xx <- powerape:::draw_x_panel(du, 500)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  dat <- as.data.frame(xx$X[, -1L])
  names(dat) <- c("d", "size", "sizebar", "T4")
  dat$y <- y; dat$id <- xx$id
  ## same start values on both sides: identical IRLS paths, so the
  ## comparison isolates the vcov convention (the ginteff single-fit lesson)
  fit <- glm(y ~ d + size + sizebar + T4, binomial("probit"), data = dat,
             start = xx$start)
  ft <- powerape:::fit_index_model(xx$X, y, "probit", start = xx$start)
  expect_equal(unname(ft$fit$coefficients), unname(coef(fit)),
               tolerance = 1e-8)
  V <- powerape:::vcov_cluster(ft$fit, xx$X, xx$id, "probit")
  Vs <- sandwich::vcovCL(fit, cluster = dat$id, type = "HC0", cadjust = TRUE)
  ## sandwich's estfun carries glm's last-iteration working weights; the
  ## analytic scores differ by convergence slack (~1e-4 relative; same
  ## band as the balanced-panel marginaleffects comparison)
  expect_equal(unname(V), unname(Vs), tolerance = 1e-3)
})

test_that("retention sweeps rebuild through the spec (ape_robust)", {
  du <- unbal_base(n_periods = 4, retention = 0.9)
  rb <- ape_robust(du, n = 150, claim = "detect",
                   vary = list(retention = c(0.7, 1)), grid_points = 2,
                   nsim = 60, seed = 3, nmax = FALSE)
  expect_identical(nrow(rb$scenarios), 2L)
  expect_true(all(is.finite(rb$scenarios$power)))
})

test_that("unbalanced elicitation guards teach", {
  expect_error(ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                             n_periods = 4, retention = .9, p_periods = 1,
                             rho = .2, baseline = .3),
               "not both")
  expect_error(ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                             n_periods = 1, rho = .2, baseline = .3),
               "at least two periods")
  expect_error(ape_dgp_panel(focal = pa_var("t", "binary", p = .5, icc = 1),
                             n_periods = 4, p_periods = c(.5, .5),
                             rho = .2, baseline = .3),
               "p_periods")
})

test_that("power_statement names the Wooldridge specification", {
  du <- unbal_base(n_periods = 5, retention = 0.85)
  pw <- ape_power(du, 150, claim = "detect", nsim = 60, seed = 21)
  txt <- paste(capture.output(power_statement(pw)), collapse = " ")
  expect_match(txt, "Wooldridge, 2019", fixed = TRUE)
  expect_match(txt, "observed periods", fixed = TRUE)
})
