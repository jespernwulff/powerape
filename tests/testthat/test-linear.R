# The linear probability analysis on the standard routes:
# set_analysis(), the engine's OLS path, where the linear coefficient
# centres, and the wiring into searches, sweeps and reports.

lin_base <- function(target = 0.10, ...) {
  set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.4),
                  covariates = list(pa_var("age", "normal"),
                                    pa_var("female", "binary", p = 0.5)),
                  baseline = 0.30, signal = 0.20, n_int = 1e5, ...), target)
}

# the confounded case: a treatment related to a covariate that drives the outcome
lin_confounded <- function(target = 0.10, correlation = 0.3, signal = 0.3) {
  set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.3),
                  covariates = list(pa_var("size", "normal")),
                  correlation = correlation, baseline = 0.15, signal = signal,
                  n_int = 1e5), target)
}

test_that("set_analysis returns a copy and refuses what it cannot analyse", {
  d <- lin_base()
  dl <- set_analysis(d, "linear")
  expect_null(d$analysis)                               # the original is untouched
  expect_identical(dl$analysis, "linear")
  expect_null(set_analysis(dl, "index")$analysis)       # and back
  expect_identical(set_analysis(dl, "index"), d)
  expect_error(set_analysis(d, "ols"), "should be one of")
  expect_error(set_analysis(list(), "linear"))
  div <- set_ape(ape_dgp_iv(focal = pa_var("training", "binary", p = 0.4),
                            instruments = pa_var("subsidy", "binary", p = 0.5),
                            endogeneity = 0.3, iv_strength = 0.3,
                            baseline = 0.30, n_int = 4e4, seed_int = 2L), 0.10)
  expect_error(set_analysis(div, "linear"), "two-stage least squares")
  div$analysis <- "linear"                              # set by hand: the engine refuses too
  expect_error(ape_power(div, n = 500, claim = "detect", nsim = 20), "No linear analysis")
  ## difference-in-differences designs: the same as the builder's argument
  a <- set_att(ape_dgp_did(baseline = 0.2, gap = 0.2, trend = 0.1, trends = "linear",
                           analysis = "linear"), 0.10)
  b <- set_analysis(set_att(ape_dgp_did(baseline = 0.2, gap = 0.2, trend = 0.1,
                                        trends = "linear"), 0.10), "linear")
  expect_identical(a, b)
})

test_that("the engine's linear estimate and SE equal lm + sandwich HC1 on replayed data", {
  skip_if_not_installed("sandwich")
  replay <- function(d, n, seed, col, nrep = 3) {
    bt <- powerape:::beta_true(d)
    powerape:::with_seed(seed, {
      out <- matrix(NA_real_, nrep, 2)
      for (r in seq_len(nrep)) {
        xx <- powerape:::draw_x(d, n)
        y <- rbinom(n, 1L, d$G(drop(xx$X %*% bt)))
        fit <- lm(y ~ xx$X - 1)
        out[r, ] <- c(coef(fit)[[col]],
                      sqrt(sandwich::vcovHC(fit, type = "HC1")[col, col]))
      }
      out
    })
  }
  designs <- list(
    binary = list(d = lin_base(), col = 2L),
    continuous = list(d = set_ape(ape_dgp(model = "logit",
                                          focal = pa_var("rd", "normal", mean = 1, sd = 2),
                                          covariates = list(pa_var("size", "normal")),
                                          correlation = 0.3, baseline = 0.25,
                                          signal = 0.2, n_int = 1e5), 0.04), col = 2L),
    aie = list(d = set_aie(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                                   moderator = pa_var("female", "binary", p = 0.55),
                                   covariates = list(pa_var("age", "normal")),
                                   baseline = 0.30, signal = 0.2, n_int = 1e5),
                           target = 0.06, main_focal = 0.10, main_moderator = 0.05),
               col = 4L),                                # the interaction, not the focal main effect
    aie_cont = list(d = set_aie(ape_dgp(focal = pa_var("rd", "normal"),
                                        moderator = pa_var("public", "binary", p = 0.4),
                                        baseline = 0.30, n_int = 1e5),
                                target = 0.04, main_focal = 0.05, main_moderator = 0.05),
                    col = 4L))
  for (k in names(designs)) {
    d <- designs[[k]]$d
    p <- ape_power(set_analysis(d, "linear"), n = 700, claim = "detect", nsim = 20, seed = 17)
    rp <- replay(d, 700, 17, designs[[k]]$col)
    expect_equal(p$draws$est[1:3], rp[, 1], tolerance = 1e-10, info = k)
    expect_equal(p$draws$se[1:3], rp[, 2], tolerance = 1e-8, info = k)
  }
})

test_that("binary focal without covariates: linear and probit estimates coincide study by study", {
  d <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.4), baseline = 0.30), 0.10)
  pi_ <- ape_power(d, n = 800, claim = "detect", nsim = 200, seed = 5)
  pl <- ape_power(set_analysis(d, "linear"), n = 800, claim = "detect", nsim = 200, seed = 5)
  ## both are the difference of the two group means
  expect_equal(pl$draws$est, pi_$draws$est, tolerance = 1e-6)
  ## HC1 is the model-based SE times sqrt(n / (n - 2)) in this saturated model
  expect_equal(pl$draws$se, pi_$draws$se * sqrt(800 / 798), tolerance = 1e-3)
  expect_equal(pl$power, pi_$power, tolerance = 0.02)
  ## the linear run reports where its estimates centre; the index run does not
  expect_false(is.null(pl$est_mean))
  expect_null(pi_$est_mean)
  se_cf <- sqrt(0.4 * 0.6 / (800 * 0.4) + 0.3 * 0.7 / (800 * 0.6))
  expect_equal(mean(pl$draws$se), se_cf, tolerance = 0.01)
})

test_that("where the linear coefficient centres: equal to the truth in three cases", {
  cen <- function(d) powerape:::linear_centre(set_analysis(d, "linear"))
  ## (a) binary focal independent of the covariates
  expect_equal(cen(lin_base()), 0.10, tolerance = 2e-3)
  ## (b) focal and covariates jointly normal: the projection coefficient is
  ## the average derivative (Stein's lemma), under probit and under logit
  for (m in c("probit", "logit")) {
    d <- set_ape(ape_dgp(model = m, focal = pa_var("rd", "normal", mean = 1, sd = 2),
                         covariates = list(pa_var("size", "normal"), pa_var("s2", "normal")),
                         correlation = 0.4, baseline = 0.25, signal = 0.30, n_int = 1e5),
                 0.04)
    expect_equal(cen(d), 0.04, tolerance = 2e-3)
  }
  ## (c) binary-by-binary AIE, focal and moderator independent of the covariate
  d <- set_aie(ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
                       moderator = pa_var("female", "binary", p = 0.55),
                       covariates = list(pa_var("age", "normal")),
                       baseline = 0.30, signal = 0.20, n_int = 1e5),
               target = 0.06, main_focal = 0.10, main_moderator = 0.05)
  expect_equal(cen(d), 0.06, tolerance = 5e-3)
  ## deterministic, and the caller's random-number state is left alone
  set.seed(1); u1 <- runif(1)
  set.seed(1); c1 <- cen(lin_base()); u2 <- runif(1)
  expect_identical(c1, cen(lin_base()))
  expect_identical(u1, u2)
})

test_that("where the linear coefficient centres: away from the truth under confounding", {
  ## a treatment correlated with a covariate that drives the outcome: the
  ## linear control leaves nonlinear confounding in the coefficient
  d1 <- set_analysis(lin_confounded(correlation = 0.3, signal = 0.3), "linear")
  d2 <- set_analysis(lin_confounded(correlation = 0.6, signal = 0.6), "linear")
  c1 <- powerape:::linear_centre(d1)
  c2 <- powerape:::linear_centre(d2)
  expect_gt(c1 / 0.10, 1.08)            # about 12% above the APE
  expect_lt(c1 / 0.10, 1.18)
  expect_gt(c2 / 0.10, 1.45)            # about 66% above
  ## the simulated estimates sit at the centre, not at the truth
  p <- ape_power(d1, n = 2000, claim = "detect", nsim = 400, seed = 9)
  expect_lt(abs(p$est_mean - c1), 4 * p$est_mcse)
  expect_gt(abs(p$est_mean - 0.10), 4 * p$est_mcse)
  ## while the probit analysis of the same design is centred on the APE
  pi_ <- ape_power(set_analysis(d1, "index"), n = 2000, claim = "detect", nsim = 400,
                   seed = 9, keep_draws = TRUE)
  expect_lt(abs(mean(pi_$draws$est) - 0.10), 4 * sd(pi_$draws$est) / sqrt(400))
  ## the design print shows both numbers
  out <- capture.output(print(d1))
  expect_true(any(grepl("linear probability model \\(OLS\\)", out)))
  expect_true(any(grepl(sprintf("coefficient is %+.4f (the true APE is +0.1000)", c1), out,
                        fixed = TRUE)))
  expect_true(any(grepl("linear estimates average", capture.output(print(p)))))
})

test_that("panel route: pooled OLS with unit means is the unit-fixed-effects estimate", {
  d <- set_ape(ape_dgp_panel(focal = pa_var("treat", "binary", p = 0.5, icc = 0.5),
                             covariates = list(pa_var("size", "normal", icc = 0.6)),
                             n_periods = 4, rho = 0.3, cre_share = 0.5,
                             baseline = 0.30, signal = 0.10, n_int = 4e4), 0.08)
  set.seed(11)
  xx <- powerape:::draw_x(d, 150)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  fo <- powerape:::fit_ols(xx$X, y)
  fe <- lm(y ~ factor(xx$id) + xx$X[, 2] + xx$X[, 3])
  expect_equal(fo$coefficients[[2]], unname(coef(fe)["xx$X[, 2]"]), tolerance = 1e-10)
  ## the engine runs it with unit-clustered SEs
  dl <- set_analysis(d, "linear")
  p <- ape_power(dl, n = 300, claim = "detect", nsim = 100, seed = 12)
  out <- capture.output(print(p))
  expect_true(any(grepl("linear probability \\(OLS\\)", out)))
  expect_true(any(grepl("unit-clustered SEs", out)))
  expect_equal(mean(p$draws$se), sd(p$draws$est), tolerance = 0.25)
  expect_lt(abs(p$est_mean - powerape:::linear_centre(dl)), 4 * p$est_mcse + 0.003)
})

test_that("searches, sweeps and reports under the linear analysis", {
  d <- lin_base()
  dl <- set_analysis(d, "linear")
  ## `se` is fixed and asking for robust SEs is not worth a warning
  expect_no_warning(a <- ape_power(dl, n = 600, claim = "detect", nsim = 60, seed = 2,
                                   se = "robust"))
  expect_identical(a$draws, ape_power(dl, n = 600, claim = "detect", nsim = 60, seed = 2)$draws)
  ## the required n agrees with the probit analysis here (same centre)
  nl <- ape_n(dl, power = 0.80, claim = "detect", nsim = 500, seed = 3, confirm = FALSE)
  ni <- ape_n(d, power = 0.80, claim = "detect", nsim = 500, seed = 3, confirm = FALSE)
  expect_lt(abs(nl$n / ni$n - 1), 0.10)
  expect_true(any(grepl("linear probability \\(OLS\\), robust \\(HC1\\) SEs",
                        capture.output(print(nl)))))
  md <- ape_mde(dl, n = 800, power = 0.80, claim = "detect", nsim = 300, seed = 4,
                confirm = FALSE)
  expect_true(any(grepl("analysis: linear probability \\(OLS\\), robust \\(HC1\\) SEs",
                        capture.output(print(md)))))
  expect_identical(md$dgp$analysis, "linear")
  cv <- ape_curve(dl, n = c(300, 900), claim = "detect", nsim = 80, seed = 5)
  expect_identical(cv$model, "linear probability (OLS)")
  expect_true("est_mean" %in% names(cv$results))
  ## ape_n refuses when the linear centre is on the other side of the boundary:
  ## a normal focal with a correlated binary covariate centres 2.5% low
  de <- set_analysis(set_ape(ape_dgp(focal = pa_var("rd", "normal"),
                                     covariates = list(pa_var("public", "binary", p = 0.3)),
                                     correlation = 0.5, baseline = 0.2, signal = 0.4,
                                     n_int = 1e5), 0.05), "linear")
  cen <- powerape:::linear_centre(de)
  expect_lt(cen, 0.0495)
  expect_error(ape_n(de, power = 0.8, claim = "minimum", sesoi = 0.0495, nsim = 100),
               "different sides")
  expect_no_error(powerape:::linear_need_side(de, "detect", NULL, "ape_n()"))
  ## a robustness sweep keeps the analysis in every scenario: in the
  ## confounded design the linear sweep differs from the probit sweep
  dc <- lin_confounded()
  rl <- suppressWarnings(ape_robust(set_analysis(dc, "linear"), n = 600, claim = "detect",
                                    nsim = 150, seed = 6, vary = list(baseline = c(0.12, 0.18)),
                                    nmax = FALSE))
  ri <- suppressWarnings(ape_robust(dc, n = 600, claim = "detect", nsim = 150, seed = 6,
                                    vary = list(baseline = c(0.12, 0.18)), nmax = FALSE))
  expect_identical(rl$analysis, "linear")
  expect_null(ri$analysis)
  ## scenario 2 is the base design: its cell equals a direct run of the
  ## same analysis at the scenario's seed, and the two analyses differ
  expect_equal(rl$scenarios$power[2],
               ape_power(set_analysis(dc, "linear"), n = 600, claim = "detect",
                         nsim = 150, seed = 8)$power)
  expect_equal(ri$scenarios$power[2],
               ape_power(dc, n = 600, claim = "detect", nsim = 150, seed = 8)$power)
  expect_false(isTRUE(all.equal(rl$scenarios$power, ri$scenarios$power)))
  expect_true(any(grepl("linear probability analysis", capture.output(print(rl)))))
  ## reports: the statement would describe an analysis that was not run
  pw <- ape_power(dl, n = 600, claim = "detect", nsim = 60, seed = 2)
  expect_error(power_statement(pw), "linear probability")
  expect_error(power_statement(nl), "linear probability")
  expect_error(power_statement(md), "linear probability")
  expect_s3_class(power_at(pw, c(0.01, 0.05)), "data.frame")
  al <- ape_alpha(pw, cost = 4)
  expect_error(power_statement(al), "linear probability")
})

test_that("the linear analysis loses no study to separation", {
  d <- set_ape(ape_dgp(focal = pa_var("t", "binary", p = 0.1), baseline = 0.05), -0.04)
  pi_ <- suppressWarnings(ape_power(d, n = 300, claim = "detect", nsim = 200, seed = 6))
  pl <- suppressWarnings(ape_power(set_analysis(d, "linear"), n = 300, claim = "detect",
                                   nsim = 200, seed = 6))
  expect_gt(pi_$separated, 0.05)
  expect_equal(pl$separated, 0)
  expect_lt(pl$outcomes[["failed"]], pi_$outcomes[["failed"]])
})
