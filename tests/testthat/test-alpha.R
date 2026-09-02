# Justified alpha (1.9.0): sample-size rules for `alpha`, re-thresholding of
# stored draws, and the Maier-Lakens error-cost optimizer with its .05 cap.

make_world <- function(target = 0.10) {
  set_ape(ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                  covariates = list(pa_var("z", "normal")),
                  baseline = 0.30, signal = 0.10), target)
}

# A powerape_power object whose draws are an exact normal quantile sample:
# power(alpha) then tracks the closed form 1 - pnorm(qnorm(1 - alpha) -
# (target - sesoi) / se) to within the step size 1 / nsim.
synthetic_draws <- function(pw, se0, nsim = 4000L) {
  est <- pw$target + se0 * qnorm(ppoints(nsim))
  pw$draws <- list(est = est, se = rep(se0, nsim), ok = rep(TRUE, nsim))
  pw$nsim <- nsim
  pw
}

test_that("claim_levels evaluates a sample-size rule at n and guards it", {
  lv <- powerape:::claim_levels("minimum", function(n) 0.5 / sqrt(n), n = 400)
  expect_equal(lv$alpha, 0.025)
  expect_equal(lv$conf_claim, 0.95)
  expect_equal(powerape:::claim_levels("detect", function(n) 0.01, n = 10)$conf_claim, 0.99)
  expect_error(powerape:::claim_levels("minimum", function(n) 0.01), "needs `n`")
  expect_error(powerape:::claim_levels("minimum", function(n) 0.7, n = 100), "\\(0, 0.5\\)")
  expect_error(powerape:::claim_levels("minimum", function(n) 0.01, conf = 0.95, n = 100),
               "not both")
  ## scalar path unchanged
  expect_equal(powerape:::claim_levels("minimum", 0.05)$conf_claim, 0.90)
})

test_that("a constant rule reproduces the scalar run bit for bit", {
  d <- make_world()
  fix <- ape_power(d, 400, claim = "minimum", sesoi = 0.03, nsim = 120, seed = 1)
  rule <- ape_power(d, 400, claim = "minimum", sesoi = 0.03,
                    alpha = function(n) 0.05, nsim = 120, seed = 1)
  expect_identical(rule$power, fix$power)
  expect_identical(rule$outcomes, fix$outcomes)
  expect_identical(rule$draws$est, fix$draws$est)
  expect_true(is.function(rule$alpha_rule))
  expect_null(fix$alpha_rule)
  expect_equal(rule$alpha, 0.05)
  ## a rule that lands on the conservative convention equals the conf override
  r2 <- ape_power(d, 400, claim = "minimum", sesoi = 0.03,
                  alpha = function(n) 0.5 / sqrt(n), nsim = 120, seed = 1)
  c95 <- ape_power(d, 400, claim = "minimum", sesoi = 0.03, conf = 0.95,
                   nsim = 120, seed = 1)
  expect_identical(r2$power, c95$power)
  expect_equal(r2$alpha, 0.025)
  expect_output(print(r2), "by sample-size rule")
})

test_that("ape_curve carries the rule's alpha per n and ape_n designs jointly", {
  d <- make_world()
  rule <- function(n) 0.05 * sqrt(100 / n)
  cv <- ape_curve(d, n = c(400, 1600), claim = "detect", alpha = rule,
                  nsim = 100, seed = 2)
  expect_equal(cv$results$alpha, rule(c(400, 1600)))
  expect_true(cv$alpha_rule)
  expect_output(print(cv), "sample-size rule")
  cv0 <- ape_curve(d, n = c(400, 1600), claim = "detect", nsim = 100, seed = 2)
  expect_equal(cv0$results$alpha, c(0.05, 0.05))
  ## a stricter constant rule needs more observations than the default
  n_fix <- ape_n(d, power = 0.8, claim = "detect", nsim = 200, seed = 3,
                 confirm = FALSE)
  n_rule <- ape_n(d, power = 0.8, claim = "detect", alpha = function(n) 0.01,
                  nsim = 200, seed = 3, confirm = FALSE)
  expect_gt(n_rule$n, n_fix$n)
  expect_equal(n_rule$alpha, 0.01)
  expect_equal(n_rule$conf, 0.99)
  expect_true(all(n_rule$history$alpha == 0.01))
  expect_true(is.function(n_rule$alpha_rule))
  expect_match(unclass(power_statement(n_rule)), "set as a function of the sample size")
  expect_output(print(n_rule), "by sample-size rule")
})

test_that("ape_mde and ape_robust resolve the rule at their fixed n", {
  d <- make_world()
  rule <- function(n) 0.5 / sqrt(n)     # .025 at n = 400
  m <- ape_mde(d, n = 400, claim = "detect", alpha = rule, nsim = 150,
               seed = 4, confirm = FALSE)
  expect_equal(m$alpha, 0.025)
  expect_equal(m$conf, 0.975)
  expect_true(is.function(m$alpha_rule))
  m0 <- ape_mde(d, n = 400, claim = "detect", nsim = 150, seed = 4,
                confirm = FALSE)
  expect_gt(m$mde, m0$mde)
  expect_match(unclass(power_statement(m)), "function of the sample size")
  rb <- ape_robust(d, n = 400, claim = "detect", alpha = rule,
                   vary = list(baseline = c(0.25, 0.35)), grid_points = 2,
                   nsim = 60, seed = 5, nmax = FALSE)
  expect_equal(rb$alpha, 0.025)
  expect_true(is.function(rb$alpha_rule))
})

test_that("power_at re-thresholds the stored draws exactly and monotonically", {
  d <- make_world()
  pw <- ape_power(d, 400, claim = "minimum", sesoi = 0.03, nsim = 200, seed = 6)
  pa <- power_at(pw, c(0.01, 0.025, 0.05, 0.10))
  expect_equal(pa$power[3], pw$power)
  expect_equal(pa$alpha, c(0.01, 0.025, 0.05, 0.10))
  expect_equal(pa$conf, 1 - 2 * pa$alpha)
  expect_true(all(diff(pa$power) >= 0))
  expect_equal(rowSums(pa[, c("minimum", "detect_only", "inconclusive",
                              "equivalence", "failed")]), rep(1, 4))
  pw0 <- ape_power(d, 400, claim = "minimum", sesoi = 0.03, nsim = 60, seed = 6,
                   keep_draws = FALSE)
  expect_null(pw0$draws)
  expect_error(power_at(pw0, 0.05), "keep_draws")
  expect_error(power_at(pw, 0.6))
})

test_that("ape_alpha reproduces the closed-form optimum on synthetic draws", {
  d <- make_world()
  pw <- ape_power(d, 2000, claim = "minimum", sesoi = 0.05, nsim = 60, seed = 7)
  se0 <- 0.015                       # precise design: optimum well below .05
  syn <- synthetic_draws(pw, se0)
  ja <- ape_alpha(syn, cost = 4, prior = 1)
  w_an <- function(a) (4 * a + (1 - (1 - pnorm(qnorm(1 - a) - 0.05 / se0)))) / 5
  opt <- optimize(w_an, c(1e-4, 0.4))
  expect_false(ja$capped)
  expect_lt(abs(ja$alpha - opt$minimum), 0.005)
  expect_lt(abs(ja$wcer - opt$objective), 5e-4)
  expect_true(ja$range[1] <= ja$alpha && ja$alpha <= ja$range[2])
  expect_equal(nrow(ja$grid), length(ja$grid$alpha))
  ## balance at Cohen's 4:1 with power .80 at alpha .05 returns the convention
  se_c <- 0.05 / (qnorm(0.95) + qnorm(0.80))
  bal <- ape_alpha(synthetic_draws(pw, se_c), cost = 4, prior = 1, error = "balance")
  expect_lt(abs(bal$alpha - 0.05), 0.006)
  expect_lt(abs(bal$beta - 0.20), 0.02)
})

test_that("the cap binds when the unconstrained optimum exceeds .05", {
  d <- make_world()
  pw <- ape_power(d, 1000, claim = "minimum", sesoi = 0.05, nsim = 60, seed = 8)
  syn <- synthetic_draws(pw, 0.03)   # imprecise design: optimum above .05
  ja <- ape_alpha(syn, cost = 1, prior = 1)
  expect_true(ja$capped)
  expect_equal(ja$alpha, 0.05)
  expect_gt(ja$alpha_uncapped, 0.05)
  expect_lt(ja$wcer_uncapped, ja$wcer)
  expect_output(print(ja), "BINDS")
  expect_match(unclass(power_statement(ja)), "lies above the cap")
  expect_warning(ape_alpha(syn, cap = 0.10), "justification")
  expect_error(ape_alpha(pw[names(pw) != "draws"]))
  expect_error(ape_alpha(42))
})

test_that("ape_alpha agrees with JustifyAlpha's optimizer on the same draws", {
  skip_if_not_installed("JustifyAlpha")
  d <- make_world()
  pw <- ape_power(d, 1500, claim = "minimum", sesoi = 0.05, nsim = 600, seed = 9)
  assign(".powerape_pf", function(alpha) power_at(pw, pmin(alpha, 0.4999))$power, envir = globalenv())
  on.exit(rm(".powerape_pf", envir = globalenv()), add = TRUE)
  ja <- JustifyAlpha::optimal_alpha(".powerape_pf(alpha = x)", costT1T2 = 4,
                                    priorH1H0 = 1)
  aa <- ape_alpha(pw, cost = 4, prior = 1)
  ## same objective on the same draws: the minimized weighted error rates
  ## agree to Monte Carlo flatness, the optima to the flat region's width
  expect_lt(abs(ja$errorrate - aa$wcer_uncapped), 2e-3)
  expect_lt(abs(ja$alpha - aa$alpha_uncapped), 0.03)
})

test_that("the realized size at the boundary tracks the chosen alpha", {
  d <- make_world()
  pw <- ape_power(d, 600, claim = "minimum", sesoi = 0.03, nsim = 300, seed = 10)
  ja <- ape_alpha(pw, cost = 4, prior = 1, size = TRUE, nsim_size = 400, seed = 11)
  expect_equal(ja$size$truth, 0.03)
  expect_lt(abs(ja$size$size - ja$alpha), 3 * ja$size$mcse + 0.01)
  expect_output(print(ja), "realized size")
  st <- unclass(power_statement(ja))
  expect_match(st, "Maier")
  expect_match(st, "realized error rate")
  expect_s3_class(power_statement(ja), "powerape_statement")
  ## detection: boundary is zero, in the hypothesized direction
  pd <- ape_power(make_world(-0.10), 600, claim = "detect", nsim = 200, seed = 12)
  jd <- ape_alpha(pd, size = TRUE, nsim_size = 300, seed = 13)
  expect_equal(jd$size$truth, 0)
  expect_lt(jd$size$size, 0.10)
  expect_silent(plot(ja))
})
