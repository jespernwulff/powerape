# The difference-in-differences route in the robustness sweep, the
# justified-alpha functions and the methods paragraph.

di_base <- function(t = 0.07, gap = 0.10, ...) {
  set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = gap, trend = 0.08,
                      rho = 0.4, ...), t)
}
di_sub <- function() {
  set_att(ape_dgp_did(model = "logit", baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                      trend = c(0.10, 0.05), rho = 0.5,
                      subgroup = pa_var("female", "binary", p = 0.4)), c(0.10, 0.03))
}

test_that("ape_robust sweeps a diff-in-diff design and re-pins the effect", {
  d <- di_base()
  r <- suppressWarnings(ape_robust(d, n = 1200, claim = "detect", nsim = 150, seed = 1,
                                   vary = list(rho = c(0.2, 0.6)), nmax = FALSE))
  expect_s3_class(r, "powerape_robust")
  expect_identical(r$estimand, "att")
  expect_equal(r$scenarios$rho, c(0.2, 0.4, 0.6))
  expect_equal(r$scenarios$implied_effect, rep(0.07, 3))
  expect_false("valid" %in% names(r$scenarios))
  ## the middle scenario is the base design: a direct run at its seed
  p2 <- suppressWarnings(ape_power(d, 1200, "detect", nsim = 150, seed = 3))
  expect_equal(r$scenarios$power[2], p2$power)
  expect_output(print(r), "robustness sweep -- detect claim, ATT", fixed = TRUE)
  ## every input of the route can be varied; rho only in a panel
  r2 <- suppressWarnings(ape_robust(d, n = 1200, claim = "detect", nsim = 60, seed = 2,
                                    vary = list(baseline = c(0.15, 0.35), gap = c(0, 0.2),
                                                trend = c(0.04, 0.12), p_treated = c(0.3, 0.5)),
                                    grid_points = 2, nmax = FALSE))
  expect_equal(nrow(r2$scenarios), 16L)
  expect_true(all(is.finite(r2$scenarios$power)))
  dx <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.08,
                            design = "cross-sections"), 0.07)
  expect_error(ape_robust(dx, n = 1200, claim = "detect", vary = list(rho = c(0.2, 0.6))),
               "Cannot vary `rho`")
  ## the linear analysis is carried into every scenario
  dl <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.08, rho = 0.4,
                            trends = "linear", analysis = "linear"), 0.07)
  rl <- suppressWarnings(ape_robust(dl, n = 1000, claim = "detect", nsim = 80, seed = 3,
                                    vary = list(gap = c(0, 0.2)), grid_points = 2, nmax = FALSE))
  expect_identical(rl$analysis, "linear")
  expect_output(print(rl), "linear probability analysis (OLS), unit-clustered SEs", fixed = TRUE)
})

test_that("a sweep over the trends scale flags the scenarios the analysis is biased in", {
  d <- di_base(gap = 0.20)
  r <- suppressWarnings(ape_robust(d, n = 1500, claim = "detect", nsim = 150, seed = 2,
                                   vary = list(trends = c("index", "linear")), nmax_power = 0.8))
  sc <- r$scenarios
  expect_identical(sc$trends, c("index", "linear"))
  expect_identical(sc$valid, c(TRUE, FALSE))
  expect_equal(sc$implied_effect, c(0.07, 0.07))
  ## in the linear world the logit analysis centres on the index reading of
  ## the same rates, which is lower here
  dlin <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.20, trend = 0.08,
                              rho = 0.4, trends = "linear"), 0.07)
  expect_lt(abs(sc$est_mean[2] - powerape:::did_readings(dlin)$index), 0.012)
  expect_lt(sc$est_mean[2], sc$est_mean[1])
  ## the biased scenario is shown but never ranked
  expect_true(r$worst$valid)
  expect_identical(r$nmax_scenario$trends, "index")
  out <- capture.output(print(r))
  expect_true(any(grepl("valid = FALSE", out, fixed = TRUE)))
  expect_true(any(grepl("trends = index", out, fixed = TRUE)))
  expect_false(any(grepl("infeasible", out, fixed = TRUE)))
  expect_error(ape_robust(d, n = 1500, claim = "detect", vary = list(trends = c("index", "probit"))),
               "must be c")
})

test_that("MDE mode and subgroup designs run; unsupported requests stop", {
  d <- di_base()
  m <- suppressWarnings(ape_robust(d, n = 1500, claim = "detect", nsim = 120, seed = 3, mode = "mde",
                                   vary = list(baseline = c(0.15, 0.35)), grid_points = 2))
  expect_true(all(m$scenarios$mde > 0.04 & m$scenarios$mde < 0.2))
  expect_output(print(m), "minimum detectable ATT", fixed = TRUE)
  ds <- di_sub()
  rs <- suppressWarnings(ape_robust(ds, n = 2500, claim = "detect", nsim = 80, seed = 4,
                                    vary = list(rho = c(0.3, 0.7)), grid_points = 2, nmax = FALSE))
  expect_identical(rs$estimand, "datt")
  expect_equal(rs$scenarios$implied_effect, rep(-0.07, 2))
  ms <- suppressWarnings(ape_robust(ds, n = 2500, claim = "detect", nsim = 80, seed = 5, mode = "mde",
                                    vary = list(rho = c(0.3, 0.7)), grid_points = 2))
  expect_true(all(ms$scenarios$mde < 0))
  expect_error(ape_robust(ds, n = 2500, claim = "detect", vary = list(baseline = c(0.2, 0.3))),
               "Cannot vary `baseline`")
  expect_error(ape_robust(d, n = 1500, claim = "detect", pin = "coefficients",
                          vary = list(rho = c(0.2, 0.6))), "not available for")
  dm <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.20, trend = 0.08, rho = 0.4,
                            trends = "index", analysis = "linear"), 0.07)
  expect_error(ape_robust(dm, n = 1500, claim = "detect", vary = list(rho = c(0.2, 0.6))),
               "needs an analysis that is valid")
  d0 <- ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.08, rho = 0.4)
  expect_error(ape_robust(d0, n = 1500, claim = "detect", vary = list(rho = c(0.2, 0.6))),
               "call set_att")
})

test_that("power_at and ape_alpha read a diff-in-diff run; the boundary world is re-pinned", {
  d <- di_base()
  pw <- suppressWarnings(ape_power(d, 1500, "minimum", sesoi = 0.02, nsim = 300, seed = 5))
  pa <- power_at(pw, c(0.01, 0.05))
  expect_equal(pa$power[2], pw$power)
  expect_lt(pa$power[1], pa$power[2])
  ja <- ape_alpha(pw, cost = 4, size = TRUE, nsim_size = 500, seed = 6)
  expect_s3_class(ja, "powerape_alpha")
  expect_equal(ja$size$truth, 0.02)
  expect_lt(abs(ja$size$size - ja$alpha), 4 * ja$size$mcse + 0.005)
  expect_output(print(ja), "true ATT = 0.020", fixed = TRUE)
  ## the claim-boundary world of each claim
  expect_equal(powerape:::boundary_dgp(d, "detect", NULL, 0.07)$target_est, 0)
  expect_equal(powerape:::boundary_dgp(d, "minimum", 0.02, -0.07)$target_est, -0.02)
  ds <- di_sub()
  bs <- powerape:::boundary_dgp(ds, "minimum", 0.02, -0.07)
  expect_equal(bs$target_est, -0.02)
  expect_equal(bs$att_sub, c(0.10, 0.08))          # the reference effect is held
  ## a result that is a rejection rate is not read as power
  dm <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.20, trend = 0.08, rho = 0.4,
                            trends = "index", analysis = "linear"), 0.07)
  pm <- suppressWarnings(ape_power(dm, 800, "detect", nsim = 40, seed = 7))
  expect_error(power_at(pm, 0.01), "rejection rate")
  expect_error(ape_alpha(pm), "rejection rate")
  expect_error(power_statement(pm), "rejection rate")
})

test_that("power_statement describes the design, the world and the analysis", {
  d <- di_base(covariates = list(pa_var("z", "normal")), signal = 0.2, n_int = 4e4)
  pw <- suppressWarnings(ape_power(d, 1500, "minimum", sesoi = 0.02, nsim = 100, seed = 8))
  s <- power_statement(pw)
  expect_s3_class(s, "powerape_statement")
  has <- function(txt, s_ = s) expect_true(grepl(txt, s_, fixed = TRUE), info = txt)
  has("average treatment effect on the treated (ATT) in a difference-in-differences design")
  has("a treated group (50% of the units) and a control group")
  has("observed in one period before and one period after one common treatment date")
  has("as a panel of the same units, with a latent within-unit correlation of 0.40")
  has("0.250 and 0.330 in the control group and 0.350 and")
  has("untreated trends parallel on the logit index")
  has("time-constant covariates z (continuous, mean 0, SD 1), independent of group membership")
  has("The assumed true ATT was 0.070 (7.0 percentage points).")
  has("n = 1500 units observed in 2 periods (3,000 unit-period observations)")
  has("a pooled logit regression of the outcome on group, period indicators and the treatment indicator, and the covariates")
  has("unit-clustered standard errors")
  has("a group-by-before/after cell with no events or no non-events")
  has("the probability of concluding a meaningful effect was")
  expect_false(grepl("focal", s, fixed = TRUE))
  ## required n and minimum detectable effect
  nn <- suppressWarnings(ape_n(d, power = 0.8, claim = "detect", nsim = 150, seed = 9, confirm = FALSE))
  has("The required number of units for 80% power", power_statement(nn))
  md <- suppressWarnings(ape_mde(d, n = 1500, claim = "detect", nsim = 150, seed = 10, confirm = FALSE))
  has("the minimum detectable ATT at 80% power", power_statement(md))
  ## linear analysis, repeated cross-sections, several periods
  dx <- set_att(ape_dgp_did(model = "probit", baseline = 0.20, gap = 0.20, trend = 0.10,
                            design = "cross-sections", n_pre = 2, n_post = 3,
                            trends = "linear", analysis = "linear"), 0.08)
  sx <- power_statement(suppressWarnings(ape_power(dx, 2000, "detect", nsim = 60, seed = 11)))
  has("observed in 2 periods before and 3 periods after", sx)
  has("as repeated cross-sections, with 60% of the observations after the treatment date", sx)
  has("untreated trends parallel in probabilities", sx)
  has("n = 2000 observations", sx)
  has("ordinary least squares of the outcome on group, period indicators and the treatment indicator (a linear probability model)", sx)
  has("heteroskedasticity-robust (HC1) standard errors", sx)
  expect_false(grepl("maximum-likelihood", sx, fixed = TRUE))
  expect_false(grepl("separation", sx, fixed = TRUE))
  ## subgroup difference
  ss <- power_statement(suppressWarnings(ape_power(di_sub(), 2500, "detect", nsim = 60, seed = 12)))
  has("difference between two subgroups' average treatment effects on the treated (by female)", ss)
  has("were 0.100 (female = 0) and 0.030 (female = 1), a true difference of -0.070", ss)
  has("fully interacted with female", ss)
  has("within a subgroup", ss)
  ## the justification paragraph names the estimand and the analysis
  sa <- power_statement(ape_alpha(pw, cost = 4))
  has("the planning value (ATT = 0.070)", sa)
  has("Each simulated study was analysed by a pooled logit regression", sa)
  has("At n = 1500 units the optimum", sa)
  ## standard designs keep their paragraph
  d1 <- set_ape(ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30), 0.10)
  s1 <- power_statement(suppressWarnings(ape_power(d1, 800, "detect", nsim = 60, seed = 13)))
  has("average partial effect (APE) of treat", s1)
  expect_false(grepl("difference-in-differences", s1, fixed = TRUE))
})
