# Difference-in-differences designs: calibration of
# the two worlds, the effect-on-the-treated estimator, closed-form anchors
# for the saturated two-by-two design, and the wiring into the power
# functions.

did_base <- function(model = "logit", design = "panel", rho = 0.5,
                     trends = "index", target = 0.07, ...) {
  d <- ape_dgp_did(model = model, baseline = 0.20, gap = 0.20, trend = 0.10,
                   design = design, rho = if (design == "panel") rho else 0,
                   trends = trends, ...)
  if (!is.null(target)) d <- set_att(d, target)
  d
}

## ---- closed forms for the saturated design (no covariates) ------------------

did_links <- function(link) {
  if (link == "probit") list(G = pnorm, Gi = qnorm, g = dnorm)
  else list(G = plogis, Gi = qlogis, g = dlogis)
}

# delta-method SE of c11 - G(Gi(c10) + Gi(c01) - Gi(c00)), independent cells
# with shares pt (treated) and pp (post)
did_se_cs <- function(link, r, n, pt = 0.5, pp = 0.5) {
  L <- did_links(link)
  gc <- L$g(L$Gi(r$cf))
  v <- function(p, share) p * (1 - p) / (n * share)
  sqrt(v(r$c11, pt * pp) +
         gc^2 * (v(r$c10, pt * (1 - pp)) / L$g(L$Gi(r$c10))^2 +
                   v(r$c01, (1 - pt) * pp) / L$g(L$Gi(r$c01))^2 +
                   v(r$c00, (1 - pt) * (1 - pp)) / L$g(L$Gi(r$c00))^2))
}

# P(Y_pre = 1, Y_post = 1) under a Gaussian copula with latent correlation rho
did_p_joint <- function(p0, p1, rho) {
  if (rho == 0) return(p0 * p1)
  h <- qnorm(p0)
  k <- qnorm(p1)
  integrate(function(a) dnorm(a) * pnorm((h - sqrt(rho) * a) / sqrt(1 - rho)) *
              pnorm((k - sqrt(rho) * a) / sqrt(1 - rho)),
            -Inf, Inf, rel.tol = 1e-10)$value
}

# the same estimator on a panel of n units: the two rates of a group come
# from the same units, so their covariance enters
did_se_panel <- function(link, r, n, rho, pt = 0.5) {
  L <- did_links(link)
  gc <- L$g(L$Gi(r$cf))
  n1 <- n * pt
  n0 <- n * (1 - pt)
  g11 <- 1
  g10 <- -gc / L$g(L$Gi(r$c10))
  g01 <- -gc / L$g(L$Gi(r$c01))
  g00 <- gc / L$g(L$Gi(r$c00))
  cov1 <- (did_p_joint(r$c10, r$c11, rho) - r$c10 * r$c11) / n1
  cov0 <- (did_p_joint(r$c00, r$c01, rho) - r$c00 * r$c01) / n0
  sqrt(g11^2 * r$c11 * (1 - r$c11) / n1 + g10^2 * r$c10 * (1 - r$c10) / n1 +
         g01^2 * r$c01 * (1 - r$c01) / n0 + g00^2 * r$c00 * (1 - r$c00) / n0 +
         2 * g11 * g10 * cov1 + 2 * g01 * g00 * cov0)
}

## ---- worlds ------------------------------------------------------------------

test_that("the index world is calibrated in closed form without covariates", {
  for (m in c("logit", "probit")) {
    L <- did_links(m)
    d <- did_base(model = m, target = NULL)
    expect_equal(d$did$th, 0)
    expect_equal(L$G(d$did$a), 0.20, tolerance = 1e-12)
    expect_equal(L$G(d$did$a + d$did$b), 0.40, tolerance = 1e-12)
    expect_equal(L$G(d$did$a + d$did$l), 0.30, tolerance = 1e-12)
    ## counterfactual treated rate: parallel on the index
    expect_equal(d$rates$cf, L$G(L$Gi(0.40) + L$Gi(0.30) - L$Gi(0.20)),
                 tolerance = 1e-12)
    d <- set_att(d, 0.07)
    expect_equal(true_att(d), 0.07, tolerance = 1e-9)
    expect_equal(d$rates$c11, d$rates$cf + 0.07, tolerance = 1e-12)
    expect_identical(d$estimand, "att")
  }
  ## the worked numbers: .20 -> .30 and .40 -> .60 read as .0667 (logit),
  ## .0745 (probit) and .1000 (linear)
  expect_equal(0.60 - did_base("logit", target = NULL)$rates$cf, 0.0667, tolerance = 1e-3)
  expect_equal(0.60 - did_base("probit", target = NULL)$rates$cf, 0.0745, tolerance = 1e-3)
})

test_that("the linear world has untreated rates parallel in probabilities", {
  for (m in c("logit", "probit")) {
    L <- did_links(m)
    d <- did_base(model = m, trends = "linear", target = 0.10)
    expect_equal(d$rates$cf, 0.50, tolerance = 1e-12)
    expect_true(abs(d$did$th) > 0.05)             # an untreated group-by-period term
    cc <- d$did
    expect_equal(L$G(cc$a + cc$b + cc$l + cc$th) - L$G(cc$a + cc$b),
                 L$G(cc$a + cc$l) - L$G(cc$a), tolerance = 1e-10)
    expect_equal(true_att(d), 0.10, tolerance = 1e-9)
    ## what the two analyses centre on in each other's world
    rd <- powerape:::did_readings(d)
    expect_equal(rd$linear, 0.10, tolerance = 1e-10)
    expect_equal(rd$index, 0.60 - L$G(L$Gi(0.40) + L$Gi(0.30) - L$Gi(0.20)),
                 tolerance = 1e-10)
  }
})

test_that("covariates: rates and the ATT hold on the group-conditional draws", {
  for (tr in c("index", "linear")) {
    d <- ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                     trends = tr, rho = 0.3,
                     covariates = list(pa_var("size", "normal"),
                                       pa_var("public", "binary", p = 0.3)),
                     correlation = 0.3, signal = 0.20, n_int = 1e5)
    d <- set_att(d, 0.06)
    id <- powerape:::did_integration(d)
    cc <- d$did
    expect_equal(mean(pnorm(cc$a + id$q0)), 0.25, tolerance = 1e-7)
    expect_equal(mean(pnorm(cc$a + cc$b + id$q1)), 0.35, tolerance = 1e-7)
    expect_equal(mean(pnorm(cc$a + cc$l + id$q0)), 0.33, tolerance = 1e-7)
    cf <- mean(pnorm(cc$a + cc$b + cc$l + cc$th + id$q1))
    expect_equal(cf, d$rates$cf, tolerance = 1e-7)
    if (tr == "linear") expect_equal(cf, 0.25 + 0.10 + 0.08, tolerance = 1e-7)
    else expect_equal(cc$th, 0)
    expect_equal(true_att(d), 0.06, tolerance = 1e-7)
    expect_equal(d$signal_implied, 0.20, tolerance = 1e-6)
    ## the covariates differ by group (correlation 0.3), so the two
    ## conditional index distributions are not the same
    expect_true(abs(mean(id$q1) - mean(id$q0)) > 0.05)
  }
})

test_that("guards fire with teaching messages", {
  expect_error(ape_dgp_did(baseline = 0.9, gap = 0.2), "treated group's pre-period rate")
  expect_error(ape_dgp_did(baseline = 0.2, trend = -0.3), "control group's post-period rate")
  ## linear trends would put the untreated treated rate above 1
  expect_error(ape_dgp_did(baseline = 0.7, gap = 0.2, trend = 0.2, trends = "linear"),
               "Linear parallel trends cannot hold")
  expect_no_error(ape_dgp_did(baseline = 0.7, gap = 0.2, trend = 0.2))
  expect_warning(ape_dgp_did(baseline = 0.3, design = "cross-sections", rho = 0.4),
                 "ignored")
  expect_error(ape_dgp_did(baseline = 0.3, signal = 0.2), "at least one covariate")
  d <- did_base(target = NULL)
  expect_error(set_att(d, 0.60), "infeasible")
  expect_error(ape_power(d, n = 200, claim = "detect", nsim = 20), "set_att")
  ## the cross-sectional pins and truths refuse the design, and set_att
  ## refuses the other routes
  expect_error(set_ape(d, 0.05), "set_att")
  expect_error(set_aie(d, 0.05, main_focal = 0.1, main_moderator = 0.1), "set_att")
  expect_error(true_ape(set_att(d, 0.05)), "true_att")
  expect_error(true_aie(set_att(d, 0.05)), "true_att")
  expect_error(set_att(ape_dgp(focal = pa_var("t", "binary", p = .5), baseline = .3), 0.05),
               "ape_dgp_did")
  expect_error(true_att(d), "set_att")
})

## ---- estimator ----------------------------------------------------------------

test_that("saturated identity: the ATT is ybar11 - G(Gi(ybar10) + Gi(ybar01) - Gi(ybar00))", {
  for (m in c("logit", "probit")) {
    for (des in c("panel", "cross-sections")) {
      L <- did_links(m)
      d <- did_base(model = m, design = des)
      set.seed(101)
      xx <- powerape:::draw_x_did(d, 900)
      y <- rbinom(length(xx$pr), 1L, xx$pr)
      ft <- powerape:::fit_index_model(xx$X, y, m, start = xx$start)
      expect_true(ft$ok)
      est <- powerape:::att_est(ft$fit$coefficients, xx$X, m, 4L)
      yb <- function(dd, pp) mean(y[xx$X[, 2] == dd & xx$X[, 3] == pp])
      closed <- yb(1, 1) - L$G(L$Gi(yb(1, 0)) + L$Gi(yb(0, 1)) - L$Gi(yb(0, 0)))
      expect_equal(est$ape, closed, tolerance = 1e-7)
      ## it is NOT the cross difference of the four rates
      cross <- yb(1, 1) - yb(1, 0) - yb(0, 1) + yb(0, 0)
      expect_true(abs(est$ape - cross) > 0.01)
    }
  }
})

test_that("the ATT Jacobian matches finite differences (covariates, both links)", {
  for (m in c("logit", "probit")) {
    d <- ape_dgp_did(model = m, baseline = 0.3, gap = 0.1, trend = 0.05,
                     covariates = list(pa_var("x", "normal")), signal = 0.2,
                     n_int = 4e4)
    d <- set_att(d, 0.08)
    set.seed(7)
    xx <- powerape:::draw_x_did(d, 400)
    b <- powerape:::did_beta(d) + c(0.05, -0.03, 0.02, 0.04, -0.06)
    est <- powerape:::att_est(b, xx$X, m, 4L)
    fd <- vapply(seq_along(b), function(j) {
      h <- 1e-6
      bp <- b; bp[j] <- bp[j] + h
      bm <- b; bm[j] <- bm[j] - h
      (powerape:::att_est(bp, xx$X, m, 4L)$ape -
         powerape:::att_est(bm, xx$X, m, 4L)$ape) / (2 * h)
    }, numeric(1))
    expect_equal(unname(est$jac), fd, tolerance = 1e-6)
  }
  ## a study without a treated post-period row has no estimate
  X0 <- cbind(1, c(0, 0, 1), c(0, 1, 0), 0)
  expect_null(powerape:::att_est(c(0, 0, 0, 0), X0, "logit", 4L))
})

test_that("ATT and unit-clustered SE match marginaleffects on a fixed panel", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  for (m in c("logit", "probit")) {
    d <- ape_dgp_did(model = m, baseline = 0.25, gap = 0.15, trend = 0.08,
                     rho = 0.5, covariates = list(pa_var("x", "normal")),
                     signal = 0.2, n_int = 4e4)
    d <- set_att(d, 0.08)
    set.seed(202)
    xx <- powerape:::draw_x_did(d, 800)
    dat <- data.frame(y = rbinom(length(xx$pr), 1L, xx$pr), D = xx$X[, 2],
                      post = xx$X[, 3], W = xx$X[, 4], x = xx$X[, 5], id = xx$id)
    fit <- glm(y ~ D + post + W + x, family = binomial(m), data = dat)
    X <- model.matrix(fit)
    est <- powerape:::att_est(coef(fit), X, m, 4L)
    V <- powerape:::vcov_cluster(fit, X, dat$id, m)
    se <- sqrt(drop(t(est$jac) %*% V %*% est$jac))

    Vs <- sandwich::vcovCL(fit, cluster = dat$id, type = "HC0", cadjust = TRUE)
    expect_equal(unname(V), unname(Vs), tolerance = 1e-5)
    mfx <- marginaleffects::avg_comparisons(
      fit, variables = list(W = c(0, 1)), newdata = subset(dat, W == 1),
      vcov = Vs)
    expect_equal(est$ape, mfx$estimate, tolerance = 1e-6)
    expect_equal(se, mfx$std.error, tolerance = 1e-4)
  }
})

## ---- engine against the closed forms ------------------------------------------

test_that("cross-sections: estimates centre on the ATT and SEs match the closed form", {
  for (m in c("logit", "probit")) {
    d <- did_base(model = m, design = "cross-sections")
    p <- ape_power(d, n = 2000, claim = "detect", nsim = 600, seed = 11)
    ok <- p$draws$ok
    expect_true(all(ok))
    se_cf <- did_se_cs(m, d$rates, 2000)
    expect_equal(mean(p$draws$se), se_cf, tolerance = 0.01)
    expect_equal(sd(p$draws$est), se_cf, tolerance = 0.08)
    ## right estimand: the target, not the cross difference (.10 / .096)
    mc <- sd(p$draws$est) / sqrt(600)
    expect_lt(abs(p$est_mean - 0.07), 4 * mc)
    expect_gt(abs(p$est_mean - powerape:::did_readings(d)$linear), 8 * mc)
    ## power agrees with the normal approximation at the closed-form SE
    expect_lt(abs(p$power - pnorm(0.07 / se_cf - qnorm(0.975))), 3 * p$mcse + 0.01)
    expect_identical(p$estimand, "att")
  }
})

test_that("panel: rho is live and the clustered SE matches the closed form", {
  for (m in c("logit", "probit")) {
    d5 <- did_base(model = m, design = "panel", rho = 0.5)
    d0 <- did_base(model = m, design = "panel", rho = 0)
    p5 <- ape_power(d5, n = 1000, claim = "detect", nsim = 600, seed = 12)
    p0 <- ape_power(d0, n = 1000, claim = "detect", nsim = 600, seed = 12)
    se5 <- did_se_panel(m, d5$rates, 1000, 0.5)
    se0 <- did_se_panel(m, d0$rates, 1000, 0)
    ## without dependence the panel of 1,000 units equals 2,000 independent rows
    expect_equal(se0, did_se_cs(m, d0$rates, 2000), tolerance = 1e-10)
    expect_equal(mean(p5$draws$se), se5, tolerance = 0.015)
    expect_equal(mean(p0$draws$se), se0, tolerance = 0.015)
    expect_equal(sd(p5$draws$est), se5, tolerance = 0.08)
    ## within-unit dependence buys precision: the closed forms differ by 17%
    expect_lt(se5 / se0, 0.87)
    expect_lt(mean(p5$draws$se) / mean(p0$draws$se), 0.87)
    expect_lt(abs(p5$est_mean - 0.07), 4 * sd(p5$draws$est) / sqrt(600))
    expect_gt(p5$power, p0$power)
  }
})

test_that("the copula draw keeps the marginal rates exact under both links", {
  for (m in c("logit", "probit")) {
    d <- did_base(model = m, design = "panel", rho = 0.6)
    set.seed(31)
    xx <- powerape:::draw_x_did(d, 60000)
    y <- rbinom(length(xx$pr), 1L, xx$pr)
    cell <- function(dd, pp) mean(y[xx$X[, 2] == dd & xx$X[, 3] == pp])
    r <- d$rates
    expect_equal(cell(0, 0), r$c00, tolerance = 0.03)
    expect_equal(cell(0, 1), r$c01, tolerance = 0.03)
    expect_equal(cell(1, 0), r$c10, tolerance = 0.03)
    expect_equal(cell(1, 1), r$c11, tolerance = 0.03)
    ## and the within-unit dependence is the copula's
    yt <- matrix(y, nrow = 2)
    tr <- xx$X[seq(1, nrow(xx$X), by = 2), 2] == 1
    expect_equal(mean(yt[1, tr] * yt[2, tr]),
                 did_p_joint(r$c10, r$c11, 0.6), tolerance = 0.03)
  }
})

test_that("covariates: the panel analysis stays centred on the ATT", {
  d <- ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                   rho = 0.4, covariates = list(pa_var("size", "normal")),
                   correlation = 0.3, signal = 0.25, n_int = 1e5)
  d <- set_att(d, 0.08)
  p <- ape_power(d, n = 1200, claim = "detect", nsim = 500, seed = 14)
  expect_true(mean(p$draws$ok) > 0.99)
  expect_lt(abs(p$est_mean - 0.08), 4 * sd(p$draws$est) / sqrt(500))
  ## calibrated SEs: their mean tracks the spread of the estimates
  expect_equal(mean(p$draws$se), sd(p$draws$est), tolerance = 0.10)
})

## ---- mismatched world ---------------------------------------------------------

test_that("linear world: the logit analysis centres on the index reading, not the ATT", {
  d <- did_base(model = "logit", design = "cross-sections", trends = "linear",
                target = 0.10)
  p <- ape_power(d, n = 2000, claim = "detect", nsim = 600, seed = 13)
  mc <- sd(p$draws$est) / sqrt(600)
  expect_lt(abs(p$est_mean - 0.0667), 4 * mc + 0.002)   # plug-in bias allowance
  expect_gt(abs(p$est_mean - 0.10), 8 * mc)
  out <- capture.output(print(p))
  expect_true(any(grepl("rejection rate", out)))
  expect_true(any(grepl("not power", out)))
  expect_false(any(grepl("^  power =", out)))
  ## searches need an analysis that is valid in the simulated world
  expect_error(ape_n(d, power = 0.8, claim = "detect"), "valid in the simulated world")
  expect_error(ape_mde(d, n = 2000), "valid in the simulated world")
  ## the curve labels its column
  cv <- ape_curve(d, n = c(500, 1000), claim = "detect", nsim = 100, seed = 1)
  expect_match(cv$mismatch, "linear-parallel trends and the logit analysis assumes index-parallel")
  expect_true("est_mean" %in% names(cv$results))
  expect_true(any(grepl("rejection", capture.output(print(cv)))))
})

## ---- wiring -------------------------------------------------------------------

test_that("the power functions accept the design", {
  d <- did_base(model = "logit", design = "panel", rho = 0.5, target = 0.12)
  r <- d$rates
  ## ape_curve
  cv <- ape_curve(d, n = c(300, 900), claim = "detect", nsim = 150, seed = 3)
  expect_true(cv$results$power[2] > cv$results$power[1])
  expect_null(cv$mismatch)
  expect_true(any(grepl("ATT", capture.output(print(cv)))))
  expect_silent(plot(cv))
  ## ape_n against the normal-approximation requirement at the closed-form SE
  c1 <- did_se_panel("logit", r, 1, 0.5)                 # SE scale per sqrt(unit)
  n_an <- (c1 * (qnorm(0.975) + qnorm(0.80)) / 0.12)^2
  nn <- ape_n(d, power = 0.80, claim = "detect", nsim = 600, seed = 5,
              confirm = FALSE)
  expect_lt(abs(nn$n / n_an - 1), 0.12)
  expect_identical(nn$estimand, "att")
  expect_true(any(grepl("units \\(x 2 periods\\)", capture.output(print(nn)))))
  ## ape_mde inverts it
  md <- ape_mde(d, n = nn$n, power = 0.80, claim = "detect", nsim = 400,
                seed = 6, confirm = FALSE)
  expect_identical(md$estimand, "att")
  expect_lt(abs(md$mde / 0.12 - 1), 0.15)
  expect_true(any(grepl("minimum detectable ATT", capture.output(print(md)))))
  ## minimum-effect and equivalence claims run on the ATT
  pm <- ape_power(d, n = 1500, claim = "minimum", sesoi = 0.04, nsim = 150, seed = 7)
  expect_true(pm$power > 0.3 && pm$power < 1)
  d0 <- did_base(model = "logit", design = "panel", rho = 0.5, target = 0)
  pe <- ape_power(d0, n = 3000, claim = "equivalence", sesoi = 0.08, nsim = 150, seed = 8)
  expect_true(pe$power > 0.3)
  out <- capture.output(print(pm))
  expect_true(any(grepl("n = 1500 units x 2 periods \\(3000 obs\\)", out)))
  expect_true(any(grepl("ATT", out)))
  ## `se` is fixed for the panel design and honoured for cross-sections
  expect_warning(ape_power(d, n = 300, claim = "detect", nsim = 30, seed = 1,
                           se = "robust"), "fixed by the route")
  dc <- did_base(model = "logit", design = "cross-sections", target = 0.12)
  expect_no_warning(ape_power(dc, n = 600, claim = "detect", nsim = 30, seed = 1,
                              se = "robust"))
  ## the sweep, the justified-alpha functions and the methods paragraph
  ## accept the design too (test-did-integration.R)
  pw <- ape_power(d, n = 600, claim = "detect", nsim = 50, seed = 2)
  expect_s3_class(power_at(pw, 0.01), "data.frame")
  expect_s3_class(power_statement(pw), "powerape_statement")
  expect_s3_class(power_statement(nn), "powerape_statement")
  expect_s3_class(power_statement(md), "powerape_statement")
})

test_that("results are reproducible from the seed and leave the RNG state alone", {
  d <- did_base(model = "probit", design = "panel", rho = 0.4)
  set.seed(99)
  u1 <- runif(1)
  set.seed(99)
  a <- ape_power(d, n = 400, claim = "detect", nsim = 60, seed = 21)
  u2 <- runif(1)
  b <- ape_power(d, n = 400, claim = "detect", nsim = 60, seed = 21)
  expect_identical(a$draws, b$draws)
  expect_identical(u1, u2)
})

test_that("sparse identifying cells are flagged on the four group-by-period cells", {
  ## a small treated group with a rare outcome: some studies have an empty
  ## or all-zero treated cell
  d <- set_att(ape_dgp_did(model = "logit", baseline = 0.03, gap = 0, trend = 0,
                           p_treated = 0.1, design = "cross-sections"), 0.02)
  p <- suppressWarnings(ape_power(d, n = 300, claim = "detect", nsim = 200, seed = 4))
  expect_gt(p$separated, 0.05)
  expect_true(p$outcomes[["failed"]] >= p$separated)
  expect_true(!is.na(p$min_cell) && p$min_cell < 10)
})

## ---- shares, robust SEs, direction, margin search -----------------------------

test_that("unequal shares: the closed forms hold when groups and periods are unbalanced", {
  ## cross-sections: 30% treated, 65% of the observations after
  d <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10,
                           design = "cross-sections", p_treated = 0.3, p_post = 0.65),
               0.07)
  p <- ape_power(d, n = 3000, claim = "detect", nsim = 500, seed = 41)
  se_cf <- did_se_cs("logit", d$rates, 3000, pt = 0.3, pp = 0.65)
  expect_equal(mean(p$draws$se), se_cf, tolerance = 0.015)
  ## the shares are live: the balanced closed form is 8% smaller
  expect_gt(se_cf / did_se_cs("logit", d$rates, 3000), 1.05)
  ## the covariate-free model is saturated, so the HC0 sandwich equals the
  ## model-based variance study by study (up to the IRLS convergence
  ## tolerance of the fitted cell rates)
  pr <- ape_power(d, n = 3000, claim = "detect", nsim = 500, seed = 41, se = "robust")
  expect_equal(pr$draws$se, p$draws$se, tolerance = 1e-3)
  expect_identical(pr$draws$est, p$draws$est)
  ## panel: a quarter of the units treated
  dp <- set_att(ape_dgp_did(model = "probit", baseline = 0.20, gap = 0.20, trend = 0.10,
                            rho = 0.5, p_treated = 0.25), 0.07)
  pp <- ape_power(dp, n = 1500, claim = "detect", nsim = 500, seed = 42)
  expect_equal(mean(pp$draws$se), did_se_panel("probit", dp$rates, 1500, 0.5, pt = 0.25),
               tolerance = 0.015)
})

test_that("negative effects and the equivalence margin search run on the ATT", {
  d <- set_att(ape_dgp_did(model = "logit", baseline = 0.30, gap = 0.10, trend = 0.05,
                           rho = 0.5), -0.08)
  expect_equal(true_att(d), -0.08, tolerance = 1e-9)
  expect_equal(d$rates$c11, d$rates$cf - 0.08, tolerance = 1e-12)
  p <- ape_power(d, n = 1500, claim = "minimum", sesoi = 0.02, nsim = 300, seed = 43)
  expect_lt(abs(p$est_mean + 0.08), 4 * sd(p$draws$est) / sqrt(300))
  se_cf <- did_se_panel("logit", d$rates, 1500, 0.5)
  expect_lt(abs(p$power - pnorm((0.08 - 0.02) / se_cf - qnorm(0.95))), 3 * p$mcse + 0.01)
  ## the search follows the pin's direction
  md <- ape_mde(d, n = 1500, power = 0.80, claim = "detect", nsim = 300, seed = 44,
                confirm = FALSE)
  expect_lt(md$mde, 0)
  expect_identical(md$direction, "negative")
  ## equivalence at a true effect of zero: the margin is se * (z_.95 + z_.90)
  d0 <- set_att(ape_dgp_did(model = "logit", baseline = 0.30, gap = 0.10, trend = 0.05,
                            rho = 0.5), 0)
  mg <- ape_mde(d0, n = 2000, power = 0.80, claim = "equivalence", nsim = 300,
                seed = 45, confirm = FALSE)
  expect_true(mg$is_margin)
  expect_identical(mg$estimand, "att")
  se0 <- did_se_panel("logit", d0$rates, 2000, 0.5)
  expect_lt(abs(mg$mde / (se0 * (qnorm(0.95) + qnorm(0.90))) - 1), 0.12)
})

## ---- the linear probability analysis -----------------

# closed-form SEs of the cross difference of the four cell means
did_se_cs_lin <- function(r, n, pt = 0.5, pp = 0.5) {
  sqrt(r$c11 * (1 - r$c11) / (n * pt * pp) +
         r$c10 * (1 - r$c10) / (n * pt * (1 - pp)) +
         r$c01 * (1 - r$c01) / (n * (1 - pt) * pp) +
         r$c00 * (1 - r$c00) / (n * (1 - pt) * (1 - pp)))
}
# panel: each group's change is taken within units, so the covariance of a
# unit's two outcomes enters negatively
did_se_panel_lin <- function(r, n, rho, pt = 0.5) {
  n1 <- n * pt
  n0 <- n * (1 - pt)
  cov1 <- did_p_joint(r$c10, r$c11, rho) - r$c10 * r$c11
  cov0 <- did_p_joint(r$c00, r$c01, rho) - r$c00 * r$c01
  sqrt((r$c11 * (1 - r$c11) + r$c10 * (1 - r$c10) - 2 * cov1) / n1 +
         (r$c01 * (1 - r$c01) + r$c00 * (1 - r$c00) - 2 * cov0) / n0)
}

did_lin <- function(model = "logit", design = "panel", rho = 0.5,
                    trends = "linear", target = 0.10, ...) {
  set_att(ape_dgp_did(model = model, baseline = 0.20, gap = 0.20, trend = 0.10,
                      design = design, rho = if (design == "panel") rho else 0,
                      trends = trends, analysis = "linear", ...), target)
}

test_that("OLS fit and robust variances match lm and sandwich on fixed data", {
  skip_if_not_installed("sandwich")
  d <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.15, trend = 0.08,
                           rho = 0.5, covariates = list(pa_var("x", "normal")),
                           signal = 0.2, analysis = "linear", n_int = 4e4), 0.08)
  set.seed(303)
  xx <- powerape:::draw_x_did(d, 500)
  dat <- data.frame(y = rbinom(length(xx$pr), 1L, xx$pr), D = xx$X[, 2],
                    post = xx$X[, 3], W = xx$X[, 4], x = xx$X[, 5], id = xx$id)
  fit <- lm(y ~ D + post + W + x, data = dat)
  X <- model.matrix(fit)
  fo <- powerape:::fit_ols(X, dat$y)
  expect_true(fo$ok)
  expect_equal(unname(fo$coefficients), unname(coef(fit)), tolerance = 1e-10)
  expect_equal(unname(powerape:::vcov_ols_hc1(fo, X)),
               unname(sandwich::vcovHC(fit, type = "HC1")), tolerance = 1e-8)
  ## the default of vcovCL carries Stata's regress, vce(cluster) factor
  expect_equal(unname(powerape:::vcov_ols_cluster(fo, X, dat$id)),
               unname(sandwich::vcovCL(fit, cluster = dat$id)), tolerance = 1e-8)
  ## a rank-deficient design (no treated post-period row) has no fit
  expect_false(powerape:::fit_ols(cbind(1, c(0, 0, 1, 1), c(0, 1, 0, 0), 0), c(0, 1, 1, 0))$ok)
})

test_that("linear analysis: the estimate is the cross difference of the cell means", {
  for (des in c("panel", "cross-sections")) {
    d <- did_lin(design = des)
    set.seed(404)
    xx <- powerape:::draw_x_did(d, 300)
    y <- rbinom(length(xx$pr), 1L, xx$pr)
    fo <- powerape:::fit_ols(xx$X, y)
    yb <- function(dd, pp) mean(y[xx$X[, 2] == dd & xx$X[, 3] == pp])
    expect_equal(fo$coefficients[[4]], yb(1, 1) - yb(1, 0) - yb(0, 1) + yb(0, 0),
                 tolerance = 1e-12)
    if (des == "panel") {
      ## ... and, in a balanced two-period panel, the unit-fixed-effects estimate
      fe <- lm(y ~ factor(xx$id) + xx$X[, 3] + xx$X[, 4])
      expect_equal(fo$coefficients[[4]], unname(coef(fe)[length(coef(fe))]),
                   tolerance = 1e-10)
    }
  }
})

test_that("linear analysis: a time-constant covariate is neutral in a balanced panel", {
  d <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.15, trend = 0.08,
                           rho = 0.5, covariates = list(pa_var("x", "normal")),
                           correlation = 0.4, signal = 0.3, trends = "linear",
                           analysis = "linear", n_int = 4e4), 0.08)
  set.seed(505)
  xx <- powerape:::draw_x_did(d, 400)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  X5 <- xx$X
  X4 <- xx$X[, 1:4]
  f5 <- powerape:::fit_ols(X5, y)
  f4 <- powerape:::fit_ols(X4, y)
  expect_equal(f5$coefficients[[4]], f4$coefficients[[4]], tolerance = 1e-12)
  ## same clustered variance once the (N - 1) / (N - K) factor is removed
  N <- nrow(X5)
  v5 <- powerape:::vcov_ols_cluster(f5, X5, xx$id)[4, 4] * (N - 5) / (N - 1)
  v4 <- powerape:::vcov_ols_cluster(f4, X4, xx$id)[4, 4] * (N - 4) / (N - 1)
  expect_equal(v5, v4, tolerance = 1e-10)
})

test_that("linear world, linear analysis: power with closed-form robust SEs", {
  ## cross-sections, 500 per cell: SE .0412, power .68 for an effect of .10
  d <- did_lin(design = "cross-sections")
  expect_true(powerape:::did_valid(d))
  p <- ape_power(d, n = 2000, claim = "detect", nsim = 600, seed = 51)
  se_cf <- did_se_cs_lin(d$rates, 2000)
  expect_equal(se_cf, 0.0412, tolerance = 5e-3)
  expect_equal(mean(p$draws$se), se_cf, tolerance = 0.01)
  expect_equal(sd(p$draws$est), se_cf, tolerance = 0.08)
  expect_lt(abs(p$est_mean - 0.10), 4 * sd(p$draws$est) / sqrt(600))
  expect_lt(abs(p$power - pnorm(0.10 / se_cf - qnorm(0.975))), 3 * p$mcse + 0.01)
  expect_lt(abs(pnorm(0.10 / se_cf - qnorm(0.975)) - 0.68), 0.005)
  out <- capture.output(print(p))
  expect_true(any(grepl("^  power =", out)))
  expect_true(any(grepl("linear probability \\(OLS\\)", out)))
  expect_true(any(grepl("robust \\(HC1\\) SEs", out)))
  ## `se` is fixed: asking for robust SEs changes nothing and does not warn
  expect_no_warning(pr <- ape_power(d, n = 2000, claim = "detect", nsim = 600,
                                    seed = 51, se = "robust"))
  expect_identical(pr$draws, p$draws)
  ## unequal shares
  du <- did_lin(design = "cross-sections", p_treated = 0.3, p_post = 0.65)
  pu <- ape_power(du, n = 3000, claim = "detect", nsim = 500, seed = 52)
  expect_equal(mean(pu$draws$se), did_se_cs_lin(du$rates, 3000, pt = 0.3, pp = 0.65),
               tolerance = 0.015)
  ## panel: within-unit dependence enters with a minus sign
  for (pt in c(0.5, 0.25)) {
    dp <- did_lin(design = "panel", rho = 0.5, p_treated = pt)
    pp <- ape_power(dp, n = 1200, claim = "detect", nsim = 500, seed = 53)
    se_p <- did_se_panel_lin(dp$rates, 1200, 0.5, pt = pt)
    expect_equal(mean(pp$draws$se), se_p, tolerance = 0.015)
    expect_equal(sd(pp$draws$est), se_p, tolerance = 0.08)
    expect_lt(se_p / did_se_panel_lin(dp$rates, 1200, 0, pt = pt), 0.90)
    expect_lt(abs(pp$est_mean - 0.10), 4 * sd(pp$draws$est) / sqrt(500))
  }
  expect_true(any(grepl("unit-clustered SEs", capture.output(print(pp)))))
})

test_that("the linear/linear cell does not depend on the link", {
  ## without covariates the linear world is four cell rates, whatever G is
  pl <- ape_power(did_lin(model = "logit"), n = 500, claim = "detect", nsim = 100, seed = 61)
  pb <- ape_power(did_lin(model = "probit"), n = 500, claim = "detect", nsim = 100, seed = 61)
  expect_equal(pl$draws$est, pb$draws$est, tolerance = 1e-8)
  expect_equal(pl$draws$se, pb$draws$se, tolerance = 1e-8)
})

test_that("index world, linear analysis: centred on the cross difference, not the ATT", {
  ## no covariates
  d <- did_lin(design = "cross-sections", trends = "index", target = 0.07)
  expect_false(powerape:::did_valid(d))
  cross <- powerape:::did_readings(d)$linear
  expect_gt(cross - 0.07, 0.03)
  p <- ape_power(d, n = 2000, claim = "detect", nsim = 600, seed = 71)
  mc <- sd(p$draws$est) / sqrt(600)
  expect_lt(abs(p$est_mean - cross), 4 * mc)
  expect_gt(abs(p$est_mean - 0.07), 8 * mc)
  out <- capture.output(print(p))
  expect_true(any(grepl("rejection rate", out)))
  expect_true(any(grepl("assumes linear-parallel trends", out)))
  expect_true(any(grepl("not on the ATT", capture.output(print(d)))))
  expect_error(ape_n(d, power = 0.8, claim = "detect"), "valid in the simulated world")
  expect_error(ape_mde(d, n = 2000), "valid in the simulated world")
  ## with a covariate that differs by group, in a panel and in cross-sections:
  ## still the cross difference of the four marginal rates
  for (des in c("panel", "cross-sections")) {
    dc <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                              design = des, rho = if (des == "panel") 0.4 else 0,
                              covariates = list(pa_var("size", "normal")),
                              correlation = 0.4, signal = 0.3,
                              trends = "index", analysis = "linear", n_int = 1e5), 0.08)
    crs <- powerape:::did_readings(dc)$linear
    pc <- ape_power(dc, n = 1500, claim = "detect", nsim = 500, seed = 72)
    expect_lt(abs(pc$est_mean - crs), 4 * sd(pc$draws$est) / sqrt(500))
  }
})

test_that("linear world with covariates: the linear analysis stays centred on the ATT", {
  for (des in c("panel", "cross-sections")) {
    d <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.10, trend = 0.08,
                             design = des, rho = if (des == "panel") 0.4 else 0,
                             covariates = list(pa_var("size", "normal")),
                             correlation = 0.4, signal = 0.3,
                             trends = "linear", analysis = "linear", n_int = 1e5), 0.08)
    p <- ape_power(d, n = 1500, claim = "detect", nsim = 500, seed = 73)
    expect_lt(abs(p$est_mean - 0.08), 4 * sd(p$draws$est) / sqrt(500))
    expect_equal(mean(p$draws$se), sd(p$draws$est), tolerance = 0.10)
  }
})

test_that("linear analysis: no separation, sparse cells still flagged; searches run", {
  ## rare outcome, small treated group: the logit analysis loses studies to
  ## separation, OLS loses none
  mk <- function(an) set_att(ape_dgp_did(model = "logit", baseline = 0.03, p_treated = 0.1,
                                         design = "cross-sections", trends = an,
                                         analysis = an), 0.02)
  pi_ <- suppressWarnings(ape_power(mk("index"), n = 300, claim = "detect", nsim = 200, seed = 4))
  expect_warning(pl <- ape_power(mk("linear"), n = 300, claim = "detect", nsim = 200, seed = 4),
                 "sparsity")
  expect_gt(pi_$separated, 0.05)
  expect_equal(pl$separated, 0)
  expect_lt(pl$outcomes[["failed"]], pi_$outcomes[["failed"]])
  expect_true(!is.na(pl$min_cell) && pl$min_cell < 10)
  ## ape_n against the normal-approximation requirement, and ape_mde back
  d <- did_lin(design = "panel", rho = 0.5)
  c1 <- did_se_panel_lin(d$rates, 1, 0.5)
  n_an <- (c1 * (qnorm(0.975) + qnorm(0.80)) / 0.10)^2
  nn <- ape_n(d, power = 0.80, claim = "detect", nsim = 600, seed = 81, confirm = FALSE)
  expect_lt(abs(nn$n / n_an - 1), 0.12)
  expect_true(any(grepl("linear probability \\(OLS\\), unit-clustered SEs",
                        capture.output(print(nn)))))
  md <- ape_mde(d, n = nn$n, power = 0.80, claim = "detect", nsim = 400, seed = 82,
                confirm = FALSE)
  expect_lt(abs(md$mde / 0.10 - 1), 0.15)
  cv <- ape_curve(d, n = c(300, 900), claim = "detect", nsim = 100, seed = 83)
  expect_null(cv$mismatch)
  expect_identical(cv$model, "linear probability (OLS)")
})

## ---- more than two periods, common timing -----------

# Closed-form SE with n_pre periods before and n_post after, on the pooled
# block rates. Gradient of the estimator in (c11, c10, c01, c00): (1, -1,
# -1, 1) for the linear analysis, the two-by-two index gradient otherwise.
did_grad <- function(analysis, link, r) {
  if (analysis == "linear") return(c(1, -1, -1, 1))
  L <- did_links(link)
  gc <- L$g(L$Gi(r$cf))
  c(1, -gc / L$g(L$Gi(r$c10)), -gc / L$g(L$Gi(r$c01)), gc / L$g(L$Gi(r$c00)))
}
did_se_panel_T <- function(analysis, link, r, n, rho, n_pre, n_post, pt = 0.5) {
  g <- did_grad(analysis, link, r)
  kap <- function(p, q) did_p_joint(p, q, rho) - p * q
  blk <- function(p0, p1, ng, g_post, g_pre) {
    v_post <- (p1 * (1 - p1) + (n_post - 1) * kap(p1, p1)) / (ng * n_post)
    v_pre <- (p0 * (1 - p0) + (n_pre - 1) * kap(p0, p0)) / (ng * n_pre)
    g_post^2 * v_post + g_pre^2 * v_pre + 2 * g_post * g_pre * kap(p0, p1) / ng
  }
  sqrt(blk(r$c10, r$c11, n * pt, g[1], g[2]) +
         blk(r$c00, r$c01, n * (1 - pt), g[3], g[4]))
}
did_se_cs_T <- function(analysis, link, r, n, pt = 0.5, pp = 0.5) {
  g <- did_grad(analysis, link, r)
  v <- function(p, share) p * (1 - p) / (n * share)
  sqrt(g[1]^2 * v(r$c11, pt * pp) + g[2]^2 * v(r$c10, pt * (1 - pp)) +
         g[3]^2 * v(r$c01, (1 - pt) * pp) + g[4]^2 * v(r$c00, (1 - pt) * (1 - pp)))
}
did_T <- function(model = "logit", analysis = "index", design = "panel",
                  n_pre = 2, n_post = 3, rho = 0.5, ...) {
  set_att(ape_dgp_did(model = model, baseline = 0.20, gap = 0.20, trend = 0.10,
                      design = design, rho = if (design == "panel") rho else 0,
                      trends = analysis, analysis = analysis,
                      n_pre = n_pre, n_post = n_post, ...),
          if (analysis == "linear") 0.10 else 0.07)
}

test_that("periods: arguments, layout and prints", {
  expect_error(ape_dgp_did(baseline = 0.3, n_pre = 0), "n_pre >= 1")
  expect_error(ape_dgp_did(baseline = 0.3, n_post = 1.5), "round")
  d2 <- ape_dgp_did(baseline = 0.3)
  expect_identical(d2$n_periods, 2L)
  expect_equal(d2$p_post, 0.5)
  d <- did_T()
  expect_identical(d$n_periods, 5L)
  expect_equal(did_T(design = "cross-sections")$p_post, 3 / 5)
  ## true coefficients: 0 on the second pre-period, l on the three post-periods
  bt <- powerape:::did_beta(d)
  expect_length(bt, 7L)
  expect_equal(bt[3], 0)
  expect_equal(bt[4:6], rep(d$did$l, 3))
  expect_equal(bt[7], d$did$th + d$did$de)
  ## panel: every unit in all five periods, in order
  set.seed(1)
  xx <- powerape:::draw_x_did(d, 50)
  expect_equal(dim(xx$X), c(250L, 7L))
  expect_identical(xx$focal_col, 7L)
  expect_true(all(tabulate(xx$id) == 5L))
  expect_equal(xx$post, rep(c(0, 0, 1, 1, 1), 50))
  expect_equal(rowSums(xx$X[, 3:6]), rep(c(0, 1, 1, 1, 1), 50))   # period 1 is the reference
  expect_equal(xx$X[, 7], xx$d * xx$post)
  ## the two-period layout is [1, D, post, W]
  xx2 <- powerape:::draw_x_did(did_base(), 20)
  expect_equal(unname(xx2$X[, 3]), xx2$post)
  expect_identical(xx2$focal_col, 4L)
  ## cross-sections: 60% after, spread evenly over each block's periods
  set.seed(2)
  xc <- powerape:::draw_x_did(did_T(design = "cross-sections"), 40000)
  expect_equal(mean(xc$post), 0.6, tolerance = 0.02)
  per <- 1 + drop(xc$X[, 3:6] %*% 1:4)
  expect_equal(as.numeric(table(per)) / 40000, c(0.2, 0.2, 0.2, 0.2, 0.2), tolerance = 0.05)
  expect_equal(xc$post, as.numeric(per > 2))
  ## prints
  out <- capture.output(print(d))
  expect_true(any(grepl("5 periods: 2 before, 3 after", out)))
  expect_true(any(grepl("all 5 periods per unit", out)))
  expect_true(any(grepl("two periods", capture.output(print(did_base())))))
  p <- ape_power(d, n = 200, claim = "detect", nsim = 30, seed = 3)
  expect_true(any(grepl("n = 200 units x 5 periods \\(1000 obs\\)", capture.output(print(p)))))
})

test_that("periods: the linear panel estimate is the change in period-averaged group means", {
  d <- did_T(analysis = "linear")
  set.seed(606)
  xx <- powerape:::draw_x_did(d, 200)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  fo <- powerape:::fit_ols(xx$X, y)
  yb <- function(dd, pp) mean(y[xx$d == dd & xx$post == pp])
  did_means <- (yb(1, 1) - yb(1, 0)) - (yb(0, 1) - yb(0, 0))
  expect_equal(fo$coefficients[[xx$focal_col]], did_means, tolerance = 1e-12)
  ## ... and the unit-plus-period fixed-effects estimate
  per <- factor(rep(1:5, 200))
  W <- xx$X[, xx$focal_col]
  fe <- lm(y ~ factor(xx$id) + per + W)
  expect_equal(fo$coefficients[[xx$focal_col]], unname(coef(fe)["W"]), tolerance = 1e-10)
  ## a time-constant covariate stays neutral with five periods
  dz <- set_att(ape_dgp_did(model = "probit", baseline = 0.25, gap = 0.15, trend = 0.08,
                            rho = 0.5, covariates = list(pa_var("x", "normal")),
                            correlation = 0.4, signal = 0.3, trends = "linear",
                            analysis = "linear", n_pre = 2, n_post = 3, n_int = 4e4), 0.08)
  set.seed(607)
  xz <- powerape:::draw_x_did(dz, 200)
  yz <- rbinom(length(xz$pr), 1L, xz$pr)
  expect_equal(powerape:::fit_ols(xz$X, yz)$coefficients[[7]],
               powerape:::fit_ols(xz$X[, 1:7], yz)$coefficients[[7]], tolerance = 1e-12)
})

test_that("periods: the closed forms reduce to the two-period ones", {
  r <- did_base()$rates
  expect_equal(did_se_panel_T("index", "logit", r, 1000, 0.5, 1, 1),
               did_se_panel("logit", r, 1000, 0.5), tolerance = 1e-12)
  expect_equal(did_se_cs_T("index", "probit", r, 2000, 0.3, 0.65),
               did_se_cs("probit", r, 2000, pt = 0.3, pp = 0.65), tolerance = 1e-12)
  rl <- did_lin()$rates
  expect_equal(did_se_panel_T("linear", "logit", rl, 1000, 0.5, 1, 1, pt = 0.25),
               did_se_panel_lin(rl, 1000, 0.5, pt = 0.25), tolerance = 1e-12)
  expect_equal(did_se_cs_T("linear", "logit", rl, 2000), did_se_cs_lin(rl, 2000),
               tolerance = 1e-12)
})

test_that("periods: simulated SEs match the closed forms, both analyses, both structures", {
  for (an in c("linear", "index")) {
    ## panel, two periods before and three after
    d <- did_T(analysis = an)
    p <- ape_power(d, n = 600, claim = "detect", nsim = 400, seed = 91)
    se_cf <- did_se_panel_T(an, "logit", d$rates, 600, 0.5, 2, 3)
    expect_equal(mean(p$draws$se), se_cf, tolerance = 0.015)
    expect_equal(sd(p$draws$est), se_cf, tolerance = 0.10)
    expect_lt(abs(p$est_mean - d$target_est), 4 * sd(p$draws$est) / sqrt(400))
    expect_lt(abs(p$power - pnorm(d$target_est / se_cf - qnorm(0.975))), 3 * p$mcse + 0.015)
    ## cross-sections, three before and one after, 30% treated, probit world
    dc <- did_T(model = "probit", analysis = an, design = "cross-sections",
                n_pre = 3, n_post = 1, p_treated = 0.3)
    pc <- ape_power(dc, n = 3000, claim = "detect", nsim = 400, seed = 92)
    expect_equal(mean(pc$draws$se),
                 did_se_cs_T(an, "probit", dc$rates, 3000, pt = 0.3, pp = 0.25),
                 tolerance = 0.015)
    expect_lt(abs(pc$est_mean - dc$target_est), 4 * sd(pc$draws$est) / sqrt(400))
  }
})

test_that("periods are live: the SE falls with the number of periods as the closed form says", {
  r <- did_lin()$rates
  cf <- vapply(1:3, function(k) did_se_panel_T("linear", "logit", r, 600, 0.5, k, k),
               numeric(1))
  ## one shared unit draw: two periods a side halve the variance
  expect_equal((cf[2] / cf[1])^2, 0.508, tolerance = 5e-3)
  expect_lt(cf[3], cf[2])
  sims <- vapply(1:2, function(k) {
    p <- ape_power(did_T(analysis = "linear", n_pre = k, n_post = k), n = 600,
                   claim = "detect", nsim = 300, seed = 93)
    mean(p$draws$se)
  }, numeric(1))
  expect_equal(sims, cf[1:2], tolerance = 0.015)
  ## and the required n follows
  d <- did_T(analysis = "linear", n_pre = 2, n_post = 2)
  n_an <- (did_se_panel_T("linear", "logit", d$rates, 1, 0.5, 2, 2) *
             (qnorm(0.975) + qnorm(0.80)) / 0.10)^2
  nn <- ape_n(d, power = 0.80, claim = "detect", nsim = 500, seed = 94, confirm = FALSE)
  expect_lt(abs(nn$n / n_an - 1), 0.12)
})

test_that("periods: ATT and clustered SE match marginaleffects with a period factor", {
  skip_if_not_installed("marginaleffects")
  skip_if_not_installed("sandwich")
  d <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.15, trend = 0.08,
                           rho = 0.5, covariates = list(pa_var("x", "normal")),
                           signal = 0.2, n_pre = 2, n_post = 2, n_int = 4e4), 0.08)
  set.seed(808)
  xx <- powerape:::draw_x_did(d, 500)
  dat <- data.frame(y = rbinom(length(xx$pr), 1L, xx$pr), D = xx$d,
                    per = factor(rep(1:4, 500)), W = xx$X[, 6], x = xx$X[, 7],
                    id = xx$id)
  fit <- glm(y ~ D + per + W + x, family = binomial("logit"), data = dat)
  X <- model.matrix(fit)
  expect_equal(unname(X), unname(xx$X), tolerance = 1e-12,
               ignore_attr = TRUE)                     # same design matrix
  est <- powerape:::att_est(coef(fit), X, "logit", 6L)
  V <- powerape:::vcov_cluster(fit, X, dat$id, "logit")
  Vs <- sandwich::vcovCL(fit, cluster = dat$id, type = "HC0", cadjust = TRUE)
  expect_equal(unname(V), unname(Vs), tolerance = 1e-5)
  mfx <- marginaleffects::avg_comparisons(fit, variables = list(W = c(0, 1)),
                                          newdata = subset(dat, W == 1), vcov = Vs)
  expect_equal(est$ape, mfx$estimate, tolerance = 1e-6)
  expect_equal(sqrt(drop(t(est$jac) %*% V %*% est$jac)), mfx$std.error, tolerance = 1e-4)
})

test_that("periods: a period without events fails the probit/logit study, not the OLS one", {
  ## mechanism: one period with no events sends its dummy off; the four
  ## pooled cells do not see it, the backup check on columns 2..T+2 does
  d <- did_T(design = "cross-sections", n_pre = 3, n_post = 3)
  set.seed(909)
  xx <- powerape:::draw_x_did(d, 1200)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  y[xx$X[, 4] == 1] <- 0L                           # period 3: nobody has the outcome
  expect_false(powerape:::cell_check(y, xx$d, xx$post)$sep)
  ft <- powerape:::fit_index_model(xx$X, y, "logit", start = xx$start)
  expect_true(!ft$ok || powerape:::hauck_donner(ft$fit, xx$X, 2:8))
  expect_true(powerape:::fit_ols(xx$X, y)$ok)
  ## in the engine: a rare outcome over ten short periods
  mk <- function(an) set_att(ape_dgp_did(model = "logit", baseline = 0.04,
                                         design = "cross-sections", trends = an,
                                         analysis = an, n_pre = 5, n_post = 5), 0.02)
  pi_ <- suppressWarnings(ape_power(mk("index"), n = 600, claim = "detect", nsim = 150, seed = 5))
  pl <- suppressWarnings(ape_power(mk("linear"), n = 600, claim = "detect", nsim = 150, seed = 5))
  expect_gt(pi_$outcomes[["failed"]], 0.20)         # roughly half the studies have such a period
  expect_true(pi_$separated >= 0.20)
  expect_lt(pl$outcomes[["failed"]], 0.02)
  ## the smallest-cell count is pooled over periods (about 150 x .04 = 6+)
  expect_gt(pi_$min_cell, 3)
})

## ---- a binary subgroup: the difference between two effects (section 17.3) ----

did_sg <- function(model = "logit", analysis = "index", design = "panel",
                   effects = c(0.10, 0.03), share = 0.5, rho = 0.5,
                   baseline = c(0.30, 0.20), gap = c(0.10, 0.15),
                   trend = c(0.10, 0.05), ...) {
  set_att(ape_dgp_did(model = model, baseline = baseline, gap = gap, trend = trend,
                      design = design, rho = if (design == "panel") rho else 0,
                      trends = analysis, analysis = analysis,
                      subgroup = pa_var("female", "binary", p = share), ...),
          effects)
}
# rates of one subgroup as a list the closed-form helpers understand
sg_rates <- function(d, m) lapply(d$rates, function(v) v[m])

test_that("subgroup: arguments, calibration and prints", {
  expect_error(ape_dgp_did(baseline = 0.3, subgroup = "female"), "binary pa_var")
  expect_error(ape_dgp_did(baseline = 0.3, subgroup = pa_var("age", "normal")), "binary pa_var")
  expect_error(ape_dgp_did(baseline = c(0.3, 0.2, 0.1), subgroup = pa_var("f", "binary", p = .5)),
               "one number")
  expect_error(ape_dgp_did(baseline = c(0.3, 0.2)), "length\\(baseline\\) == 1L|num1\\(baseline\\)")
  expect_error(ape_dgp_did(baseline = c(0.3, 0.9), gap = 0.2,
                           subgroup = pa_var("f", "binary", p = .5)), "Subgroup f = 1")
  d0 <- ape_dgp_did(model = "logit", baseline = c(0.30, 0.05), gap = c(0.10, 0.05),
                    trend = 0.10, subgroup = pa_var("female", "binary", p = 0.4))
  expect_equal(d0$rates$c01, c(0.40, 0.15))                  # a scalar trend applies to both
  expect_equal(d0$rates$cf, plogis(qlogis(d0$rates$c10) + qlogis(d0$rates$c01) -
                                     qlogis(d0$rates$c00)), tolerance = 1e-12)
  expect_error(set_att(d0, 0.05), "two effects")
  expect_error(set_att(d0, c(0.05, 0.9)), "infeasible for subgroup female = 1")
  expect_error(ape_power(d0, n = 500, claim = "detect", nsim = 20), "set_att")
  d <- set_att(d0, c(0.10, 0.03))
  expect_identical(d$estimand, "datt")
  expect_equal(d$target_est, -0.07)
  expect_equal(unname(true_att(d)), c(0.10, 0.03, -0.07), tolerance = 1e-9)
  expect_equal(d$rates$c11, d$rates$cf + c(0.10, 0.03), tolerance = 1e-12)
  out <- capture.output(print(d))
  expect_true(any(grepl("with a subgroup", out)))
  expect_true(any(grepl("true difference \\(DATT\\) pinned at -0.0700", out)))
  expect_true(any(grepl("triple cross difference", out)))
  ## covariates related to group and subgroup: every rate holds on its own cell
  dz <- ape_dgp_did(model = "probit", baseline = c(0.25, 0.15), gap = c(0.10, 0.05),
                    trend = c(0.08, 0.04), trends = "linear", rho = 0.3,
                    subgroup = pa_var("female", "binary", p = 0.4),
                    covariates = list(pa_var("size", "normal")), correlation = 0.25,
                    signal = 0.2, n_int = 1e5)
  dz <- set_att(dz, c(0.08, 0.02))
  qc <- powerape:::did_cells_q(dz)
  cc <- dz$did
  for (m in 1:2) {
    expect_equal(mean(pnorm(cc$a[m] + qc[[1]][[m]])), dz$rates$c00[m], tolerance = 1e-7)
    expect_equal(mean(pnorm(cc$a[m] + cc$b[m] + qc[[2]][[m]])), dz$rates$c10[m], tolerance = 1e-7)
    expect_equal(mean(pnorm(cc$a[m] + cc$l[m] + qc[[1]][[m]])), dz$rates$c01[m], tolerance = 1e-7)
    ## linear world: the untreated treated rate moves by the control trend
    expect_equal(mean(pnorm(cc$a[m] + cc$b[m] + cc$l[m] + cc$th[m] + qc[[2]][[m]])),
                 dz$rates$c10[m] + (dz$rates$c01[m] - dz$rates$c00[m]), tolerance = 1e-7)
  }
  expect_equal(unname(true_att(dz)), c(0.08, 0.02, -0.06), tolerance = 1e-7)
})

test_that("subgroup: the same eight rates read on two scales", {
  ## men: control .30 -> .40, treated .40 -> .60; women: .05 -> .15, .10 -> .30
  mk <- function(m) set_att(ape_dgp_did(model = m, baseline = c(0.30, 0.05),
                                        gap = c(0.10, 0.05), trend = c(0.10, 0.10),
                                        trends = "linear",
                                        subgroup = pa_var("female", "binary", p = 0.5)),
                            c(0.10, 0.10))
  rl <- powerape:::did_readings(mk("logit"))
  rp <- powerape:::did_readings(mk("probit"))
  expect_equal(rl$linear, 0, tolerance = 1e-10)             # no difference in points ...
  expect_equal(rl$index, -0.062, tolerance = 2e-2)          # ... six points on the logit index
  expect_equal(rp$index, -0.043, tolerance = 2e-2)
  ## and a difference of zero in points is not a zero interaction coefficient
  d <- set_att(ape_dgp_did(model = "logit", baseline = c(0.30, 0.05), gap = c(0.10, 0.05),
                           trend = c(0.10, 0.10),
                           subgroup = pa_var("female", "binary", p = 0.5)), c(0.08, 0.08))
  expect_equal(d$target_est, 0)
  bt <- powerape:::did_beta(d)
  expect_gt(abs(bt[8]), 0.03)                               # column 2T + 4 = 8
})

test_that("subgroup: design matrix and true coefficients", {
  d <- did_sg(n_pre = 2, n_post = 2)
  set.seed(3)
  xx <- powerape:::draw_x_did(d, 60)
  expect_equal(dim(xx$X), c(240L, 12L))                     # 2 x (T + 2) columns
  expect_identical(xx$w_col, 6L)
  expect_identical(xx$focal_col, 12L)
  expect_equal(unname(xx$X[, 7:12]), unname(xx$m * xx$X[, 1:6]))   # the block times the subgroup
  expect_equal(xx$X[, 12], xx$m * xx$d * xx$post)
  expect_true(all(tapply(xx$m, xx$id, function(v) length(unique(v))) == 1L))
  bt <- powerape:::did_beta(d)
  cc <- d$did
  expect_equal(bt[1:6], c(cc$a[1], cc$b[1], 0, cc$l[1], cc$l[1], cc$th[1] + cc$de[1]))
  expect_equal(bt[7:12], c(cc$a[2], cc$b[2], 0, cc$l[2], cc$l[2], cc$th[2] + cc$de[2]) - bt[1:6])
  ## the subgroup share is live
  xs <- powerape:::draw_x_did(did_sg(share = 0.25), 20000)
  expect_equal(mean(xs$m), 0.25, tolerance = 0.05)
  ## cell rates of both subgroups come out of the draw
  dd <- did_sg(design = "cross-sections")
  set.seed(4)
  xc <- powerape:::draw_x_did(dd, 200000)
  y <- rbinom(length(xc$pr), 1L, xc$pr)
  for (m in 0:1)
    expect_equal(mean(y[xc$m == m & xc$d == 1 & xc$post == 1]), dd$rates$c11[m + 1],
                 tolerance = 0.03)
})

test_that("subgroup: the pooled estimate is the difference of two separate estimates", {
  for (m in c("logit", "probit")) {
    d <- did_sg(model = m, n_pre = 1, n_post = 2)
    set.seed(11)
    xx <- powerape:::draw_x_did(d, 600)
    y <- rbinom(length(xx$pr), 1L, xx$pr)
    ft <- powerape:::fit_index_model(xx$X, y, m, start = xx$start)
    est <- powerape:::datt_est(ft$fit$coefficients, xx, m)
    ## separate fits of [1, D, f_2, f_3, W] within each subgroup
    sep <- vapply(0:1, function(s) {
      r <- xx$m == s
      f <- powerape:::fit_index_model(xx$X[r, 1:5], y[r], m)
      powerape:::att_est(f$fit$coefficients, xx$X[r, 1:5], m, 5L)$ape
    }, numeric(1))
    expect_equal(est$by, sep, tolerance = 1e-6)
    expect_equal(est$ape, sep[2] - sep[1], tolerance = 1e-6)
    ## Jacobian of the difference by finite differences
    b <- ft$fit$coefficients
    fd <- vapply(seq_along(b), function(j) {
      h <- 1e-6
      bp <- b; bp[j] <- bp[j] + h
      bm <- b; bm[j] <- bm[j] - h
      (powerape:::datt_est(bp, xx, m)$ape - powerape:::datt_est(bm, xx, m)$ape) / (2 * h)
    }, numeric(1))
    expect_equal(unname(est$jac), fd, tolerance = 1e-6)
    ## the clustered variance of the difference is the sum of the two
    ## separate ones, up to the G / (G - 1) factors
    V <- powerape:::vcov_cluster(ft$fit, xx$X, xx$id, m)
    v_pool <- drop(t(est$jac) %*% V %*% est$jac)
    v_sep <- sum(vapply(0:1, function(s) {
      r <- xx$m == s
      f <- powerape:::fit_index_model(xx$X[r, 1:5], y[r], m)
      e <- powerape:::att_est(f$fit$coefficients, xx$X[r, 1:5], m, 5L)
      drop(t(e$jac) %*% powerape:::vcov_cluster(f$fit, xx$X[r, 1:5], xx$id[r], m) %*% e$jac)
    }, numeric(1)))
    expect_equal(v_pool, v_sep, tolerance = 0.02)
  }
  ## linear analysis: the coefficient is the difference of the two
  ## differences in period-averaged group changes
  d <- did_sg(analysis = "linear", n_pre = 1, n_post = 2)
  set.seed(12)
  xx <- powerape:::draw_x_did(d, 400)
  y <- rbinom(length(xx$pr), 1L, xx$pr)
  yb <- function(s, dd, pp) mean(y[xx$m == s & xx$d == dd & xx$post == pp])
  dd_s <- function(s) (yb(s, 1, 1) - yb(s, 1, 0)) - (yb(s, 0, 1) - yb(s, 0, 0))
  expect_equal(powerape:::fit_ols(xx$X, y)$coefficients[[xx$focal_col]], dd_s(1) - dd_s(0),
               tolerance = 1e-12)
  ## no study without treated post-period rows in both subgroups
  xx$X[, xx$focal_col] <- 0
  expect_null(powerape:::datt_est(rep(0, ncol(xx$X)), xx, "logit"))
})

test_that("subgroup: simulated SEs equal the two closed forms combined", {
  comb <- function(an, link, d, n, share, design) {
    s <- c(1 - share, share)
    v <- vapply(1:2, function(m) {
      r <- sg_rates(d, m)
      if (design == "panel") did_se_panel_T(an, link, r, n * s[m], 0.5, 1, 1)^2
      else did_se_cs_T(an, link, r, n * s[m])^2
    }, numeric(1))
    sqrt(sum(v))
  }
  for (an in c("index", "linear")) {
    d <- did_sg(analysis = an)
    p <- ape_power(d, n = 2400, claim = "detect", nsim = 400, seed = 21)
    se_cf <- comb(an, "logit", d, 2400, 0.5, "panel")
    expect_equal(mean(p$draws$se), se_cf, tolerance = 0.015)
    expect_equal(sd(p$draws$est), se_cf, tolerance = 0.10)
    expect_lt(abs(p$est_mean - (-0.07)), 4 * sd(p$draws$est) / sqrt(400))
    expect_identical(p$estimand, "datt")
    expect_true(any(grepl("assumed true DATT -0.0700", capture.output(print(p)))))
  }
  ## a quarter of the units in the subgroup, cross-sections, probit
  dq <- did_sg(model = "probit", design = "cross-sections", share = 0.25)
  pq <- ape_power(dq, n = 6000, claim = "detect", nsim = 300, seed = 22)
  expect_equal(mean(pq$draws$se), comb("index", "probit", dq, 6000, 0.25, "cross-sections"),
               tolerance = 0.02)
  ## cost: with equal subgroups and the same rates, a difference as large
  ## as the effect needs about four times the units of the effect
  de <- did_sg(analysis = "linear", effects = c(0, 0.10), baseline = 0.20, gap = 0.20,
               trend = 0.10)
  one <- did_se_panel_T("linear", "logit", sg_rates(de, 2), 1, 0.5, 1, 1)^2
  both <- comb("linear", "logit", de, 1, 0.5, "panel")^2
  expect_gt(both / one, 3.8)
  expect_lt(both / one, 4.2)
})

test_that("subgroup: mismatched analysis, searches, sparse cells", {
  ## index world, linear analysis: centred on the triple cross difference
  dm <- set_att(ape_dgp_did(model = "logit", baseline = c(0.30, 0.05), gap = c(0.10, 0.05),
                            trend = c(0.10, 0.10), rho = 0.5, analysis = "linear",
                            subgroup = pa_var("female", "binary", p = 0.5)),
                c(0.091, 0.029))
  expect_false(powerape:::did_valid(dm))
  crs <- powerape:::did_readings(dm)$linear
  expect_lt(abs(crs), 0.005)
  p <- ape_power(dm, n = 2000, claim = "detect", nsim = 300, seed = 31)
  expect_lt(abs(p$est_mean - crs), 4 * sd(p$draws$est) / sqrt(300))
  out <- capture.output(print(p))
  expect_true(any(grepl("rejection rate", out)))
  expect_true(any(grepl("the true DATT is -0.0620", out)))
  expect_error(ape_n(dm, power = 0.8, claim = "detect"), "valid in the simulated world")
  ## searches in a valid design
  d <- did_sg(analysis = "linear", effects = c(0.02, 0.12), baseline = 0.20, gap = 0.20,
              trend = 0.10)
  c2 <- sum(vapply(1:2, function(m) did_se_panel_T("linear", "logit", sg_rates(d, m), 0.5, 0.5, 1, 1)^2,
                   numeric(1)))
  n_an <- c2 * (qnorm(0.975) + qnorm(0.80))^2 / 0.10^2
  nn <- ape_n(d, power = 0.80, claim = "detect", nsim = 500, seed = 32, confirm = FALSE)
  expect_lt(abs(nn$n / n_an - 1), 0.12)
  expect_identical(nn$estimand, "datt")
  md <- ape_mde(d, n = nn$n, power = 0.80, claim = "detect", nsim = 300, seed = 33,
                confirm = FALSE)
  expect_identical(md$estimand, "datt")
  expect_lt(abs(md$mde / 0.10 - 1), 0.20)
  expect_equal(md$dgp$att_sub[1], 0.02)                    # the reference effect is held
  expect_true(any(grepl("minimum detectable DATT", capture.output(print(md)))))
  expect_error(ape_mde(ape_dgp_did(baseline = 0.2, subgroup = pa_var("f", "binary", p = .5)),
                       n = 500), "pin the design first")
  ## equivalence: no meaningful difference between the subgroups
  d0 <- did_sg(analysis = "linear", effects = c(0.08, 0.08))
  pe <- ape_power(d0, n = 6000, claim = "equivalence", sesoi = 0.10, nsim = 150, seed = 34)
  expect_gt(pe$power, 0.5)
  ## the identifying cells are checked within each subgroup
  xx <- list(d = rep(c(0, 1), each = 40), post = rep(c(0, 1), 40),
             m = rep(c(0, 0, 1, 1), times = 20))
  y <- rep(c(0L, 1L), length.out = 80)
  y[] <- as.integer(seq_len(80) %% 3 == 0)
  ok <- powerape:::did_cell_check(y, xx)
  y2 <- y
  y2[xx$m == 1 & xx$d == 1 & xx$post == 1] <- 0L          # one subgroup cell without events
  expect_false(ok$sep)
  expect_true(powerape:::did_cell_check(y2, xx)$sep)
  expect_false(powerape:::cell_check(y2, xx$d, xx$post)$sep)   # invisible to the pooled cells
})
