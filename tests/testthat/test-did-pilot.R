# ape_dgp_did_pilot(): the inputs of a difference-in-differences design
# estimated from pre-treatment data.

# A large pre-treatment panel drawn from the route itself.
pilot_from_route <- function(n = 20000, rho = 0.5, n_pre = 3) {
  d0 <- set_att(ape_dgp_did(model = "logit", baseline = 0.25, gap = 0.10, trend = 0.05,
                            rho = rho, n_pre = n_pre, n_post = 1), 0.05)
  powerape:::with_seed(41, {
    xx <- powerape:::draw_x(d0, n)
    y <- rbinom(length(xx$pr), 1L, xx$pr)
    per <- rep(seq_len(n_pre + 1L), n)
    keep <- per <= n_pre
    data.frame(y = y[keep], treated = xx$X[keep, 2], firm = xx$id[keep], year = per[keep])
  })
}

test_that("the joint rate and the rho solver are exact", {
  expect_equal(powerape:::did_joint_rate(0.3, 0), 0.09)
  ## increasing in rho, approaching p as rho -> 1
  jr <- vapply(c(0.1, 0.5, 0.9, 0.99), function(r) powerape:::did_joint_rate(0.3, r), numeric(1))
  expect_true(all(diff(jr) > 0))
  expect_lt(abs(powerape:::did_joint_rate(0.3, 0.995) - 0.3), 0.02)
  ## against a direct simulation of the copula
  set.seed(42)
  a <- rnorm(4e5)
  pc <- pnorm((qnorm(0.3) - sqrt(0.5) * a) / sqrt(0.5))
  expect_equal(powerape:::did_joint_rate(0.3, 0.5), mean(pc^2), tolerance = 5e-3)
  ## the solver returns the rho that generated the joint rates
  p <- c(0.25, 0.35)
  for (r0 in c(0.2, 0.5, 0.8)) {
    m <- vapply(p, powerape:::did_joint_rate, numeric(1), rho = r0)
    expect_equal(powerape:::did_pilot_rho(p, m, c(100, 100)), r0, tolerance = 1e-6)
  }
  ## no positive dependence in the data: 0
  expect_identical(powerape:::did_pilot_rho(p, p^2 - 0.01, c(100, 100)), 0)
})

test_that("the inputs are recovered from a large pilot drawn from the route", {
  pl <- pilot_from_route()
  d <- ape_dgp_did_pilot(pl, outcome = "y", group = "treated", id = "firm", period = "year",
                         trend = 0.05, model = "logit")
  expect_s3_class(d, "powerape_dgp")
  expect_identical(d$design, "panel")
  expect_lt(abs(d$baseline - 0.25), 0.01)
  expect_lt(abs(d$gap - 0.10), 0.012)
  expect_lt(abs(d$p_treated - 0.5), 0.012)
  expect_lt(abs(d$rho - 0.5), 0.03)
  expect_equal(d$pilot$n_rows, 60000)
  expect_equal(d$pilot$n_units, 20000)
  ## rates are flat before treatment in the route, so the change between
  ## pilot periods is near zero
  expect_lt(abs(d$pilot$control_change), 0.01)
  expect_true(any(grepl("pilot: baseline, gap, treated share, rho estimated from 60000 pre-treatment rows (20000 units)",
                        capture.output(print(d)), fixed = TRUE)))
  ## the design works downstream, and the paragraph says where inputs came from
  d <- set_att(d, 0.08)
  pw <- suppressWarnings(ape_power(d, 1200, "detect", nsim = 60, seed = 1))
  expect_true(pw$power > 0 && pw$power < 1)
  expect_true(grepl("estimated from 60000 pre-treatment pilot observations", power_statement(pw), fixed = TRUE))
  ## and it can be swept: the rebuilt scenarios come from the stored inputs
  r <- suppressWarnings(ape_robust(d, n = 1200, claim = "detect", nsim = 40, seed = 2,
                                   vary = list(rho = c(0.3, 0.7)), grid_points = 2, nmax = FALSE))
  expect_equal(nrow(r$scenarios), 2L)
})

test_that("overrides, cross-sections and guards", {
  pl <- pilot_from_route(n = 3000)
  d <- ape_dgp_did_pilot(pl, "y", "treated", id = "firm", trend = 0.05, baseline = 0.20, model = "logit")
  expect_equal(d$baseline, 0.20)
  expect_identical(d$pilot$overridden, "baseline")
  expect_false(grepl("baseline", powerape:::did_pilot_line(d), fixed = TRUE))
  ## without a unit identifier: repeated cross-sections, no rho
  dx <- ape_dgp_did_pilot(pl, "y", "treated", trend = 0.05, model = "logit")
  expect_identical(dx$design, "cross-sections")
  expect_true(is.na(dx$pilot$rho))
  ## one period per unit: rho cannot be estimated
  one <- pl[pl$year == 1, ]
  expect_error(ape_dgp_did_pilot(one, "y", "treated", id = "firm", trend = 0.05), "supply `rho`")
  expect_equal(ape_dgp_did_pilot(one, "y", "treated", id = "firm", trend = 0.05, rho = 0.3)$rho, 0.3)
  ## trend is the user's to give
  expect_error(ape_dgp_did_pilot(pl, "y", "treated", id = "firm", period = "year"),
               "`trend` is required")
  expect_error(ape_dgp_did_pilot(pl, "y", "treated", id = "firm", period = "year"),
               "control group's rate changed by")
  expect_error(ape_dgp_did_pilot(pl, "y", "nope", trend = 0.05), "Not in `data`")
  bad <- pl; bad$y[1] <- 2
  expect_error(ape_dgp_did_pilot(bad, "y", "treated", trend = 0.05), "must be binary")
  bad <- pl; bad$treated[1] <- 1 - bad$treated[1]
  expect_error(ape_dgp_did_pilot(bad, "y", "treated", id = "firm", trend = 0.05), "constant within a unit")
  expect_error(ape_dgp_did_pilot(pl, "y", "treated", trend = 0.05,
                                 subgroup = pa_var("m", "binary", p = 0.5)), "subgroup")
  expect_warning(ape_dgp_did_pilot(pl[1:300, ], "y", "treated", trend = 0.05), "fewer than 500 rows")
})

# Covariates and rho ------------------------------------------------------------
# Covariates are constant over time and add within-unit dependence of their
# own; the design's rho is the latent correlation given the covariates.

pilot_with_covariate <- function(model, n = 20000, rho = 0.2, signal = 0.3, seed = 43) {
  d0 <- set_att(ape_dgp_did(model = model, baseline = 0.25, gap = 0.10, trend = 0.05,
                            rho = rho, n_pre = 2, n_post = 1,
                            covariates = list(pa_var("z", "normal")), signal = signal,
                            n_int = 1e5), 0.05)
  powerape:::with_seed(seed, {
    xx <- powerape:::draw_x(d0, n)
    y <- rbinom(length(xx$pr), 1L, xx$pr)
    per <- rep(1:3, n)
    keep <- per <= 2
    data.frame(y = y[keep], treated = xx$d[keep], firm = xx$id[keep], year = per[keep])
  })
}

test_that("with covariates, rho is the dependence given the covariates", {
  ## probit with a normal covariate: the raw latent correlation of two
  ## outcomes of a unit is (rho + s2) / (1 + s2), s2 = signal / (1 - signal)
  pl <- pilot_with_covariate("probit")
  Y <- matrix(pl$y, ncol = 2, byrow = TRUE)
  ctrl <- pl$treated[pl$year == 1] == 0
  s2 <- 0.3 / 0.7
  jr <- mean(Y[ctrl, 1] * Y[ctrl, 2])
  expect_lt(abs(jr - powerape:::did_joint_rate(0.25, (0.2 + s2) / (1 + s2))), 0.012)
  expect_gt(jr, powerape:::did_joint_rate(0.25, 0.2) + 0.015)
  ## the label says so, and only when covariates move the outcome
  dz <- ape_dgp_did(baseline = 0.25, rho = 0.2, covariates = list(pa_var("z", "normal")),
                    signal = 0.3, n_int = 4e4)
  expect_true(any(grepl("rho = 0.20 (latent within-unit correlation given the covariates)",
                        capture.output(print(dz)), fixed = TRUE)))
  d0 <- ape_dgp_did(baseline = 0.25, rho = 0.2)
  expect_true(any(grepl("rho = 0.20 (latent within-unit correlation);",
                        capture.output(print(d0)), fixed = TRUE)))
  dn <- ape_dgp_did(baseline = 0.25, rho = 0.2, covariates = list(pa_var("z", "normal")),
                    signal = 0, n_int = 4e4)
  expect_false(any(grepl("given the covariates", capture.output(print(dn)), fixed = TRUE)))
})

test_that("the joint rate with covariates: limits, grid and simulation", {
  lf <- powerape:::link_funs("probit")
  ## no covariate index: the plain joint rate
  for (r in c(0, 0.3, 0.8))
    expect_equal(powerape:::did_joint_rate_cov(lf, qnorm(0.3), 0, r),
                 powerape:::did_joint_rate(0.3, r), tolerance = 1e-8)
  ## the quantile grid reproduces averages over the draws
  set.seed(44)
  q <- 0.7 * rnorm(2e5)
  g <- powerape:::did_q_grid(q)
  expect_length(g, 2000L)
  expect_equal(mean(pnorm(-0.6 + g)^2), mean(pnorm(-0.6 + q)^2), tolerance = 1e-4)
  expect_identical(powerape:::did_q_grid(1:10), 1:10)
  ## against a direct simulation, logit link
  ll <- powerape:::link_funs("logit")
  a <- rnorm(2e5)
  p <- plogis(-1 + q)
  pc <- pnorm((qnorm(p) - sqrt(0.4) * a) / sqrt(0.6))
  expect_equal(powerape:::did_joint_rate_cov(ll, -1, q[1:20000], 0.4), mean(pc^2),
               tolerance = 0.01)
})

test_that("a pilot with covariates returns the rho given the covariates", {
  s2 <- 0.3 / 0.7
  for (model in c("probit", "logit")) {
    pl <- pilot_with_covariate(model, n = 60000)
    d <- ape_dgp_did_pilot(pl, "y", "treated", id = "firm", period = "year", trend = 0.05,
                           model = model, n_pre = 2, n_post = 1,
                           covariates = list(pa_var("z", "normal")), signal = 0.3,
                           n_int = 1e5)
    expect_lt(abs(d$rho - 0.2), 0.04)
    expect_identical(d$rho, d$pilot$rho)
    expect_identical(d$spec$rho, d$rho)
    expect_gt(d$pilot$rho_raw, d$rho + 0.05)
    if (model == "probit") {
      ## closed form under probit with a normal covariate
      expect_lt(abs(d$pilot$rho_raw - (0.2 + s2) / (1 + s2)), 0.03)
      expect_lt(abs(d$rho - (d$pilot$rho_raw - 0.3) / 0.7), 0.005)
    }
    ## the design equals one built directly with that rho
    dd <- ape_dgp_did(model = model, baseline = d$baseline, gap = d$gap, trend = 0.05,
                      p_treated = d$p_treated, design = "panel", rho = d$rho,
                      n_pre = 2, n_post = 1, covariates = list(pa_var("z", "normal")),
                      signal = 0.3, n_int = 1e5)
    same <- setdiff(names(dd), "spec")
    expect_equal(d[same], dd[same])
    expect_true(any(grepl("rho given the covariates", capture.output(print(d)), fixed = TRUE)))
    expect_true(grepl("the within-unit correlation given the covariates estimated from",
                      powerape:::did_pilot_clause(d), fixed = TRUE))
    ## without the covariates the raw correlation is the design's rho
    d0 <- ape_dgp_did_pilot(pl, "y", "treated", id = "firm", period = "year", trend = 0.05,
                            model = model, n_pre = 2, n_post = 1)
    expect_identical(d0$rho, d0$pilot$rho_raw)
    expect_equal(d0$rho, d$pilot$rho_raw)
    expect_false(any(grepl("given the covariates", capture.output(print(d0)), fixed = TRUE)))
  }
})

test_that("pilot with covariates: no conversion without signal or with a given rho; too strong a signal warns", {
  pl <- pilot_with_covariate("probit", n = 6000)
  cv <- list(pa_var("z", "normal"))
  d0 <- ape_dgp_did_pilot(pl, "y", "treated", id = "firm", trend = 0.05,
                          covariates = cv, signal = 0, n_int = 4e4)
  expect_identical(d0$rho, d0$pilot$rho_raw)
  dr <- ape_dgp_did_pilot(pl, "y", "treated", id = "firm", trend = 0.05, rho = 0.33,
                          covariates = cv, signal = 0.3, n_int = 4e4)
  expect_equal(dr$rho, 0.33)
  expect_false(powerape:::did_pilot_converted(dr$pilot))
  ## a signal that alone implies more dependence than the pilot shows
  expect_warning(
    dw <- ape_dgp_did_pilot(pl, "y", "treated", id = "firm", trend = 0.05,
                            covariates = cv, signal = 0.7, n_int = 4e4),
    "more within-unit dependence than")
  expect_identical(dw$rho, 0)
})
