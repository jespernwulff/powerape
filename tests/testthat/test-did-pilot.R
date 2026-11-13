# ape_dgp_did_pilot(): the inputs of a difference-in-differences design
# estimated from pre-treatment data (DESIGN.md section 20).

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
