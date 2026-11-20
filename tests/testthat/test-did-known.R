# Known-answer checks of the difference-in-differences route: cases where
# the answer is fixed by an identity, by another route of the package, or
# by exact enumeration.

test_that("the scales agree when there is no baseline gap or no control-group change", {
  ## with gap = 0 or trend = 0, "parallel" means the same on every scale:
  ## the linear and the index reading of the same four rates coincide
  for (lk in c("probit", "logit")) for (s in list(c(0, 0.10), c(0.20, 0))) {
    d <- set_att(ape_dgp_did(model = lk, baseline = 0.20, gap = s[1], trend = s[2],
                             design = "cross-sections", trends = "linear"), 0.10)
    rd <- powerape:::did_readings(d)
    expect_equal(rd$linear, 0.10, tolerance = 1e-10)
    expect_equal(rd$index, 0.10, tolerance = 1e-10)
    ## and the two worlds are the same world: identical rates
    di <- set_att(ape_dgp_did(model = lk, baseline = 0.20, gap = s[1], trend = s[2],
                              design = "cross-sections", trends = "index"), 0.10)
    expect_equal(unlist(di$rates), unlist(d$rates), tolerance = 1e-10)
  }
  ## with a gap and a trend they differ
  d <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0.20, trend = 0.10,
                           design = "cross-sections", trends = "linear"), 0.10)
  expect_gt(abs(powerape:::did_readings(d)$index - 0.10), 0.02)
  ## so an analysis on the "wrong" scale is unbiased when there is no gap
  dm <- set_att(ape_dgp_did(model = "logit", baseline = 0.20, gap = 0, trend = 0.10,
                            design = "cross-sections", trends = "linear", analysis = "index"), 0.10)
  p <- suppressWarnings(ape_power(dm, 2000, "detect", nsim = 800, seed = 1))
  expect_lt(abs(p$est_mean - 0.10), 4 * p$est_mcse)
})

test_that("cross-sections in a linear world are the interaction route's design", {
  ## the effect on the treated is the cross difference of the four rates:
  ## the two-way interaction effect with focal = group, moderator = period
  n <- 1500
  dd <- set_att(ape_dgp_did(model = "probit", baseline = 0.20, gap = 0.20, trend = 0.10,
                            design = "cross-sections", trends = "linear", analysis = "linear"), 0.10)
  da <- set_aie(ape_dgp(model = "probit", focal = pa_var("group", "binary", p = 0.5),
                        moderator = pa_var("post", "binary", p = 0.5), baseline = 0.20),
                target = 0.10, main_focal = 0.20, main_moderator = 0.10)
  pc <- c(0.20, 0.30, 0.40, 0.60)
  expect_equal(unname(unlist(dd$rates[c("c00", "c01", "c10", "c11")])), pc, tolerance = 1e-10)
  p_cf <- pnorm(0.10 / sqrt(sum(pc * (1 - pc) / (n * 0.25))) - qnorm(0.975))
  near <- function(p) expect_lt(abs(p$power - p_cf), 3.5 * p$mcse)
  near(suppressWarnings(ape_power(dd, n, "detect", nsim = 3000, seed = 2)))
  near(suppressWarnings(ape_power(set_analysis(da, "linear"), n, "detect", nsim = 3000, seed = 3)))
  near(suppressWarnings(ape_power(da, n, "detect", nsim = 1500, seed = 4)))
})

test_that("exact power by enumeration: both analyses, both links", {
  ## n = 80 in four cells, rates .30/.40/.40 and an effect of .30
  n <- 80
  for (cs in list(c("logit", "index"), c("probit", "index"), c("logit", "linear"))) {
    lk <- cs[1]; an <- cs[2]
    d <- set_att(ape_dgp_did(model = lk, baseline = 0.30, gap = 0.10, trend = 0.10,
                             design = "cross-sections", trends = an, analysis = an), 0.30)
    pc <- unname(unlist(d$rates[c("c00", "c01", "c10", "c11")]))
    ex <- exact_did_xs(pc, n, lk, an, draws = 300L, seed = 7L)
    p <- suppressWarnings(ape_power(d, n, "detect", nsim = 4000, seed = 8))
    expect_lt(abs(p$power - ex$power), 3.5 * sqrt(p$mcse^2 + ex$size_se^2))
    ## the share of studies without a usable estimate: the sharper check of
    ## the separation rule (none under the linear analysis)
    fs <- p$outcomes[["failed"]]
    expect_lt(abs(fs - ex$failed),
              3.5 * sqrt(max(fs, ex$failed, 1e-4) * (1 - fs) / 4000 + ex$failed_se^2) + 1e-4)
    if (an == "linear") expect_lt(ex$failed, 1e-6) else expect_gt(ex$failed, 0.005)
  }
})

test_that("a design whose two scales coincide is valid under either analysis", {
  ## no gap between the groups: the linear world and the index world are one
  ## world, so the logit analysis is valid in the "linear" world
  mk <- function(tr, an = "index", ...)
    set_att(ape_dgp_did(model = "logit", baseline = 0.20, trend = 0.10,
                        design = "cross-sections", trends = tr, analysis = an, ...), 0.08)
  dl <- mk("linear", gap = 0)
  di <- mk("index", gap = 0)
  expect_true(powerape:::did_valid(dl))
  pl <- suppressWarnings(ape_power(dl, 1500, "detect", nsim = 200, seed = 4))
  pi <- suppressWarnings(ape_power(di, 1500, "detect", nsim = 200, seed = 4))
  expect_identical(pl$power, pi$power)
  expect_false(any(grepl("rejection rate", capture.output(print(pl)))))
  expect_true(any(grepl("mean the same on both scales", capture.output(print(dl)))))
  n <- suppressWarnings(ape_n(dl, claim = "detect", nsim = 200, seed = 5))
  expect_gt(n$n, 1000)
  ## no change in the control group, and the linear analysis in an index world
  expect_true(powerape:::did_valid(set_att(ape_dgp_did(
    model = "probit", baseline = 0.20, gap = 0.15, trend = 0,
    trends = "index", analysis = "linear"), 0.08)))
  ## with a gap and a trend, or with covariates that move the outcome, the
  ## mismatch stays a mismatch
  dg <- mk("linear", gap = 0.20)
  expect_false(powerape:::did_valid(dg))
  expect_error(suppressWarnings(ape_n(dg, claim = "detect", nsim = 100, seed = 1)),
               "biased for the effect on the treated")
  dc <- mk("linear", gap = 0, covariates = list(pa_var("z", "normal")),
           signal = 0.3, n_int = 1e5)
  expect_false(powerape:::did_valid(dc))
  expect_true(any(grepl("centres near that covariate-free reading",
                        capture.output(print(dc)))))
  ## a subgroup design: the condition holds in each subgroup
  ds <- set_att(ape_dgp_did(model = "logit", baseline = c(0.2, 0.3), gap = 0,
                            trend = c(0.1, 0.05), trends = "linear",
                            subgroup = pa_var("female", "binary", p = 0.4)), c(0.08, 0.03))
  expect_true(powerape:::did_valid(ds))
  ds2 <- set_att(ape_dgp_did(model = "logit", baseline = c(0.2, 0.3), gap = c(0, 0.1),
                             trend = c(0.1, 0.05), trends = "linear",
                             subgroup = pa_var("female", "binary", p = 0.4)), c(0.08, 0.03))
  expect_false(powerape:::did_valid(ds2))
})

test_that("inputs a design does not use are named once, and argument errors are worded", {
  expect_warning(d <- ape_dgp_did(baseline = 0.3, gap = 0.1, trend = 0.05, p_post = 0.4),
                 "`p_post` is ignored")
  expect_null(d$spec$p_post)
  expect_warning(ape_dgp_did(baseline = 0.3, gap = 0.1, trend = 0.05, correlation = 0.3),
                 "`correlation` is ignored")
  ## a rebuilt copy stays silent: the sweep rebuilds every scenario
  expect_no_warning(ape_robust(set_att(d, 0.06), n = 3000, claim = "detect",
                               vary = list(rho = c(0, 0.3)), nsim = 20, seed = 1,
                               nmax = FALSE))
  expect_no_warning(ape_dgp_did(baseline = 0.3, gap = 0.1, trend = 0.05,
                                design = "cross-sections", p_post = 0.4))
  expect_error(ape_dgp_did(baseline = c(0.3, 0.2), gap = 0.1),
               "take two numbers only in a design with a `subgroup`")
  expect_error(set_att(d, c(0.05, 0.02)), "`target` is one number for this design")
})

test_that("the pilot builder orders periods by their type and keeps quiet", {
  set.seed(5)
  N <- 600
  pil <- data.frame(id = rep(1:N, each = 3), per = rep(1:3, N),
                    treated = rep(rbinom(N, 1, 0.4), each = 3))
  pil$y <- rbinom(nrow(pil), 1, 0.2 + 0.1 * pil$treated + 0.05 * (pil$per - 1))
  b <- ape_dgp_did_pilot(pil, "y", "treated", id = "id", period = "per", trend = 0.05)
  lab <- c("c_first", "b_second", "a_third")
  alt <- list(date = as.Date("2020-01-01") + 365 * (pil$per - 1),
              numtext = c("9", "10", "11")[pil$per],
              fac = factor(lab[pil$per], levels = lab))
  for (nm in names(alt)) {
    q <- pil
    q$per <- alt[[nm]]
    expect_no_warning(o <- ape_dgp_did_pilot(q, "y", "treated", id = "id",
                                             period = "per", trend = 0.05))
    expect_equal(o$pilot$control_change, b$pilot$control_change, info = nm)
    expect_equal(o$rho, b$rho, info = nm)
  }
  ## other text is ordered alphabetically, which reverses these labels
  q <- pil
  q$per <- lab[pil$per]
  expect_no_warning(o <- ape_dgp_did_pilot(q, "y", "treated", id = "id",
                                           period = "per", trend = 0.05))
  expect_equal(o$pilot$control_change, -b$pilot$control_change)
  ## a unit identifier with repeated cross-sections asked for: no warning
  ## about a rho the user did not give
  expect_no_warning(x <- ape_dgp_did_pilot(pil, "y", "treated", id = "id", period = "per",
                                           trend = 0.05, design = "cross-sections"))
  expect_identical(x$design, "cross-sections")
  expect_equal(x$rho, 0)
})
