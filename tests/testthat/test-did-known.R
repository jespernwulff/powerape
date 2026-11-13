# Known-answer checks of the difference-in-differences route: cases where
# the answer is fixed by an identity, by another route of the package, or
# by exact enumeration (DESIGN.md section 20).

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
