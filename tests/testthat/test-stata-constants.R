# Regression tests against Stata for the linear probability analysis, the
# difference-in-differences route (both analyses, with a subgroup) and the
# three-way interaction effect. The constants were computed with StataNow
# 19.5 (revision 14 Oct 2026) on 2026-11-11, with official commands only:
# regress, logit, probit, margins, didregress, xtdidregress. The data are
# the fixtures of helper-stata-fixtures.R; the first test confirms they are
# the ones Stata read. Where a standard error is compared the fixture is a
# logit, because Stata's probit variance uses the observed information and
# the package's the expected (the two coincide for logit).

st <- c(
  sum_lin_y = 244, sum_lin_z = 1.24998766830473e+02,
  sum_did_y = 323, sum_did_z = -3.24264721683877e+01,
  sum_didx_y = 528, sum_didx_period = 3024,
  sum_sub_y = 501, sum_sub_m = 482,
  sum_tw_y_bb = 612, sum_tw_y_cc = 604,
  sum_tw_c1 = 1.62427608622864e+03, sum_tw_x = -3.54411219846858e+01,
  lin_hc1_b = 1.99143782218859e-01, lin_hc1_se = 4.16081438565490e-02,
  lin_cl_b = 2.85699826204525e-01, lin_cl_se = 7.56549784041522e-02,
  did_att = 3.75498231461815e-02, did_att_se = 6.70190830682530e-02,
  did_lin_b = 5.24839743589708e-02, did_lin_se = 6.16605663659612e-02,
  did_xtdid_b = 5.24839743589742e-02, did_xtdid_se = 6.15831519394962e-02,
  didx_att = 9.04580826660472e-02, didx_att_se = 5.93718463221340e-02,
  didx_lin_b = 1.03827792495657e-01, didx_lin_se = 5.56246674240080e-02,
  didx_didreg_b = 1.03827792495657e-01, didx_didreg_se = 5.56246674240080e-02,
  sub_lin_b = 2.33826005160269e-01, sub_lin_se = 1.07726414518428e-01,
  sub_ddd_b = 2.33826005160249e-01,
  sub_att_m0 = 3.18561836770495e-02, sub_att_m1 = 2.76595636385756e-01,
  sub_datt = 2.44739452708706e-01, sub_datt_se = 1.13517706191325e-01,
  tw_bb = 4.61873389630323e-02, tw_bb_se = 9.93148656735186e-02,
  tw_cb = 2.27617058408507e-02, tw_cb_se = 3.71669983685035e-02,
  tw_cc = 3.30045449204120e-02, tw_cc_se = 2.24837646085153e-02,
  tw_bbp = 1.22983253250856e-02, tw_bbp_se = 9.58342038953850e-02,
  tw_lin_b = 1.57119855530808e-02, tw_lin_se = 3.63434024584511e-02)

st_se <- function(e, V) sqrt(drop(t(e$jac) %*% V %*% e$jac))
st_eq <- function(x, name, tol) expect_equal(unname(x), st[[name]], tolerance = tol, info = name)
## the package's fitting tolerance leaves the estimate a few parts in a
## million from Stata's (it closes with a tighter fit); OLS is exact
TOL_ML <- 2e-5
TOL_OLS <- 1e-9

test_that("the fixtures are the datasets Stata read", {
  lin <- stata_fixture("lin"); did <- stata_fixture("did"); didx <- stata_fixture("didx")
  sub <- stata_fixture("sub"); tw <- stata_fixture("tw")
  st_eq(sum(lin$y), "sum_lin_y", 1e-12); st_eq(sum(lin$z), "sum_lin_z", 1e-12)
  st_eq(sum(did$y), "sum_did_y", 1e-12); st_eq(sum(did$z), "sum_did_z", 1e-12)
  st_eq(sum(didx$y), "sum_didx_y", 1e-12); st_eq(sum(didx$period), "sum_didx_period", 1e-12)
  st_eq(sum(sub$y), "sum_sub_y", 1e-12); st_eq(sum(sub$m), "sum_sub_m", 1e-12)
  st_eq(sum(tw$y_bb), "sum_tw_y_bb", 1e-12); st_eq(sum(tw$y_cc), "sum_tw_y_cc", 1e-12)
  st_eq(sum(tw$c1), "sum_tw_c1", 1e-12); st_eq(sum(tw$x), "sum_tw_x", 1e-12)
})

test_that("linear analysis: coefficient and robust SEs equal Stata's regress", {
  f <- stata_fixture("lin")
  ## regress y d z, vce(robust)
  X <- cbind(1, f$d, f$z)
  fo <- powerape:::fit_ols(X, f$y)
  st_eq(fo$coefficients[2], "lin_hc1_b", TOL_OLS)
  st_eq(sqrt(powerape:::vcov_ols_hc1(fo, X)[2, 2]), "lin_hc1_se", TOL_OLS)
  ## regress y d m dm z, vce(cluster id)
  X <- cbind(1, f$d, f$m, f$d * f$m, f$z)
  fo <- powerape:::fit_ols(X, f$y)
  st_eq(fo$coefficients[4], "lin_cl_b", TOL_OLS)
  st_eq(sqrt(powerape:::vcov_ols_cluster(fo, X, f$id)[4, 4]), "lin_cl_se", TOL_OLS)
})

test_that("diff-in-diff panel: ATT and clustered SE equal logit + margins", {
  f <- stata_fixture("did")
  X <- cbind(1, f$d, as.numeric(f$period == 2), f$w, f$z)
  ## logit y i.d i.period i.w z, vce(cluster id)
  ## margins, dydx(w) subpop(if w == 1)
  ft <- powerape:::fit_index_model(X, f$y, "logit")
  e <- powerape:::att_est(ft$fit$coefficients, X, "logit", 4L)
  st_eq(e$ape, "did_att", TOL_ML)
  st_eq(st_se(e, powerape:::vcov_cluster(ft$fit, X, f$id, "logit")), "did_att_se", TOL_ML)
  ## regress y i.d i.period i.w z, vce(cluster id)
  fo <- powerape:::fit_ols(X, f$y)
  se <- sqrt(powerape:::vcov_ols_cluster(fo, X, f$id)[4, 4])
  st_eq(fo$coefficients[4], "did_lin_b", TOL_OLS)
  st_eq(se, "did_lin_se", TOL_OLS)
  ## xtdidregress (y z) (w), group(id) time(period): the same estimate; its
  ## SE counts the within regressors (period, w, constant), the pooled one
  ## all five columns
  st_eq(fo$coefficients[4], "did_xtdid_b", TOL_OLS)
  N <- nrow(X)
  st_eq(se * sqrt((N - 5) / (N - 3)), "did_xtdid_se", TOL_OLS)
})

test_that("diff-in-diff cross-sections, four periods: ATT, model SE, linear analysis", {
  f <- stata_fixture("didx")
  X <- cbind(1, f$d, as.numeric(f$period == 2), as.numeric(f$period == 3),
             as.numeric(f$period == 4), f$w)
  ## logit y i.d i.period i.w ; margins, dydx(w) subpop(if w == 1)
  ft <- powerape:::fit_index_model(X, f$y, "logit")
  e <- powerape:::att_est(ft$fit$coefficients, X, "logit", 6L)
  st_eq(e$ape, "didx_att", TOL_ML)
  st_eq(st_se(e, powerape:::vcov_from_glmfit(ft$fit)), "didx_att_se", TOL_ML)
  ## regress y i.d i.period i.w, vce(robust) and didregress ..., vce(robust)
  fo <- powerape:::fit_ols(X, f$y)
  se <- sqrt(powerape:::vcov_ols_hc1(fo, X)[6, 6])
  st_eq(fo$coefficients[6], "didx_lin_b", TOL_OLS)
  st_eq(se, "didx_lin_se", TOL_OLS)
  st_eq(fo$coefficients[6], "didx_didreg_b", TOL_OLS)
  st_eq(se, "didx_didreg_se", TOL_OLS)
})

test_that("subgroup: both effects, their difference and its SE equal margins + lincom", {
  f <- stata_fixture("sub")
  post <- as.numeric(f$period == 2)
  X <- cbind(1, f$d, post, f$w, f$m, f$m * f$d, f$m * post, f$m * f$w)
  xx <- list(X = X, focal_col = 8L, w_col = 4L, m = f$m, id = f$id)
  ## logit y i.w##i.m d post md mpost, vce(cluster id)
  ## margins, dydx(w) subpop(if w == 1) over(m) post ; lincom
  ft <- powerape:::fit_index_model(X, f$y, "logit")
  e <- powerape:::datt_est(ft$fit$coefficients, xx, "logit")
  st_eq(e$by[1], "sub_att_m0", TOL_ML)
  st_eq(e$by[2], "sub_att_m1", TOL_ML)
  st_eq(e$ape, "sub_datt", TOL_ML)
  st_eq(st_se(e, powerape:::vcov_cluster(ft$fit, X, f$id, "logit")), "sub_datt_se", TOL_ML)
  ## regress y d post w m md mpost mw, vce(cluster id); the official triple
  ## difference, didregress (y) (mw), group(d m) time(period), is the same
  fo <- powerape:::fit_ols(X, f$y)
  st_eq(fo$coefficients[8], "sub_lin_b", TOL_OLS)
  st_eq(sqrt(powerape:::vcov_ols_cluster(fo, X, f$id)[8, 8]), "sub_lin_se", TOL_OLS)
  st_eq(fo$coefficients[8], "sub_ddd_b", TOL_OLS)
})

test_that("three-way effect and SE equal margins contrasts for every type pair", {
  f <- stata_fixture("tw")
  X3 <- function(a, b) cbind(1, a, b, a * b, f$w, a * f$w, b * f$w, a * b * f$w, f$x)
  run <- function(y, a, b, ta, tb, link = "logit") {
    X <- X3(a, b)
    ft <- powerape:::fit_index_model(X, y, link)
    e <- powerape:::aie3_est(ft$fit$coefficients, X, link, ta, tb)
    list(e = e, fit = ft$fit, X = X, se = st_se(e, powerape:::vcov_from_glmfit(ft$fit)))
  }
  ## logit y_bb i.f##i.m##i.w x ; margins r.f#r.m#r.w
  r <- run(f$y_bb, f$f, f$m, "binary", "binary")
  st_eq(r$e$ape, "tw_bb", TOL_ML); st_eq(r$se, "tw_bb_se", TOL_ML)
  ## logit y_cb c.c1##i.m##i.w x ; margins r.m#r.w, dydx(c1)
  r <- run(f$y_cb, f$c1, f$m, "normal", "binary")
  st_eq(r$e$ape, "tw_cb", TOL_ML); st_eq(r$se, "tw_cb_se", TOL_ML)
  ## logit y_cc c.c1##c.c2##i.w x ; margins r.w, expression(<dPr/dc1>) dydx(c2)
  ## (Stata differentiates the expression numerically)
  r <- run(f$y_cc, f$c1, f$c2, "normal", "normal")
  st_eq(r$e$ape, "tw_cc", 1e-4); st_eq(r$se, "tw_cc_se", 1e-4)
  ## probit y_bbp i.f##i.m##i.w x ; margins r.f#r.m#r.w
  r <- run(f$y_bbp, f$f, f$m, "binary", "binary", "probit")
  ## this estimate is a tenth of its SE, so compare in SE units
  expect_lt(abs(r$e$ape - st[["tw_bbp"]]), 1e-4 * st[["tw_bbp_se"]])
  ## Stata's probit SE uses the observed information: with it the package's
  ## Jacobian reproduces the SE; under the package's expected information
  ## the SE is within one percent
  eta <- drop(r$X %*% r$fit$coefficients); q <- 2 * f$y_bbp - 1
  lam <- q * dnorm(q * eta) / pnorm(q * eta)
  V_oim <- solve(crossprod(r$X * sqrt(lam * (lam + eta))))
  st_eq(st_se(r$e, V_oim), "tw_bbp_se", TOL_ML)
  expect_lt(abs(r$se / st[["tw_bbp_se"]] - 1), 0.01)
  ## linear analysis: regress y_cb c1 m c1m w c1w mw c1mw x, vce(robust)
  X <- X3(f$c1, f$m)
  fo <- powerape:::fit_ols(X, f$y_cb)
  st_eq(fo$coefficients[8], "tw_lin_b", TOL_OLS)
  st_eq(sqrt(powerape:::vcov_ols_hc1(fo, X)[8, 8]), "tw_lin_se", TOL_OLS)
})
