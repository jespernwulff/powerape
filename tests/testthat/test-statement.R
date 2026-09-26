test_that("power_statement renders APE power, ape_n, and AIE objects", {
  d <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
               covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
               baseline = 0.30, signal = 0.15)
  d <- set_ape(d, 0.10)

  pw <- ape_power(d, n = 600, claim = "minimum", sesoi = 0.03,
                  nsim = 200, seed = 2)
  st <- power_statement(pw)
  expect_s3_class(st, "powerape_statement")
  txt <- unclass(st)
  expect_match(txt, "average partial effect")
  expect_match(txt, "minimum-effect claim")
  expect_match(txt, "0.100")
  expect_match(txt, "n = 600")
  expect_match(txt, "Riesthuis")
  expect_match(txt, "probit")
  expect_output(print(st), "powerape package")

  an <- ape_n(d, power = 0.8, claim = "detect", nsim = 400, seed = 3)
  expect_match(unclass(power_statement(an)), "required total sample size")

  d2 <- ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                moderator = pa_var("female", "binary", p = 0.55),
                baseline = 0.30)
  d2 <- set_aie(d2, 0.08, main_focal = 0.10, main_moderator = 0.05)
  pw2 <- ape_power(d2, n = 1500, claim = "detect", nsim = 150, seed = 4)
  txt2 <- unclass(power_statement(pw2))
  expect_match(txt2, "average interaction effect")
  expect_match(txt2, "moderator")
  expect_match(txt2, "main-effect APEs")
})

test_that("power_statement counts panel n in units, not observations", {
  dp <- ape_dgp_panel("probit",
                      focal = pa_var("treat", "binary", p = 0.5, icc = 1),
                      n_periods = 3, rho = 0.2, baseline = 0.30,
                      n_int = 2e4, seed_int = 7L)
  dp <- set_ape(dp, 0.10)
  pw <- ape_power(dp, n = 80, claim = "detect", nsim = 60, seed = 5)
  txt <- unclass(power_statement(pw))
  expect_match(txt, "n = 80 units observed over 3 periods")
  expect_match(txt, "240 unit-period observations")
  expect_no_match(txt, "total sample size")
})

test_that("power_statement states the dependence, marginals, and persistence it priced", {
  mk <- function(r) set_ape(ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                                    covariates = list(pa_var("age", "normal", mean = 45, sd = 12)),
                                    correlation = r, baseline = 0.30, signal = 0.15), 0.10)
  st <- function(r) unclass(power_statement(
    ape_power(mk(r), n = 600, claim = "detect", nsim = 30, seed = 1)))
  s0 <- st(NULL); s2 <- st(0.2); s7 <- st(0.7)
  expect_match(s0, "mutually independent")
  expect_no_match(s0, "Gaussian-copula")          # no dependence was specified
  expect_match(s2, "latent (Gaussian-copula) correlation of 0.20", fixed = TRUE)
  expect_match(s7, "correlation of 0.70", fixed = TRUE)
  expect_match(s0, "age (continuous, mean 45, SD 12)", fixed = TRUE)
  ## worlds that price differently must not read the same
  mask <- function(s) gsub("[0-9.]+", "#", s)
  expect_false(identical(mask(s0), mask(s2)))

  mkp <- function(icc) set_ape(ape_dgp_panel(
    focal = pa_var("treat", "binary", p = 0.5, icc = icc),
    covariates = list(pa_var("size", "normal", icc = 0.6)),
    n_periods = 3, rho = 0.2, baseline = 0.30, signal = 0.1, n_int = 2e4), 0.10)
  sp <- function(icc) unclass(power_statement(
    ape_power(mkp(icc), n = 60, claim = "detect", nsim = 30, seed = 2)))
  expect_match(sp(1), "time-constant")
  expect_match(sp(0), "redrawn each period")
  expect_match(sp(1), "within-unit persistence (icc) 0.60", fixed = TRUE)

  ## ape_n keeps its standard-error type, and the statement says so
  an <- ape_n(mk(NULL), power = 0.8, claim = "detect", nsim = 200, seed = 3,
              se = "robust", confirm = FALSE)
  expect_identical(an$se, "robust")
  expect_match(unclass(power_statement(an)), "sandwich")
  expect_output(print(an), "robust SEs")

  ## empirical route: the resampling mode is stated
  set.seed(4)
  pilot <- data.frame(treat = rbinom(600, 1, 0.5), age = rnorm(600, 45, 12))
  de <- set_ape(ape_dgp_empirical(data = pilot, focal = "treat",
                                  baseline = 0.30, signal = 0.1), 0.10)
  expect_match(unclass(power_statement(
    ape_power(de, n = 400, claim = "detect", nsim = 30, seed = 5))),
    "resampled jointly")
})

test_that("power_statement rejects other objects", {
  expect_error(power_statement(42))
})
