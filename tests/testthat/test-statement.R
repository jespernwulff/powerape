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

test_that("power_statement rejects other objects", {
  expect_error(power_statement(42))
})
