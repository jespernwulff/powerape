test_that("ape_n confirmation tightens the answer and reports itself", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  an <- ape_n(d, power = 0.80, claim = "detect", nsim = 1000, seed = 11)
  expect_true(an$confirmed)
  expect_true(any(an$history$stage == "confirm"))
  expect_identical(an$nsim_confirm, 4000L)
  ## the returned power/MCSE come from the high-precision confirmation run
  expect_lt(an$mcse, 0.009)
  ## exact finite-sample power at the returned n is close to the goal
  ex <- exact_power_sat(an$n, 0.5, 0.30, 0.40, 0.95, "detect")
  expect_gt(ex, 0.80 - 0.015)
  expect_output(print(an), "confirmed")
})

test_that("confirm = FALSE restores the fast single-stage behaviour", {
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  an <- ape_n(d, power = 0.80, claim = "detect", nsim = 800, seed = 12,
              confirm = FALSE)
  expect_false(an$confirmed)
  expect_false(any(an$history$stage == "confirm"))
  expect_true(is.na(an$nsim_confirm))
  expect_output(print(an), "larger nsim")
})

test_that("the confirmation rule accepts within max(.005, 1.5 MCSE) and pushes n up beyond it", {
  ## power_once is mocked with a scripted sequence, so the rule itself is
  ## tested: the seeded run above accepts in its first round, where the
  ## tolerance never binds (validation-gaps review, 2026-10-14)
  d <- set_ape(ape_dgp("probit", focal = pa_var("t", "binary", p = 0.5),
                       baseline = 0.30), 0.10)
  scripted <- function(seq_pw) {
    i <- 0L
    ## the real formals: a `...` mock would let the caller's `se = se`
    ## partially match `sesoi`
    function(dgp, n, claim, sesoi, conf = NULL, nsim, seed,
             enforce = TRUE, se = "model", alpha = NULL,
             keep_draws = FALSE, separation = "fail") {
      i <<- i + 1L
      p <- seq_pw[min(i, length(seq_pw))]
      list(power = p, mcse = sqrt(p * (1 - p) / nsim), alpha = 0.05,
           separated = 0, min_cell = NA_real_)
    }
  }
  ## search hits .80; the first confirmation (nsim 4,000) is .02 short, beyond
  ## the .0098 tolerance; the pushed candidate confirms
  local({
    testthat::local_mocked_bindings(power_once = scripted(c(0.80, 0.78, 0.80)),
                                    .package = "powerape")
    an <- ape_n(d, power = 0.80, claim = "detect", nsim = 1000, seed = 1)
    h <- an$history[an$history$stage == "confirm", ]
    expect_identical(nrow(h), 2L)
    expect_gt(h$n[2L], h$n[1L])
    expect_true(an$confirmed)
  })
  ## .0095 short sits inside 1.5 MCSE (.0097 at .7905): accepted at once
  local({
    testthat::local_mocked_bindings(power_once = scripted(c(0.80, 0.7905)),
                                    .package = "powerape")
    an <- ape_n(d, power = 0.80, claim = "detect", nsim = 1000, seed = 1)
    expect_identical(sum(an$history$stage == "confirm"), 1L)
    expect_true(an$confirmed)
  })
})
