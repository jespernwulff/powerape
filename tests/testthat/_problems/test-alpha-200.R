# Extracted from test-alpha.R:200

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "powerape", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
make_world <- function(target = 0.10) {
  set_ape(ape_dgp("probit", focal = pa_var("treat", "binary", p = 0.5),
                  covariates = list(pa_var("z", "normal")),
                  baseline = 0.30, signal = 0.10), target)
}
synthetic_draws <- function(pw, se0, nsim = 4000L) {
  est <- pw$target + se0 * qnorm(ppoints(nsim))
  pw$draws <- list(est = est, se = rep(se0, nsim), ok = rep(TRUE, nsim))
  pw$nsim <- nsim
  pw
}

# test -------------------------------------------------------------------------
skip_if_not_installed("JustifyAlpha")
d <- make_world()
N <- 20000L
se0 <- 0.0201
pw <- ape_power(d, 2000, claim = "minimum", sesoi = 0.05, nsim = 40, seed = 20)
pw$draws <- list(est = pw$target + se0 * qnorm(ppoints(N)), se = rep(se0, N),
                   ok = rep(TRUE, N))
pw$nsim <- N
assign(".powerape_pf_smooth", function(alpha) 1 - pnorm(qnorm(1 - alpha) - 0.05 / se0),
         envir = globalenv())
on.exit(rm(".powerape_pf_smooth", envir = globalenv()), add = TRUE)
cases <- expand.grid(error = c("minimize", "balance"), cost = c(1, 4), prior = c(1, 3),
                       stringsAsFactors = FALSE)
d_alpha <- d_w <- d_beta <- numeric(nrow(cases))
for (i in seq_len(nrow(cases))) {
    aa <- suppressWarnings(ape_alpha(pw, cost = cases$cost[i], prior = cases$prior[i],
                                     error = cases$error[i], cap = 0.4999))
    ja <- JustifyAlpha::optimal_alpha(".powerape_pf_smooth(alpha = x)",
                                      costT1T2 = cases$cost[i], priorH1H0 = cases$prior[i],
                                      error = cases$error[i], verbose = FALSE)
    d_alpha[i] <- abs(aa$alpha - ja$alpha)
    d_w[i] <- abs(aa$wcer - ja$errorrate)
    d_beta[i] <- abs(aa$beta - ja$beta)
  }
expect_true(all(d_alpha <= 5e-4 + 1e-3 + 1e-6))
expect_true(all(d_w[cases$error == "minimize"] < 1e-4))
expect_true(all(d_beta < 2e-3))
