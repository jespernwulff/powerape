# The linear probability analysis ----------------------------------------------
# Every design is analysed by its index
# model (probit/logit) unless its `analysis` field says "linear": then the
# simulated studies are fitted by ordinary least squares on the same design
# matrix, the estimate is the coefficient of the focal variable (or of the
# focal-by-moderator product for AIE designs, or of the treatment indicator
# for difference-in-differences designs), and the standard errors are
# always robust: HC1 for independent observations, clustered by unit when
# units are observed repeatedly.

is_linear <- function(dgp) identical(dgp$analysis, "linear")

# Standard errors of the linear analysis, as named in prints.
linear_se_label <- function(dgp) {
  if (unit_design(dgp)) ", unit-clustered SEs" else ", robust (HC1) SEs"
}

# Design-matrix column whose coefficient is the linear analysis' estimate.
# `xx` is a draw from draw_x(); AIE designs keep the interaction in column 4.
linear_col <- function(dgp, xx) {
  if (is_did(dgp)) xx$focal_col
  else if (is_aie3(dgp)) 8L          # the triple product (aie3.R)
  else if (identical(dgp$estimand, "aie") || !is.null(dgp$moderator)) 4L
  else xx$focal_col
}

#' Choose the analysis the simulated studies receive
#'
#' By default every simulated study is analysed with the design's index
#' model: a probit or logit fit, with the average partial effect (or
#' average interaction effect) and its delta-method standard error.
#' `set_analysis(dgp, "linear")` returns a copy of the design that is
#' analysed with a **linear probability model** instead: ordinary least
#' squares on the same regressors, where the estimate is the coefficient
#' of the focal variable (APE designs), of the focal-by-moderator product
#' (AIE designs), of the triple product (three-way designs), or of the
#' treatment indicator (difference-in-differences designs). One design can
#' so be priced both ways, for comparison.
#'
#' The world and the pinned effect do not change: the truth remains the
#' average partial effect (or AIE, or ATT) set with [set_ape()],
#' [set_aie()], or [set_att()]. The linear coefficient is another
#' estimator of that effect and need not centre on it. In large samples
#' it equals the average partial effect when a binary focal variable is
#' independent of the covariates (as under randomization) and when the
#' focal variable and the covariates are jointly normal; it equals the
#' average interaction effect when a binary focal variable and a binary
#' moderator are independent of the covariates. Otherwise it differs, and
#' it can differ substantially when a binary focal variable is related to
#' covariates that also move the outcome: the linear control then leaves
#' part of the nonlinear confounding in the coefficient. (With a 30%
#' treated group, a 15% baseline, one covariate with `signal = 0.3` and a
#' latent correlation of 0.3 between the two, the linear coefficient
#' centres about 12% above the average partial effect.) The printed design
#' shows the value the linear analysis centres on in large samples next to
#' the truth, and [ape_power()] reports the mean of the simulated
#' estimates.
#'
#' @section Standard errors and separation:
#' A linear probability model is heteroskedastic by construction, so its
#' standard errors are always robust and the `se` argument of the power
#' functions is ignored: HC1 for independent observations (as Stata's
#' `regress, vce(robust)`) and clustered by unit for panel designs (with
#' the small-sample factor of Stata's `regress, vce(cluster)`). Ordinary
#' least squares has no separation, so no simulated study is discarded for
#' a cell without events or without non-events; the sparse-cell warning
#' still applies.
#'
#' @section Limits:
#' Instrumental-variable designs are refused: their linear counterpart is
#' two-stage least squares, which is not available. [power_statement()]
#' describes a linear analysis for difference-in-differences designs
#' only. Before searching, [ape_n()]
#' checks that the value the linear analysis centres on lies on the same
#' side of the claim's boundary as the true effect, and stops otherwise.
#'
#' @param dgp A `powerape_dgp`.
#' @param analysis `"index"` (the design's probit or logit model, the
#'   default of every design) or `"linear"` (linear probability model).
#'
#' @return A copy of `dgp` with the analysis set.
#' @examples
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
#'              covariates = list(pa_var("age", "normal")),
#'              baseline = 0.30, signal = 0.15)
#' d <- set_ape(d, target = 0.10)
#' dl <- set_analysis(d, "linear")
#' dl
#' \donttest{
#' ape_power(d, n = 800, claim = "detect", nsim = 300, seed = 1)
#' ape_power(dl, n = 800, claim = "detect", nsim = 300, seed = 1)
#' }
#' @export
set_analysis <- function(dgp, analysis = c("index", "linear")) {
  stopifnot(inherits(dgp, "powerape_dgp"))
  analysis <- match.arg(analysis)
  if (identical(dgp$route, "iv") && analysis == "linear")
    stop(paste("A linear analysis is not available for instrumental-variable",
               "designs: its counterpart there is two-stage least squares,",
               "a different estimator from the control-function probit this",
               "route simulates."), call. = FALSE)
  if (is_did(dgp)) {
    dgp$analysis <- analysis
    dgp$spec$analysis <- analysis
  } else {
    ## standard designs carry no field under the default analysis
    dgp$analysis <- if (analysis == "linear") "linear" else NULL
  }
  dgp
}

# Large-sample value of the linear analysis' estimate: OLS of the true
# outcome probabilities on the design matrix in one large draw with a fixed
# seed (deterministic; the caller's RNG state is preserved). Using the
# probabilities instead of simulated outcomes removes the Bernoulli noise.
linear_centre <- function(dgp, n_big = NULL) {
  if (is.null(dgp$target_est))
    stop("This DGP has no target effect yet.", call. = FALSE)
  if (identical(dgp$route, "iv"))
    stop("No linear analysis for instrumental-variable designs.", call. = FALSE)
  n_big <- n_big %||% (if (unit_design(dgp)) 1e5 else 4e5)
  with_seed(dgp$seed_int %||% 20260812L, {
    xx <- draw_x(dgp, n_big)
    pr <- if (identical(dgp$route, "panel")) xx$pr
          else dgp$G(drop(xx$X %*% beta_true(dgp)))
    fo <- fit_ols(xx$X, pr)
    if (!isTRUE(fo$ok)) NA_real_ else fo$coefficients[[linear_col(dgp, xx)]]
  })
}

# Before a search: the linear analysis must centre on the same side of the
# claim's boundary as the truth, or "the n that reaches the goal" does not
# exist (or answers another question).
linear_need_side <- function(dgp, claim, sesoi, what) {
  if (!is_linear(dgp) || is_did(dgp)) return(invisible(TRUE))
  cen <- linear_centre(dgp)
  tgt <- dgp$target_est
  sgn <- if (tgt < 0) -1 else 1
  ok <- switch(claim,
               detect = sgn * cen > 0,
               minimum = sgn * cen > sesoi,
               equivalence = abs(cen) < sesoi)
  if (!isTRUE(ok))
    stop(sprintf(paste(
      "%s: the linear analysis centres on %+.4f in large samples, the true",
      "%s is %+.4f, and the two are on different sides of the %s claim's",
      "boundary, so no sample size delivers this claim for the right",
      "reason. Compare the two analyses with ape_power() instead."),
      what, cen, est_label(dgp$estimand %||% "ape"), tgt, claim), call. = FALSE)
  invisible(TRUE)
}
