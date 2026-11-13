#' powerape: Power Analysis for Average Partial Effects
#'
#' Simulation-based power analysis in which the effect size of interest is
#' the average partial effect (APE) or the average interaction effect (AIE)
#' from a nonlinear model such as probit or logit, including
#' smallest-effect-size-of-interest (SESOI) powering via confidence-interval
#' criteria for detection, minimum-effect, and equivalence claims
#' (Riesthuis, 2024) and robustness sweeps over contextual assumptions
#' (Hancock and Feng, 2026).
#'
#' Designs are built by [ape_dgp()] (cross-sections, with an optional
#' moderator or two), [ape_dgp_empirical()] and [ape_dgp_from_fit()] (pilot
#' data), [ape_dgp_panel()] (panels), [ape_dgp_iv()] (endogenous focal
#' variables), and [ape_dgp_did()] and [ape_dgp_did_pilot()]
#' (difference-in-differences, with the effect on the treated in
#' probability points as the target). The effect is pinned with
#' [set_ape()], [set_aie()], [set_aie3()] or [set_att()], and
#' [set_analysis()] switches a design to a linear probability analysis for
#' comparison. [ape_power()], [ape_curve()], [ape_n()], [ape_mde()] and
#' [ape_robust()] compute power, required sample sizes, detectable effects
#' and sensitivity; [ape_alpha()] justifies the error rate, and
#' [power_statement()] writes the methods paragraph.
#'
#' @keywords internal
#' @importFrom stats aggregate binomial dlogis dnorm glm.fit integrate isoreg lm.fit plogis pnorm poisson qlogis qnorm rbinom rnorm rpois runif sd uniroot var
#' @importFrom graphics abline arrows lines mtext points
#' @importFrom utils packageVersion
"_PACKAGE"
