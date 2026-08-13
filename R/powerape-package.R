#' powerape: Power Analysis for Average Partial Effects
#'
#' Simulation-based power analysis in which the effect size of interest is
#' the average partial effect (APE) or the average interaction effect (AIE)
#' from a nonlinear model such as probit or logit, including
#' smallest-effect-size-of-interest (SESOI) powering via confidence-interval
#' criteria for detection, minimum-effect, and equivalence claims
#' (Riesthuis, 2024) and robustness sweeps over contextual assumptions
#' (Hancock and Feng, 2025).
#'
#' @keywords internal
#' @importFrom stats aggregate binomial dlogis dnorm glm.fit isoreg plogis pnorm qlogis qnorm rbinom rnorm sd uniroot var
#' @importFrom graphics abline arrows lines mtext points
#' @importFrom utils packageVersion
"_PACKAGE"
