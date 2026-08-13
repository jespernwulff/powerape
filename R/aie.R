# Average interaction effects: pinning and truth evaluation --------------------
# Definitions follow the ginteff package (Wulff; R port of Radean's Stata
# command): first differences for binary variables, derivatives for
# continuous ones, averaged over the covariate distribution.

# True AIE on the integration set, in the centered parameterization.
truth_aie <- function(dgp, bf, bm, bi, id) {
  G <- dgp$G
  g <- dgp$g
  gp <- gprime_fun(dgp$link)
  b0 <- dgp$beta0
  xz <- id$idxz
  td <- dgp$focal$type
  tm <- dgp$moderator$type
  xd <- id$xd_c
  md <- id$md_c
  if (td == "binary" && tm == "binary") {
    mean(G(b0 + bf + bm + bi + xz)) - mean(G(b0 + bm + xz)) -
      mean(G(b0 + bf + xz)) + mean(G(b0 + xz))
  } else if (td != "binary" && tm == "binary") {
    (bf + bi) * mean(g(b0 + bf * xd + bm + bi * xd + xz)) -
      bf * mean(g(b0 + bf * xd + xz))
  } else if (td == "binary" && tm != "binary") {
    (bm + bi) * mean(g(b0 + bf + bm * md + bi * md + xz)) -
      bm * mean(g(b0 + bm * md + xz))
  } else {
    idx <- b0 + bf * xd + bm * md + bi * xd * md + xz
    mean(gp(idx) * (bf + bi * md) * (bm + bi * xd) + g(idx) * bi)
  }
}

#' Pin the true AIE of a moderated DGP (inversion step)
#'
#' For a DGP with a moderator, solves all three effect coefficients from
#' targets in effect-size units: the two main effects from their
#' *conditional-at-reference APEs* (the APE of each variable with the other
#' held at its reference value -- 0 for binary, the mean for continuous),
#' and the interaction coefficient from the target average interaction
#' effect (AIE), holding the mains fixed. The AIE is `ginteff`'s estimand:
#' the double difference for binary-by-binary, the difference of average
#' derivatives for mixed pairs, and the average cross-partial for
#' continuous-by-continuous.
#'
#' All three anchors are explicit -- there are no silent defaults. Errors
#' report the attainable range when a target is infeasible (for
#' binary-by-binary, the four counterfactual cell rates must stay in (0,1);
#' otherwise the AIE plateaus in the interaction coefficient).
#'
#' @param dgp A `powerape_dgp` from [ape_dgp()] with a `moderator`.
#' @param target Assumed true AIE (the planning value), in APE units.
#' @param main_focal,main_moderator Conditional-at-reference APEs anchoring
#'   the two main effects.
#'
#' @return The `dgp` with `beta_focal`, `beta_mod`, `beta_int`, and
#'   `target_est` filled in; `estimand` becomes `"aie"`.
#' @examples
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5),
#'              moderator = pa_var("female", "binary", p = 0.55),
#'              baseline = 0.30)
#' d <- set_aie(d, target = 0.04, main_focal = 0.10, main_moderator = 0.05)
#' d
#' @export
set_aie <- function(dgp, target, main_focal, main_moderator) {
  stopifnot(inherits(dgp, "powerape_dgp"))
  if (is.null(dgp$moderator))
    stop("This DGP has no moderator: add one in ape_dgp(), or pin an APE with set_ape().",
         call. = FALSE)
  if (missing(main_focal) || missing(main_moderator))
    stop(paste("set_aie() needs explicit main-effect anchors: `main_focal` and",
               "`main_moderator` are the conditional-at-reference APEs (the APE",
               "of each variable with the other held at its reference value)."),
         call. = FALSE)
  stopifnot(is.numeric(target), length(target) == 1L, is.finite(target),
            is.numeric(main_focal), length(main_focal) == 1L,
            is.numeric(main_moderator), length(main_moderator) == 1L)

  G <- dgp$G
  id <- integration_draw(dgp)
  idxz <- id$idxz
  p0 <- mean(G(dgp$beta0 + idxz))

  ## 1. Main effects, each at the other variable's reference (interaction
  ##    term drops out, so these are ordinary single-variable inversions).
  bf <- if (dgp$focal$type == "binary") {
    solve_binary_shift(dgp, main_focal, idxz, p0, "main-effect APE (focal)")
  } else {
    solve_cont_main(dgp, main_focal, id$xd_c, idxz, pa_var_sd(dgp$focal),
                    "main-effect APE (focal)")
  }
  bm <- if (dgp$moderator$type == "binary") {
    solve_binary_shift(dgp, main_moderator, idxz, p0, "main-effect APE (moderator)")
  } else {
    solve_cont_main(dgp, main_moderator, id$md_c, idxz, pa_var_sd(dgp$moderator),
                    "main-effect APE (moderator)")
  }

  ## 2. Interaction coefficient from the AIE target, mains held fixed.
  s_d <- if (dgp$focal$type == "binary") 1 else pa_var_sd(dgp$focal)
  s_m <- if (dgp$moderator$type == "binary") 1 else pa_var_sd(dgp$moderator)
  bi <- solve_interaction(function(b) truth_aie(dgp, bf, bm, b, id),
                          target, s_d * s_m, "AIE")

  dgp$beta_focal <- bf
  dgp$beta_mod <- bm
  dgp$beta_int <- bi
  dgp$main_focal <- main_focal
  dgp$main_moderator <- main_moderator
  dgp$estimand <- "aie"
  dgp$target_est <- target
  dgp$p1 <- NULL
  dgp
}

#' Recompute the true AIE implied by a pinned DGP
#'
#' @param dgp A `powerape_dgp` after [set_aie()].
#' @return The true AIE (numeric scalar), evaluated on the deterministic
#'   integration set; equals the [set_aie()] target up to root-finding
#'   tolerance.
#' @export
true_aie <- function(dgp) {
  if (!identical(dgp$estimand, "aie"))
    stop("This DGP has no pinned AIE - call set_aie() first.", call. = FALSE)
  truth_aie(dgp, dgp$beta_focal, dgp$beta_mod, dgp$beta_int, integration_draw(dgp))
}
