# Inversion helpers: from targets in effect-size units to coefficients --------

# Binary shift coefficient: solve mean(G(beta0 + b + idxz)) - p_ref = target.
solve_binary_shift <- function(dgp, target, idxz, p_ref, label) {
  eps <- 1e-4
  if (target <= -p_ref + eps || target >= (1 - p_ref) - eps)
    stop(sprintf(paste(
      "Target %s %.4f is infeasible for this DGP: with baseline",
      "P(Y=1 | reference) = %.3f the attainable range is (%.3f, %.3f).",
      "Adjust the target or the baseline."),
      label, target, p_ref, -p_ref, 1 - p_ref), call. = FALSE)
  G <- dgp$G
  uniroot(function(b) mean(G(dgp$beta0 + b + idxz)) - p_ref - target,
          c(-8, 8), extendInt = "upX", tol = 1e-9)$root
}

# Continuous main effect: delta(b) = b * mean(g(beta0 + b*x_c + idxz)) has
# delta(0) = 0 and plateaus as |b| grows; probe a grid on the target's sign
# side, error if the target exceeds the plateau, else bracket and root-find.
solve_cont_main <- function(dgp, target, x_c, idxz, scale, label) {
  g <- dgp$g
  f <- function(b) b * mean(g(dgp$beta0 + b * x_c + idxz))
  s <- if (target < 0) -1 else 1
  bgrid <- s * exp(seq(log(1e-3 / scale), log(200 / scale), length.out = 100))
  vals <- vapply(bgrid, f, numeric(1))
  amax <- max(abs(vals))
  if (abs(target) >= amax * (1 - 1e-3))
    stop(sprintf(paste(
      "Target %s %.4f is infeasible for this DGP: the attainable range is",
      "about (%.4f, %.4f). The average derivative plateaus as the",
      "coefficient grows because the index disperses; adjust the target,",
      "the baseline, or the focal variable's spread."),
      label, target, -amax, amax), call. = FALSE)
  i <- which(abs(vals) >= abs(target))[1L]
  lo <- if (i == 1L) 0 else bgrid[i - 1L]
  uniroot(function(b) f(b) - target, sort(c(lo, bgrid[i])), tol = 1e-9)$root
}

# Interaction coefficient: f(bi) need not pass through the origin (the
# mechanical Ai-Norton interaction is nonzero at bi = 0), so scan a signed
# grid for the crossing nearest bi = 0; error with the attainable range if
# the target is never reached.
solve_interaction <- function(f, target, scale, label) {
  gr <- exp(seq(log(1e-3 / scale), log(50 / scale), length.out = 60))
  bs <- c(-rev(gr), 0, gr)
  vals <- vapply(bs, f, numeric(1))
  h <- vals - target
  cross <- which(h[-1L] * h[-length(h)] <= 0)
  if (!length(cross))
    stop(sprintf(paste(
      "Target %s %.4f is infeasible for this DGP: the attainable range is",
      "about (%.4f, %.4f). Adjust the target, the main effects, or the",
      "baseline."),
      label, target, min(vals), max(vals)), call. = FALSE)
  mid <- abs((bs[cross] + bs[cross + 1L]) / 2)
  j <- cross[which.min(mid)]
  uniroot(function(b) f(b) - target, c(bs[j], bs[j + 1L]), tol = 1e-9)$root
}

#' Pin the true APE of a DGP (inversion step)
#'
#' Solves for the focal coefficient such that the data-generating process in
#' `dgp` has average partial effect exactly `target` (in probability
#' points): the average discrete change for a binary focal variable, the
#' average derivative for a continuous one. Root-finding runs on the
#' deterministic integration set. Errors if the target is infeasible --
#' for a binary focal the APE must lie inside `(-p0, 1 - p0)` given baseline
#' rate p0; for a continuous focal the APE plateaus in the coefficient, and
#' the attainable range is reported. For DGPs with a moderator, use
#' [set_aie()] instead.
#'
#' @param dgp A `powerape_dgp` from [ape_dgp()] or [ape_dgp_empirical()],
#'   without a moderator.
#' @param target Assumed true APE (the planning value), e.g. `0.05` for five
#'   probability points. Always explicit -- there is no default truth.
#'
#' @return The `dgp` with the focal coefficient and `target_est` filled in
#'   (and, for a binary focal, the implied `P(Y = 1 | focal = 1)`).
#' @examples
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' d
#' @export
set_ape <- function(dgp, target) {
  stopifnot(inherits(dgp, "powerape_dgp"),
            is.numeric(target), length(target) == 1L, is.finite(target))
  if (!is.null(dgp$moderator))
    stop("This DGP has a moderator: pin the estimand with set_aie().", call. = FALSE)
  if (identical(dgp$route, "panel")) {
    id <- panel_integration(dgp)
    root <- if (dgp$focal$type == "binary") {
      solve_binary_panel(dgp, target, id$q)
    } else {
      solve_cont_panel(dgp, target, id$xd_c, id$q)
    }
    dgp$p1 <- if (dgp$focal$type == "binary") {
      mean(dgp$mix_P(dgp$beta0 + root + id$q))
    } else NULL
    dgp$beta_focal <- root
    dgp$estimand <- "ape"
    dgp$target_est <- target
    return(dgp)
  }
  G <- dgp$G
  if (dgp$focal$type == "binary") {
    if (dgp$route == "parametric" && dgp$k == 0L) {
      idxz <- 0
      p0 <- G(dgp$beta0)
    } else {
      idxz <- integration_draw(dgp)$idxz
      p0 <- mean(G(dgp$beta0 + idxz))
    }
    root <- solve_binary_shift(dgp, target, idxz, p0, "APE")
    dgp$p1 <- mean(G(dgp$beta0 + root + idxz))
  } else {
    id <- integration_draw(dgp)
    root <- solve_cont_main(dgp, target, id$xd_c, id$idxz,
                            pa_var_sd(dgp$focal), "APE")
    dgp$p1 <- NULL
  }
  dgp$beta_focal <- root
  dgp$estimand <- "ape"
  dgp$target_est <- target
  dgp
}

#' Recompute the true APE implied by a pinned DGP
#'
#' @param dgp A `powerape_dgp` after [set_ape()].
#' @return The true APE (numeric scalar), evaluated on the deterministic
#'   integration set; equals the [set_ape()] target up to root-finding
#'   tolerance.
#' @export
true_ape <- function(dgp) {
  if (!identical(dgp$estimand, "ape"))
    stop("This DGP has no pinned APE - call set_ape() first.", call. = FALSE)
  if (identical(dgp$route, "panel")) return(true_ape_panel(dgp))
  G <- dgp$G
  if (dgp$focal$type == "binary") {
    if (dgp$route == "parametric" && dgp$k == 0L)
      return(G(dgp$beta0 + dgp$beta_focal) - G(dgp$beta0))
    idxz <- integration_draw(dgp)$idxz
    mean(G(dgp$beta0 + dgp$beta_focal + idxz) - G(dgp$beta0 + idxz))
  } else {
    id <- integration_draw(dgp)
    dgp$beta_focal * mean(dgp$g(dgp$beta0 + dgp$beta_focal * id$xd_c + id$idxz))
  }
}
