# Robust power / n_max mode (Hancock & Feng, 2025) -----------------------------

# Rebuild a DGP from its stored constructor spec with some inputs replaced.
rebuild_dgp <- function(dgp, spec) {
  suppressWarnings(
    switch(dgp$builder %||% "parametric",
           from_fit = do.call(ape_dgp_from_fit, spec),
           empirical = do.call(ape_dgp_empirical, spec),
           panel = do.call(ape_dgp_panel, spec),
           do.call(ape_dgp, spec))
  )
}

#' Robustness sweep over contextual assumptions (worst-case power and n_max)
#'
#' Power computed from a single guess at the contextual inputs -- baseline
#' rate, nuisance signal, covariate dependence -- is fragile (Hancock &
#' Feng, 2025). `ape_robust()` re-runs the power analysis over a scenario
#' grid of those inputs and reports per-scenario power, the worst case, and
#' (optionally) `n_max`: the sample size that reaches the target power in
#' the worst scenario -- the insurance-premium n.
#'
#' Two pinning modes (DESIGN.md section 3.3):
#' * `pin = "ape"` (default): the inversion is re-run in every scenario, so
#'   the true effect stays at the planning value and only the *precision*
#'   channel varies. Scenarios where the target is infeasible are reported
#'   as such, not dropped silently.
#' * `pin = "coefficients"`: the effect coefficients stay at their
#'   base-case values while the contextual calibration is redone; the
#'   implied true effect drifts and is reported per scenario. Claim-boundary
#'   guards are relaxed here on purpose -- showing that a scenario collapses
#'   power is the point.
#'
#' @inheritParams ape_power
#' @param vary Named list of scenario values for contextual inputs.
#'   Parametric route: `baseline`, `signal`, `correlation` (scalar);
#'   empirical route: `baseline`, `signal`. A length-2 element is expanded
#'   to `grid_points` equally spaced values; longer vectors are used as-is.
#'   Scenarios are the full factorial grid.
#' @param pin `"ape"` (re-invert per scenario) or `"coefficients"` (hold
#'   coefficients, let the implied effect drift).
#' @param grid_points Grid size used to expand length-2 `vary` elements.
#' @param nmax Run the n_max search in the worst scenario (default TRUE).
#' @param nmax_power Target power for the n_max search (default 0.90).
#'
#' @return A `powerape_robust` object: `scenarios` (inputs, implied effect,
#'   power, MCSE), the worst scenario, marginal mean power per input, and
#'   the n_max result.
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' ape_robust(d, n = 700, claim = "detect",
#'            vary = list(baseline = c(0.20, 0.40)),
#'            nsim = 300, seed = 1, nmax = FALSE)
#' }
#' @export
ape_robust <- function(dgp, n, claim = c("minimum", "detect", "equivalence"),
                       sesoi = NULL, conf = 0.95, nsim = 800, seed = NULL,
                       vary, pin = c("ape", "coefficients"),
                       grid_points = 3, nmax = TRUE, nmax_power = 0.90) {
  claim <- match.arg(claim)
  pin <- match.arg(pin)
  check_coherence(dgp, claim, sesoi, conf)
  stopifnot(is.list(vary), length(vary) >= 1L)
  if (is.null(names(vary)) || !all(nzchar(names(vary))))
    stop("`vary` must be a fully named list.", call. = FALSE)
  allowed <- switch(dgp$builder %||% "parametric",
                    from_fit = "baseline",
                    empirical = c("baseline", "signal"),
                    panel = c("baseline", "signal", "rho"),
                    c("baseline", "signal", "correlation"))
  bad <- setdiff(names(vary), allowed)
  if (length(bad))
    stop(sprintf("Cannot vary %s for this DGP route; varyable inputs: %s.",
                 paste0("`", bad, "`", collapse = ", "),
                 paste(allowed, collapse = ", ")), call. = FALSE)
  if ("signal" %in% names(vary) && !is.null(dgp$spec$gamma))
    stop("Cannot vary `signal` when the DGP was built with explicit `gamma`.",
         call. = FALSE)
  if ("correlation" %in% names(vary) && is.matrix(dgp$spec$correlation))
    stop("Varying `correlation` needs a scalar correlation in the base DGP.",
         call. = FALSE)

  levels_list <- lapply(vary, function(v) {
    stopifnot(is.numeric(v), length(v) >= 2L)
    if (length(v) == 2L) seq(v[1L], v[2L], length.out = grid_points) else v
  })
  grid <- expand.grid(levels_list, KEEP.OUT.ATTRS = FALSE)

  is_aie <- identical(dgp$estimand, "aie")
  rows <- vector("list", nrow(grid))
  dgps <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    spec <- dgp$spec
    for (nm in names(grid)) spec[[nm]] <- grid[i, nm]
    d_i <- tryCatch(rebuild_dgp(dgp, spec), error = function(e) e)
    if (!inherits(d_i, "error")) {
      if (pin == "ape") {
        d_i <- tryCatch(
          if (is_aie) {
            set_aie(d_i, dgp$target_est, main_focal = dgp$main_focal,
                    main_moderator = dgp$main_moderator)
          } else {
            set_ape(d_i, dgp$target_est)
          },
          error = function(e) e)
      } else {
        d_i$beta_focal <- dgp$beta_focal
        d_i$beta_mod <- dgp$beta_mod
        d_i$beta_int <- dgp$beta_int
        d_i$estimand <- dgp$estimand
        d_i$target_est <- dgp$target_est
        d_i$main_focal <- dgp$main_focal
        d_i$main_moderator <- dgp$main_moderator
        implied <- if (is_aie) true_aie(d_i) else true_ape(d_i)
        d_i$target_est <- implied
      }
    }
    if (inherits(d_i, "error")) {
      rows[[i]] <- data.frame(grid[i, , drop = FALSE],
                              implied_effect = NA_real_,
                              power = NA_real_, mcse = NA_real_,
                              note = conditionMessage(d_i),
                              stringsAsFactors = FALSE)
      next
    }
    si <- if (is.null(seed)) NULL else seed + i
    pw <- power_once(d_i, n, claim, sesoi, conf, nsim, si, enforce = FALSE)
    rows[[i]] <- data.frame(grid[i, , drop = FALSE],
                            implied_effect = d_i$target_est,
                            power = pw$power, mcse = pw$mcse, note = "",
                            stringsAsFactors = FALSE)
    dgps[[i]] <- d_i
  }
  scenarios <- do.call(rbind, rows)
  rownames(scenarios) <- NULL

  ok_idx <- which(!is.na(scenarios$power))
  if (!length(ok_idx))
    stop("No scenario was feasible; widen or shift the `vary` ranges.", call. = FALSE)
  worst_i <- ok_idx[which.min(scenarios$power[ok_idx])]

  marginals <- lapply(names(vary), function(nm) {
    ag <- aggregate(scenarios$power[ok_idx],
                    by = list(value = scenarios[[nm]][ok_idx]), FUN = mean)
    names(ag) <- c("value", "mean_power")
    ag
  })
  names(marginals) <- names(vary)

  nmax_res <- NULL
  if (nmax) {
    nmax_res <- tryCatch(
      ape_n(dgps[[worst_i]], power = nmax_power, claim = claim, sesoi = sesoi,
            conf = conf, nsim = nsim,
            seed = if (is.null(seed)) NULL else seed + nrow(grid) + 1L),
      error = function(e) e)
    if (inherits(nmax_res, "error")) {
      warning("n_max search failed in the worst scenario: ",
              conditionMessage(nmax_res))
      nmax_res <- NULL
    }
  }

  structure(list(scenarios = scenarios,
                 worst = scenarios[worst_i, , drop = FALSE],
                 marginals = marginals, n = as.integer(n), claim = claim,
                 sesoi = sesoi, conf = conf, nsim = as.integer(nsim),
                 pin = pin, estimand = dgp$estimand %||% "ape",
                 target = dgp$target_est,
                 nmax = nmax_res, nmax_power = nmax_power),
            class = "powerape_robust")
}

#' @export
print.powerape_robust <- function(x, ...) {
  cat(sprintf("powerape robustness sweep -- %s claim, %s, pin = %s\n",
              x$claim, toupper(x$estimand), x$pin))
  cat(sprintf("  n = %d, nsim = %d per scenario, %d scenario(s) over: %s\n",
              x$n, x$nsim, nrow(x$scenarios),
              paste(names(x$marginals), collapse = ", ")))
  ok <- !is.na(x$scenarios$power)
  cat(sprintf("  power range [%.3f, %.3f]; worst scenario:\n",
              min(x$scenarios$power[ok]), max(x$scenarios$power[ok])))
  print(x$worst[, setdiff(names(x$worst), "note"), drop = FALSE],
        row.names = FALSE)
  if (any(!ok))
    cat(sprintf("  %d infeasible scenario(s); see `$scenarios$note`.\n", sum(!ok)))
  cat("  scenarios:\n")
  print(x$scenarios[, setdiff(names(x$scenarios), "note"), drop = FALSE],
        row.names = FALSE)
  if (!is.null(x$nmax))
    cat(sprintf("  n_max: n = %d reaches %.0f%% power in the worst scenario (achieved %.3f, MCSE %.3f).\n",
                x$nmax$n, 100 * x$nmax_power, x$nmax$power, x$nmax$mcse))
  invisible(x)
}
