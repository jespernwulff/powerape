# Robust power / n_max mode (Hancock & Feng, 2025) -----------------------------

# Rebuild a DGP from its stored constructor spec with some inputs replaced.
rebuild_dgp <- function(dgp, spec) {
  suppressWarnings(
    switch(dgp$builder %||% "parametric",
           from_fit = do.call(ape_dgp_from_fit, spec),
           empirical = do.call(ape_dgp_empirical, spec),
           panel = do.call(ape_dgp_panel, spec),
           iv = do.call(ape_dgp_iv, spec),
           do.call(ape_dgp, spec))
  )
}

# The claim's distance from its null boundary at an assumed effect `t`.
claim_distance <- function(claim, t, sesoi) {
  switch(claim, detect = abs(t), minimum = abs(t) - sesoi,
         equivalence = sesoi - abs(t))
}

#' Robustness sweep over contextual assumptions (worst-case power and n_max)
#'
#' Power computed from a single guess at the contextual inputs -- baseline
#' rate, nuisance signal, covariate dependence -- is fragile (Hancock &
#' Feng, 2025). `ape_robust()` re-runs the power analysis over a scenario
#' grid of those inputs and reports per-scenario power, the worst case, and
#' (optionally) `n_max`: the sample size that reaches the target power in
#' every scenario -- the insurance-premium n.
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
#' **n_max.** The insurance-premium n is the largest requirement across
#' scenarios. Power at the user's `n` cannot rank the scenarios for that
#' purpose -- at a generous `n` every scenario's power is near 1 -- so the
#' scenarios are ranked by a criterion that does not saturate: the
#' normal-approximation requirement implied by each scenario's simulated
#' standard error and its effect's distance from the claim boundary.
#' [ape_n()] then runs in the least favorable scenario and in any scenario
#' within 5% of it on that criterion (up to three), and `n_max` is the
#' largest confirmed answer. When the implied effect of a feasible
#' scenario sits on or beyond the claim boundary (possible with `pin =
#' "coefficients"`), no n reaches the target there, and `n_max` is not
#' reported.
#'
#' @inheritParams ape_power
#' @param vary Named list of scenario values for contextual inputs.
#'   Parametric route: `baseline`, `signal`, `correlation` (scalar);
#'   empirical route: `baseline`, `signal`. A length-2 element is expanded
#'   to `grid_points` equally spaced values; longer vectors are used as-is.
#'   Scenarios are the full factorial grid.
#' @param pin `"ape"` (re-invert per scenario) or `"coefficients"` (hold
#'   coefficients, let the implied effect drift). Ignored in MDE mode,
#'   which always re-solves.
#' @param mode `"power"` (default) reports power at the pinned effect per
#'   scenario; `"mde"` reports the **minimum detectable effect** per
#'   scenario via [ape_mde()] -- the worst case is then the MDE largest
#'   in magnitude, answering "what is the smallest effect this design
#'   finds even under the least favorable contextual assumptions?".
#' @param power Target power for MDE mode (default 0.80; unused in power
#'   mode).
#' @param grid_points Grid size used to expand length-2 `vary` elements.
#' @param nmax Run the n_max search (default TRUE; power mode only).
#' @param nmax_power Target power for the n_max search (default 0.90, the
#'   default of [ape_n()]).
#' @param direction MDE mode only: `"positive"` or `"negative"` search
#'   direction, as in [ape_mde()]; default: the sign of the base DGP's
#'   pinned effect.
#'
#' @return A `powerape_robust` object: `scenarios` (inputs, implied effect,
#'   power, MCSE, share of simulated studies with separation), the worst
#'   scenario, marginal mean power per input, and the n_max result
#'   (`nmax`, the [ape_n()] result in the scenario that sets it;
#'   `nmax_scenario`, that scenario's inputs; `nmax_candidates`, every
#'   scenario searched with its answer).
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
                       sesoi = NULL, alpha = 0.05, conf = NULL, nsim = 800,
                       seed = NULL, vary, pin = c("ape", "coefficients"),
                       mode = c("power", "mde"), power = 0.80,
                       grid_points = 3, nmax = TRUE, nmax_power = 0.90,
                       se = c("model", "robust"), direction = NULL,
                       separation = c("fail", "keep")) {
  claim <- match.arg(claim)
  pin <- match.arg(pin)
  mode <- match.arg(mode)
  se <- match.arg(se)
  separation <- match.arg(separation)
  if (!identical(se, "model") && dgp$route %in% c("panel", "iv")) {
    warning(sprintf(paste("`se` is fixed by the route: panel designs use",
                          "unit-clustered SEs, IV designs the stacked robust",
                          "sandwich; `se = \"%s\"` is ignored."), se))
    se <- "model"
  }
  stopifnot(is.numeric(n), length(n) == 1L, n >= 20)
  ## n is fixed across scenarios, so a sample-size rule resolves once
  lv <- claim_levels(claim, alpha, conf, n)
  conf <- lv$conf_claim
  rule <- is.function(alpha)
  if (mode == "mde" && claim == "equivalence")
    stop(paste("MDE mode searches the effect for detect/minimum claims; for",
               "the equivalence analog run ape_mde(claim = \"equivalence\")",
               "per scenario directly."), call. = FALSE)
  if (mode == "power") check_coherence(dgp, claim, sesoi, conf)
  if (mode == "mde" && !is.null(dgp$moderator) &&
      (is.null(dgp$main_focal) || is.null(dgp$main_moderator)))
    stop("MDE mode for an AIE design needs the DGP pinned with set_aie() first.",
         call. = FALSE)
  direction <- direction %||%
    (if (!is.null(dgp$target_est) && dgp$target_est < 0) "negative" else "positive")
  direction <- match.arg(direction, c("positive", "negative"))
  stopifnot(is.list(vary), length(vary) >= 1L)
  if (is.null(names(vary)) || !all(nzchar(names(vary))))
    stop("`vary` must be a fully named list.", call. = FALSE)
  allowed <- switch(dgp$builder %||% "parametric",
                    from_fit = "baseline",
                    empirical = c("baseline", "signal"),
                    panel = c("baseline", "signal", "rho", "retention"),
                    iv = c("baseline", "signal", "endogeneity", "iv_strength"),
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
  crit <- rep(NA_real_, nrow(grid))  # SE-based requirement ranking (power mode)
  for (i in seq_len(nrow(grid))) {
    spec <- dgp$spec
    for (nm in names(grid)) spec[[nm]] <- grid[i, nm]
    d_i <- tryCatch(rebuild_dgp(dgp, spec), error = function(e) e)
    if (mode == "mde") {
      if (inherits(d_i, "error")) {
        rows[[i]] <- data.frame(grid[i, , drop = FALSE], mde = NA_real_,
                                power = NA_real_, mcse = NA_real_,
                                separated = NA_real_,
                                note = conditionMessage(d_i),
                                stringsAsFactors = FALSE)
        next
      }
      si <- if (is.null(seed)) NULL else seed + i
      m_i <- tryCatch(
        suppressWarnings(
          ape_mde(d_i, n, power = power, claim = claim, sesoi = sesoi,
                  conf = conf, nsim = nsim, seed = si,
                  main_focal = dgp$main_focal,
                  main_moderator = dgp$main_moderator, confirm = FALSE,
                  se = se, direction = direction, separation = separation)),
        error = function(e) e)
      if (inherits(m_i, "error")) {
        rows[[i]] <- data.frame(grid[i, , drop = FALSE], mde = NA_real_,
                                power = NA_real_, mcse = NA_real_,
                                separated = NA_real_,
                                note = conditionMessage(m_i),
                                stringsAsFactors = FALSE)
      } else {
        rows[[i]] <- data.frame(grid[i, , drop = FALSE], mde = m_i$mde,
                                power = m_i$power, mcse = m_i$mcse,
                                separated = m_i$separated, note = "",
                                stringsAsFactors = FALSE)
      }
      next
    }
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
                              separated = NA_real_,
                              note = conditionMessage(d_i),
                              stringsAsFactors = FALSE)
      next
    }
    si <- if (is.null(seed)) NULL else seed + i
    pw <- power_once(d_i, n, claim, sesoi, conf, nsim, si, enforce = FALSE,
                     se = se, keep_draws = TRUE, separation = separation)
    rows[[i]] <- data.frame(grid[i, , drop = FALSE],
                            implied_effect = d_i$target_est,
                            power = pw$power, mcse = pw$mcse,
                            separated = pw$separated, note = "",
                            stringsAsFactors = FALSE)
    dgps[[i]] <- d_i
    ## requirement criterion: (SE scale c = se * sqrt(n)) / (distance to
    ## the claim boundary); (c (z + z_power) / dist)^2 is the normal-
    ## approximation n, so the ranking by c / dist is free of alpha and of
    ## the target power -- and, unlike power at n, it never saturates
    ok_d <- pw$draws$ok
    dist <- claim_distance(claim, d_i$target_est, sesoi)
    if (any(ok_d))
      crit[i] <- if (dist > 0) mean(pw$draws$se[ok_d]) * sqrt(n) / dist else Inf
  }
  scenarios <- do.call(rbind, rows)
  rownames(scenarios) <- NULL

  key <- if (mode == "mde") abs(scenarios$mde) else scenarios$power
  ok_idx <- which(!is.na(key))
  if (!length(ok_idx))
    stop("No scenario was feasible; widen or shift the `vary` ranges.", call. = FALSE)
  worst_i <- if (mode == "mde") {
    ok_idx[which.max(key[ok_idx])]
  } else {
    ## lowest power; ties (e.g. every scenario saturated at 1) go to the
    ## scenario with the least favorable precision
    cr_ok <- crit[ok_idx]
    cr_ok[is.na(cr_ok)] <- -Inf
    ok_idx[order(scenarios$power[ok_idx], -cr_ok)[1L]]
  }

  marginals <- lapply(names(vary), function(nm) {
    ag <- aggregate(key[ok_idx],
                    by = list(value = scenarios[[nm]][ok_idx]), FUN = mean)
    names(ag) <- c("value", if (mode == "mde") "mean_mde" else "mean_power")
    if (mode == "mde" && identical(direction, "negative")) ag$mean_mde <- -ag$mean_mde
    ag
  })
  names(marginals) <- names(vary)

  nmax_res <- NULL
  nmax_i <- NA_integer_
  nmax_cand <- NULL
  if (nmax && mode == "power") {
    cr <- crit[ok_idx]
    unsolvable <- ok_idx[!is.na(cr) & !is.finite(cr)]
    if (length(unsolvable)) {
      vals <- unlist(grid[unsolvable[1L], , drop = FALSE])
      warning(sprintf(paste(
        "n_max not reported: in %d scenario(s) the implied effect sits on or",
        "beyond the claim boundary (e.g. %s), so no sample size reaches the",
        "target power there."),
        length(unsolvable),
        paste(sprintf("%s = %s", names(vals), format(vals)), collapse = ", ")),
        call. = FALSE)
    } else if (all(is.na(cr))) {
      warning(paste("n_max not reported: no scenario produced usable estimates",
                    "at this n to rank the scenarios by."), call. = FALSE)
    } else {
      o <- ok_idx[order(-cr)]
      o <- o[!is.na(crit[o])]
      cand <- o[crit[o] >= 0.95 * crit[o[1L]]]
      cand <- cand[seq_len(min(3L, length(cand)))]
      runs <- lapply(seq_along(cand), function(k) {
        sk <- if (is.null(seed)) NULL else seed + nrow(grid) + 1L + 1000L * (k - 1L)
        tryCatch(
          ape_n(dgps[[cand[k]]], power = nmax_power, claim = claim, sesoi = sesoi,
                alpha = if (rule) alpha else lv$alpha,
                conf = if (rule) NULL else conf, nsim = nsim, seed = sk,
                se = se, separation = separation),
          error = function(e) e)
      })
      failed <- vapply(runs, inherits, logical(1), "error")
      if (all(failed)) {
        warning("n_max search failed in the least favorable scenario: ",
                conditionMessage(runs[[1L]]))
      } else {
        nn <- vapply(runs, function(r) if (inherits(r, "error")) NA_real_ else r$n,
                     numeric(1))
        j <- which.max(nn)
        nmax_res <- runs[[j]]
        nmax_i <- cand[j]
        nmax_cand <- data.frame(grid[cand, , drop = FALSE],
                                criterion = crit[cand], n = nn,
                                row.names = NULL)
        if (any(failed))
          warning("n_max search failed in a near-least-favorable scenario: ",
                  conditionMessage(runs[[which(failed)[1L]]]))
      }
    }
  }

  structure(list(scenarios = scenarios,
                 worst = scenarios[worst_i, , drop = FALSE],
                 marginals = marginals, n = as.integer(n), claim = claim,
                 sesoi = sesoi, conf = conf, alpha = lv$alpha,
                 alpha_rule = if (rule) alpha else NULL,
                 nsim = as.integer(nsim), se = se, separation = separation,
                 pin = pin, mode = mode, goal = power,
                 direction = if (mode == "mde") direction else NULL,
                 estimand = if (!is.null(dgp$moderator)) "aie"
                            else dgp$estimand %||% "ape",
                 target = dgp$target_est,
                 nmax = nmax_res, nmax_power = nmax_power,
                 nmax_scenario = if (is.na(nmax_i)) NULL
                                 else scenarios[nmax_i, , drop = FALSE],
                 nmax_candidates = nmax_cand),
            class = "powerape_robust")
}

#' @export
print.powerape_robust <- function(x, ...) {
  if (identical(x$mode, "mde")) {
    cat(sprintf("powerape robustness sweep -- minimum detectable %s, %s claim%s\n",
                toupper(x$estimand), x$claim,
                if (identical(x$direction, "negative")) ", decreases" else ""))
    cat(sprintf("  n = %d, %.0f%% target power, nsim = %d per scenario, %d scenario(s) over: %s\n",
                x$n, 100 * x$goal, x$nsim, nrow(x$scenarios),
                paste(names(x$marginals), collapse = ", ")))
    ok <- !is.na(x$scenarios$mde)
    cat(sprintf("  MDE range [%.4f, %.4f]; worst (largest-MDE) scenario:\n",
                min(x$scenarios$mde[ok]), max(x$scenarios$mde[ok])))
    print(x$worst[, setdiff(names(x$worst), "note"), drop = FALSE],
          row.names = FALSE)
    if (any(!ok))
      cat(sprintf("  %d infeasible scenario(s); see `$scenarios$note`.\n", sum(!ok)))
    cat("  scenarios:\n")
    print(x$scenarios[, setdiff(names(x$scenarios), "note"), drop = FALSE],
          row.names = FALSE)
    return(invisible(x))
  }
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
  if (!is.null(x$nmax)) {
    sc <- x$nmax_scenario
    inputs <- names(x$marginals)
    sc_lab <- if (is.null(sc)) "" else
      paste0(" (", paste(sprintf("%s = %s", inputs, format(unlist(sc[inputs]))),
                         collapse = ", "), ")")
    cat(sprintf("  n_max: n = %d reaches %.0f%% power in every scenario; set by the least favorable one%s, achieved %.3f (MCSE %.3f).\n",
                x$nmax$n, 100 * x$nmax_power, sc_lab, x$nmax$power, x$nmax$mcse))
  }
  if (any(x$scenarios$separated[ok] > 0.01, na.rm = TRUE))
    cat(sprintf("  note: separation in more than 1%% of simulated studies in %d scenario(s) (column `separated`; separation = \"%s\").\n",
                sum(x$scenarios$separated[ok] > 0.01, na.rm = TRUE), x$separation))
  invisible(x)
}
