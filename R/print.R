# Print methods ----------------------------------------------------------------

pa_var_label <- function(v) {
  base <- if (v$type == "binary") sprintf("%s (binary, p = %.2f", v$name, v$p)
          else sprintf("%s (normal, mean %.2f, sd %.2f", v$name, v$mean, v$sd)
  if (!is.null(v$icc)) base <- sprintf("%s, icc = %.2f", base, v$icc)
  paste0(base, ")")
}

#' @export
print.powerape_dgp <- function(x, ...) {
  cat(sprintf("powerape DGP -- %s%s\n", x$model,
              if (x$route == "empirical") " (empirical covariates)" else ""))
  focal_note <- if (x$route == "empirical") {
    if (!is.null(x$emp_xd)) " [resampled jointly with covariates]"
    else " [drawn independently of covariates]"
  } else ""
  cat(sprintf("  focal: %s%s\n", pa_var_label(x$focal), focal_note))
  if (!is.null(x$moderator))
    cat(sprintf("  moderator: %s (focal-by-moderator interaction included)\n",
                pa_var_label(x$moderator)))
  if (x$k > 0L) {
    labs <- if (x$route == "empirical") {
      paste(x$covariates, collapse = ", ")
    } else {
      paste(vapply(x$covariates, pa_var_label, character(1)), collapse = ", ")
    }
    cat("  covariates:", labs, "\n")
  } else {
    cat("  covariates: none\n")
  }
  if (identical(x$route, "panel"))
    cat(sprintf("  panel: %s periods per unit; rho = %.2f (latent unit share), cre_share = %.2f; CRE %s w/ Mundlak means%s, unit-clustered SEs\n",
                panel_len_label(x), x$rho, x$cre_share, x$model,
                if (isTRUE(x$unbalanced)) " over observed periods + T-cohort dummies (Wooldridge 2019)" else ""))
  if (identical(x$route, "iv"))
    cat(sprintf("  iv: endogenous focal (%s first stage; instruments: %s); endogeneity rho = %.2f, first-stage strength = %.2f; control-function probit, stacked robust SEs\n",
                x$first_stage,
                paste(vapply(x$instruments, function(v) v$name, character(1)),
                      collapse = ", "),
                x$endogeneity, x$iv_strength))
  if (x$route == "empirical")
    cat(sprintf("  pilot rows: %d (resampled with replacement)\n", x$n_emp))
  if (identical(x$builder, "from_fit"))
    cat(sprintf("  pilot model: %s (nuisance coefficients lifted; focal coefficient re-solved)\n",
                x$formula_str))
  cat(sprintf("  baseline P(Y=1 | reference) = %.3f, nuisance signal (latent pseudo-R2) = %.3f\n",
              x$baseline, x$signal_implied))
  if (!is.null(x$target_est)) {
    if (identical(x$estimand, "aie")) {
      cat(sprintf("  true AIE pinned at %+.4f (beta_int = %.4f)\n",
                  x$target_est, x$beta_int))
      cat(sprintf("  main APEs at reference: %s %+.3f (beta = %.4f), %s %+.3f (beta = %.4f)\n",
                  x$focal$name, x$main_focal, x$beta_focal,
                  x$moderator$name, x$main_moderator, x$beta_mod))
    } else if (x$focal$type == "binary") {
      cat(sprintf("  true APE pinned at %+.4f (beta_focal = %.4f, implied P(Y=1 | focal=1) = %.3f)\n",
                  x$target_est, x$beta_focal, x$p1))
    } else {
      cat(sprintf("  true APE pinned at %+.4f per unit of %s (beta_focal = %.4f, reference = %.2f)\n",
                  x$target_est, x$focal$name, x$beta_focal, x$focal_ref))
    }
  } else {
    cat("  target effect not set yet - call set_ape() or set_aie().\n")
  }
  invisible(x)
}

#' @export
print.powerape_power <- function(x, ...) {
  claim_lab <- switch(x$claim,
    minimum = sprintf("minimum-effect claim (CI lower bound > %.3f)", x$sesoi),
    detect = "detection claim (CI excludes 0, directional)",
    equivalence = sprintf("equivalence claim (CI within +/-%.3f)", x$sesoi))
  cat(sprintf("powerape -- %s\n", claim_lab))
  n_lab <- if (identical(x$dgp$route, "panel")) {
    sprintf("n = %d units x %s periods (%s obs)", x$n,
            panel_len_label(x$dgp), panel_obs_label(x$dgp, x$n))
  } else {
    sprintf("n = %d", x$n)
  }
  se_lab <- if (identical(x$dgp$route, "iv")) ", stacked robust SEs"
            else if (identical(x$se, "robust")) ", robust SEs"
            else ""
  cat(sprintf("  %s, %s, assumed true %s %+.4f, %.0f%% CI, nsim = %d%s\n",
              x$model, n_lab, toupper(x$estimand), x$target, 100 * x$conf,
              x$nsim, se_lab))
  cat(sprintf("  power = %.3f (MCSE %.3f)\n", x$power, x$mcse))
  o <- x$outcomes
  cat("  outcomes:", paste(sprintf("%s %.3f", gsub("_", "-", names(o)), o),
                           collapse = " | "), "\n")
  if (o[["failed"]] > 0)
    cat(sprintf("  note: %.1f%% of replications failed to fit; they count against power.\n",
                100 * o[["failed"]]))
  invisible(x)
}

#' @export
print.powerape_curve <- function(x, ...) {
  cat(sprintf("powerape power curve -- %s claim, %s (%s), nsim = %d per point\n",
              x$claim, x$model, toupper(x$estimand), x$nsim))
  print(x$results, row.names = FALSE)
  invisible(x)
}

#' @export
print.powerape_n <- function(x, ...) {
  cat(sprintf("powerape required sample size -- %s claim\n", x$claim))
  lab <- if (isTRUE(x$confirmed)) "confirmed" else "achieved"
  unit_lab <- if (identical(x$dgp$route, "panel")) {
    sprintf("n = %d units (x %s periods)", x$n, panel_len_label(x$dgp))
  } else {
    sprintf("n = %d", x$n)
  }
  cat(sprintf("  %s for %.0f%% target power (%s %.3f, MCSE %.3f)\n",
              unit_lab, 100 * x$goal, lab, x$power, x$mcse))
  cat(sprintf("  assumed true %s %+.4f, sesoi %s, %.0f%% CI, %s\n",
              toupper(x$estimand), x$target,
              if (is.null(x$sesoi)) "-" else sprintf("%.3f", x$sesoi),
              100 * x$conf, x$model))
  n_search <- sum(x$history$stage == "search")
  n_conf <- sum(x$history$stage == "confirm")
  if (isTRUE(x$confirm_requested)) {
    if (isTRUE(x$confirmed)) {
      cat(sprintf("  search: %d step(s); confirmed in %d round(s) at nsim = %d.\n",
                  n_search, n_conf, x$nsim_confirm))
    } else {
      cat(sprintf("  search: %d step(s); confirmation (%d round(s), nsim = %d) fell short -- treat as a lower bound.\n",
                  n_search, n_conf, x$nsim_confirm))
    }
  } else {
    cat(sprintf("  search: %d step(s); confirm with ape_power() at a larger nsim.\n",
                n_search))
  }
  invisible(x)
}
