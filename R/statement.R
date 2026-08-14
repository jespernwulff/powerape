# Citable power statement ------------------------------------------------------

describe_dgp <- function(d) {
  var_desc <- function(v) {
    if (v$type == "binary") sprintf("%s (binary, prevalence %.2f)", v$name, v$p)
    else sprintf("%s (continuous, mean %.1f, SD %.1f)", v$name, v$mean, v$sd)
  }
  parts <- if (identical(d$route, "panel")) {
    sprintf(paste0("a correlated random effects %s model over %d periods per ",
                   "unit (latent unit-effect share rho = %.2f, %.0f%% of it ",
                   "loaded on regressor means), estimated by pooled %s with ",
                   "Mundlak means and unit-clustered standard errors, with ",
                   "focal variable %s"),
            d$model, d$n_periods, d$rho, 100 * d$cre_share, d$model,
            var_desc(d$focal))
  } else if (identical(d$route, "iv")) {
    sprintf(paste0("a probit model with an endogenous focal variable %s ",
                   "(error correlation %.2f with the first stage), ",
                   "instrumented by %s (%s first stage, instrument strength ",
                   "%.2f as a share of the focal's %svariance) and estimated ",
                   "by control-function probit with stacked ",
                   "method-of-moments robust standard errors (the estimator ",
                   "of Stata's cfprobit)"),
            var_desc(d$focal), d$endogeneity,
            paste(vapply(d$instruments, var_desc, character(1)),
                  collapse = ", "),
            d$first_stage, d$iv_strength,
            if (d$focal$type == "binary") "latent " else "")
  } else {
    sprintf("a %s model with focal variable %s", d$model, var_desc(d$focal))
  }
  if (!is.null(d$moderator))
    parts <- paste0(parts,
                    sprintf(", moderator %s, and their interaction",
                            var_desc(d$moderator)))
  covs <- if (identical(d$builder, "from_fit")) {
    sprintf(paste0("covariate rows and nuisance coefficients taken from a ",
                   "fitted pilot model (%s; %d pilot observations, ",
                   "resampled with replacement)"),
            d$formula_str, d$n_emp)
  } else if (d$route == "empirical") {
    sprintf("covariates resampled with replacement from %d pilot observations (%s)",
            d$n_emp, paste(d$covariates, collapse = ", "))
  } else if (d$k > 0L) {
    sprintf("parametric covariates (%s; Gaussian-copula dependence)",
            paste(vapply(d$covariates, function(v) v$name, character(1)),
                  collapse = ", "))
  } else {
    "no additional covariates"
  }
  sprintf(paste0("%s; %s; baseline outcome rate %.3f with the focal%s at ",
                 "reference; nuisance covariates contribute a latent ",
                 "pseudo-R-squared of %.2f"),
          parts, covs, d$baseline,
          if (is.null(d$moderator)) "" else " and moderator",
          d$signal_implied)
}

claim_text <- function(claim, sesoi, conf) {
  cl <- sprintf("%.0f%%", 100 * conf)
  switch(claim,
    minimum = sprintf(paste0("the minimum-effect claim (the %s confidence ",
                             "interval's lower bound exceeding the smallest ",
                             "effect size of interest, %.3f)"), cl, sesoi),
    detect = sprintf("the detection claim (the %s confidence interval excluding zero, directional)", cl),
    equivalence = sprintf("the equivalence claim (the %s confidence interval lying within +/-%.3f)", cl, sesoi))
}

#' Render a power analysis as a citable methods paragraph
#'
#' Turns a [ape_power()] or [ape_n()] result into a self-contained methods
#' paragraph stating the estimand, the full data-generating assumptions,
#' the claim and CI convention, and the Monte Carlo precision -- the
#' transparency artifact for grant applications and preregistrations.
#'
#' @param x A `powerape_power` or `powerape_n` object.
#' @return A character string of class `powerape_statement` (printed
#'   wrapped).
#' @examples
#' \donttest{
#' d <- ape_dgp(focal = pa_var("treat", "binary", p = 0.5), baseline = 0.30)
#' d <- set_ape(d, target = 0.10)
#' pw <- ape_power(d, n = 800, claim = "minimum", sesoi = 0.03,
#'                 nsim = 500, seed = 1)
#' power_statement(pw)
#' }
#' @export
power_statement <- function(x) {
  stopifnot(inherits(x, "powerape_power") || inherits(x, "powerape_n") ||
              inherits(x, "powerape_mde"))
  d <- x$dgp
  is_aie <- identical(x$estimand, "aie")
  est_long <- if (is_aie) "average interaction effect (AIE)" else "average partial effect (APE)"
  est_short <- if (is_aie) "AIE" else "APE"
  subject <- if (is_aie) sprintf("%s-by-%s", d$focal$name, d$moderator$name) else d$focal$name

  s1 <- sprintf(paste0("We conducted a simulation-based power analysis for the ",
                       "%s of %s using the powerape package (version %s), ",
                       "following the confidence-interval approach of ",
                       "Riesthuis (2024)."),
                est_long, subject, as.character(packageVersion("powerape")))
  s2 <- sprintf("The assumed data-generating process was %s.", describe_dgp(d))
  s3 <- if (inherits(x, "powerape_mde")) {
    NULL
  } else {
    sprintf("The assumed true %s was %.3f (%.1f percentage points)%s.",
            est_short, x$target, 100 * x$target,
            if (is_aie) sprintf(paste0(", with conditional-at-reference ",
                                       "main-effect APEs of %.3f (%s) and ",
                                       "%.3f (%s)"),
                                d$main_focal, d$focal$name,
                                d$main_moderator, d$moderator$name)
            else "")
  }
  is_panel <- identical(d$route, "panel")
  n_text <- if (is_panel) {
    sprintf("n = %d units observed over %d periods (%d unit-period observations)",
            x$n, d$n_periods, x$n * d$n_periods)
  } else {
    sprintf("n = %d", x$n)
  }
  s4 <- if (inherits(x, "powerape_mde")) {
    conf_frag <- if (isTRUE(x$confirmed)) {
      sprintf(paste0(", confirmed by a high-precision verification run ",
                     "(simulated power %.3f, Monte Carlo SE %.3f, %d ",
                     "replications)"),
              x$power, x$mcse, x$nsim_confirm)
    } else {
      sprintf(" (simulated power %.3f, Monte Carlo SE %.3f)", x$power, x$mcse)
    }
    if (isTRUE(x$is_margin)) {
      sprintf(paste0("At a sample size of %s, the tightest equivalence ",
                     "bounds establishable with %.0f%% power were +/-%.3f ",
                     "(%.1f percentage points)%s."),
              n_text, 100 * x$goal, x$mde, 100 * x$mde, conf_frag)
    } else {
      sprintf(paste0("At a sample size of %s, the minimum detectable %s at ",
                     "%.0f%% power for %s was %.3f (%.1f percentage ",
                     "points)%s."),
              n_text, est_short, 100 * x$goal,
              claim_text(x$claim, x$sesoi, x$conf), x$mde, 100 * x$mde,
              conf_frag)
    }
  } else if (inherits(x, "powerape_n")) {
    size_noun <- if (is_panel) "number of units" else "total sample size"
    if (isTRUE(x$confirmed)) {
      sprintf(paste0("The required %s for %.0f%% power for %s ",
                     "was %s, confirmed by a high-precision verification ",
                     "run (simulated power %.3f, Monte Carlo SE %.3f, %d ",
                     "replications)."),
              size_noun, 100 * x$goal, claim_text(x$claim, x$sesoi, x$conf),
              n_text, x$power, x$mcse, x$nsim_confirm)
    } else {
      sprintf(paste0("The required %s for %.0f%% power for %s ",
                     "was %s (simulated power %.3f, Monte Carlo SE %.3f, ",
                     "%d replications per search step)."),
              size_noun, 100 * x$goal, claim_text(x$claim, x$sesoi, x$conf),
              n_text, x$power, x$mcse, x$nsim)
    }
  } else {
    sprintf(paste0("At a sample size of %s, simulated power for %s ",
                   "was %.3f (Monte Carlo SE %.3f; %d replications)."),
            n_text, claim_text(x$claim, x$sesoi, x$conf), x$power, x$mcse, x$nsim)
  }
  s5 <- if (inherits(x, "powerape_power") && !is.null(x$sesoi)) {
    o <- x$outcomes
    sprintf(paste0("Across replications, the probability of concluding a ",
                   "meaningful effect was %.3f, of detection without ",
                   "meaningfulness %.3f, of an inconclusive result %.3f, ",
                   "and of equivalence %.3f."),
            o[["minimum"]], o[["detect_only"]], o[["inconclusive"]],
            o[["equivalence"]])
  } else {
    NULL
  }
  est_desc <- if (identical(d$route, "iv")) {
    paste0("two-step control-function estimation with stacked ",
           "method-of-moments robust standard errors and delta-method ",
           "Wald confidence intervals")
  } else if (identical(x$se %||% "model", "robust")) {
    paste0("maximum likelihood with heteroskedasticity-robust (sandwich) ",
           "standard errors and delta-method Wald confidence intervals")
  } else {
    "maximum likelihood with delta-method Wald confidence intervals"
  }
  s6 <- if (inherits(x, "powerape_power")) {
    sprintf(paste0("Estimation used %s; replications that failed to ",
                   "converge (%.1f%%) counted against the claim."),
            est_desc, 100 * x$outcomes[["failed"]])
  } else {
    sprintf("Estimation used %s.", est_desc)
  }

  structure(paste(c(s1, s2, s3, s4, s5, s6), collapse = " "),
            class = "powerape_statement")
}

#' @export
print.powerape_statement <- function(x, ...) {
  cat(strwrap(unclass(x), width = 78), sep = "\n")
  invisible(x)
}
