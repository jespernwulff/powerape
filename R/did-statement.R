# Methods paragraph for difference-in-differences designs ----------------------
# power_statement() dispatches here for designs built with ape_dgp_did().
# Both analyses are described: the pooled
# probit/logit with the effect on the treated as an average partial effect,
# and the linear probability analysis.

# A result whose analysis assumes another parallel-trends scale than the
# simulated world has is a rejection rate, not power: the functions that
# re-read or report stored power stop on it.
did_power_only <- function(dgp, what) {
  if (is_did(dgp) && !did_valid(dgp))
    stop(sprintf(paste(
      "%s needs power, and this result is a rejection rate: %s, so the",
      "analysis is biased for the effect on the treated. Use a design whose",
      "`trends` and `analysis` are on the same scale."),
      what, did_mismatch_text(dgp)), call. = FALSE)
  invisible(TRUE)
}

# "one period" / "3 periods"
did_n_periods <- function(k) if (k == 1L) "one period" else sprintf("%d periods", k)

# The assumed world of a difference-in-differences design, in words.
describe_did <- function(d) {
  panel <- identical(d$design, "panel")
  r <- d$rates
  sub <- did_has_subgroup(d)
  what <- if (panel) "units" else "observations"
  s_design <- sprintf(paste0("a difference-in-differences design with a treated ",
                             "group (%.0f%% of the %s) and a control group, ",
                             "observed in %s before and %s after one common ",
                             "treatment date"),
                      100 * d$p_treated, what, did_n_periods(d$n_pre),
                      did_n_periods(d$n_post))
  s_data <- if (panel) {
    sprintf(paste0("as a panel of the same units, with a latent within-unit ",
                   "correlation of %.2f%s shared by all periods (Gaussian copula)"),
            d$rho, if (did_cov_dep(d)) " given the covariates," else "")
  } else {
    sprintf(paste0("as repeated cross-sections, with %.0f%% of the observations ",
                   "after the treatment date"), 100 * d$p_post)
  }
  rate_txt <- function(m) {
    sprintf(paste0("%.3f and %.3f in the control group and %.3f and %.3f in ",
                   "the treated group (%.3f without treatment)"),
            r$c00[m], r$c01[m], r$c10[m], r$c11[m], r$cf[m])
  }
  s_rates <- if (sub) {
    nm <- d$subgroup$name
    sprintf(paste0("a binary subgroup variable %s (prevalence %.2f); outcome ",
                   "rates before and after treatment of %s for %s = 0, and of ",
                   "%s for %s = 1"),
            nm, d$subgroup$p, rate_txt(1L), nm, rate_txt(2L), nm)
  } else {
    sprintf("outcome rates before and after treatment of %s", rate_txt(1L))
  }
  s_trend <- sprintf("untreated trends parallel %s%s",
                     if (identical(d$trends, "index")) sprintf("on the %s index", d$model)
                     else "in probabilities",
                     if (sub) " within each subgroup" else "")
  var_desc <- function(v) {
    sprintf("%s (%s)", v$name,
            if (v$type == "binary") sprintf("binary, prevalence %.2f", v$p)
            else sprintf("continuous, mean %.3g, SD %.3g", v$mean, v$sd))
  }
  s_cov <- if (d$k > 0L) {
    sig <- if (identical(d$link, "logit") && isTRUE(d$s2 > 0)) {
      sprintf(paste0("contributing %.2f of the latent-index variance at unit ",
                     "error variance (a McKelvey-Zavoina R-squared of %.3f on ",
                     "the logistic scale)"), d$signal_implied, mz_r2_logit(d$s2))
    } else {
      sprintf("contributing a latent pseudo-R-squared of %.2f", d$signal_implied)
    }
    rr <- d$spec$correlation
    off <- if (is.matrix(rr)) rr[upper.tri(rr)] else rr
    dep <- if (is.null(rr) || !length(off) || all(abs(off) < 1e-12)) {
      "independent of group membership"
    } else if (is.matrix(rr)) {
      sprintf(paste0("with a latent (Gaussian-copula) correlation matrix over ",
                     "group membership and the covariates (off-diagonal ",
                     "entries from %.2f to %.2f)"), min(off), max(off))
    } else {
      sprintf(paste0("with an exchangeable latent (Gaussian-copula) correlation ",
                     "of %.2f among group membership and the covariates"), rr)
    }
    sprintf("%scovariates %s, %s, %s",
            if (panel) "time-constant " else "",
            paste(vapply(d$covariates, var_desc, character(1)), collapse = ", "),
            dep, sig)
  } else {
    "no covariates"
  }
  paste(c(paste(s_design, s_data), s_rates, s_trend, s_cov, did_pilot_clause(d)),
        collapse = "; ")
}

# How the simulated studies of a difference-in-differences design were
# analysed, in words (one clause, no full stop).
did_analysis_text <- function(x) {
  d <- x$dgp
  sub <- did_has_subgroup(d)
  regs <- paste0("group, period indicators and the treatment indicator",
                 if (d$k > 0L) ", and the covariates" else "",
                 if (sub) sprintf(", fully interacted with %s", d$subgroup$name) else "")
  if (is_linear(d)) {
    sprintf(paste0("ordinary least squares of the outcome on %s (a linear ",
                   "probability model), with the %s as the coefficient of %s ",
                   "and %s standard errors"),
            regs,
            if (sub) "difference between the subgroups' effects" else "effect on the treated",
            if (sub) sprintf("the treatment indicator's product with %s", d$subgroup$name)
            else "the treatment indicator",
            if (unit_design(d)) "unit-clustered" else "heteroskedasticity-robust (HC1)")
  } else {
    se_txt <- if (unit_design(d)) "unit-clustered standard errors"
              else if (identical(x$se %||% "model", "robust")) "heteroskedasticity-robust (sandwich) standard errors"
              else "model-based standard errors"
    sprintf(paste0("a pooled %s regression of the outcome on %s, with %s ",
                   "computed as %s the average partial effect of the treatment ",
                   "indicator over the treated post-treatment observations ",
                   "(Wooldridge, 2023), and delta-method Wald confidence ",
                   "intervals from %s"),
            d$model, regs,
            if (sub) "the estimate" else "the effect on the treated",
            if (sub) "the difference between the two subgroups in" else "",
            se_txt)
  }
}

# power_statement() for ape_power(), ape_n() and ape_mde() results on a
# difference-in-differences design.
did_statement <- function(x) {
  d <- x$dgp
  rule <- !is.null(x$alpha_rule)
  sub <- did_has_subgroup(d)
  lin <- is_linear(d)
  est_short <- if (sub) "difference in ATTs" else "ATT"
  est_long <- if (sub) {
    sprintf(paste0("difference between two subgroups' average treatment ",
                   "effects on the treated (by %s)"), d$subgroup$name)
  } else {
    "average treatment effect on the treated (ATT)"
  }
  s1 <- sprintf(paste0("We conducted a simulation-based power analysis for the ",
                       "%s in a difference-in-differences design, in ",
                       "probability points, using the powerape package ",
                       "(version %s), following the confidence-interval ",
                       "approach of Riesthuis (2024)."),
                est_long, as.character(packageVersion("powerape")))
  s2 <- sprintf("The assumed data-generating process was %s.", describe_did(d))
  s3 <- if (inherits(x, "powerape_mde")) {
    NULL
  } else if (sub) {
    sprintf(paste0("The assumed effects on the treated were %.3f (%s = 0) and ",
                   "%.3f (%s = 1), a true difference of %.3f (%.1f percentage ",
                   "points)."),
            d$att_sub[1], d$subgroup$name, d$att_sub[2], d$subgroup$name,
            x$target, 100 * x$target)
  } else {
    sprintf("The assumed true ATT was %.3f (%.1f percentage points).",
            x$target, 100 * x$target)
  }
  n_text <- if (unit_design(d)) {
    sprintf("n = %d units observed in %d periods (%s unit-period observations)",
            x$n, d$n_periods, format(x$n * d$n_periods, big.mark = ","))
  } else {
    sprintf("n = %d observations", x$n)
  }
  s4 <- statement_result(x, n_text, est_short, rule,
                         if (unit_design(d)) "number of units" else "total sample size")
  s5 <- statement_outcomes(x)
  keep <- identical(x$separation, "keep")
  sepd <- x$separated %||% 0
  cell <- cell_label(d$estimand)
  an <- did_analysis_text(x)
  s6 <- if (lin) {
    ## ordinary least squares has no separation
    if (inherits(x, "powerape_power")) {
      sprintf(paste0("Each simulated study was analysed by %s; replications ",
                     "whose fit failed (%.1f%%) counted against the claim."),
              an, 100 * x$outcomes[["failed"]])
    } else {
      sprintf(paste0("Each simulated study was analysed by %s; replications ",
                     "whose fit failed counted against the claim."), an)
    }
  } else if (inherits(x, "powerape_power")) {
    other <- x$outcomes[["failed"]] - (if (keep) 0 else sepd)
    if (keep) {
      sprintf(paste0("Each simulated study was analysed by %s; replications ",
                     "whose model fit failed (%.1f%%) counted against the ",
                     "claim, and replications with separation (a %s with no ",
                     "events or no non-events; %.1f%%) were retained with ",
                     "their Wald intervals."),
              an, 100 * other, cell, 100 * sepd)
    } else {
      sprintf(paste0("Each simulated study was analysed by %s; replications ",
                     "without a maximum-likelihood estimate, from a failed fit ",
                     "or from separation (a %s with no events or no ",
                     "non-events), counted against the claim (%.1f%% of ",
                     "replications%s)."),
              an, cell, 100 * x$outcomes[["failed"]],
              if (sepd > 0) sprintf(", %.1f%% with separation", 100 * sepd) else "")
    }
  } else {
    sprintf(paste0("Each simulated study was analysed by %s; replications ",
                   "without a maximum-likelihood estimate (failed fits%s) ",
                   "counted against the claim%s."),
            an,
            if (keep) "" else sprintf(" and separation, a %s with no events or no non-events", cell),
            if (keep) ", while replications with separation were retained with their Wald intervals" else "")
  }
  structure(paste(c(s1, s2, s3, s4, s5, s6), collapse = " "),
            class = "powerape_statement")
}
