# Exact finite-sample power of the engine's decision rule in the saturated
# case (binary focal, no covariates), by triple binomial enumeration.
# Shared by test-exact-power.R and test-ape-n-confirm.R; mirrored in
# validation/run-validation.R (V1).
#
# `separation` mirrors the engine's convention for an arm with no events or
# no non-events (separation; DESIGN.md section 16): "fail" (the engine's
# default since 1.11.0) scores it as a failed fit; "keep" scores its
# degenerate Wald interval, in which the boundary arm contributes zero
# variance. A constant focal column and constant y (both arms at the same
# boundary) are failed fits under both conventions. Under "keep" a
# completely separated sample (one arm all 0, the other all 1) is scored
# with a zero SE, as the probit engine's ~1e-7 SE effectively does; its
# mass is negligible in the designs tested here.
exact_power_sat <- function(n, pd, p0, p1, conf, claim, sesoi = NA,
                            separation = c("fail", "keep")) {
  separation <- match.arg(separation)
  z <- qnorm(1 - (1 - conf) / 2)
  total <- 0
  for (n1 in 0:n) {
    pn1 <- dbinom(n1, n, pd)
    if (pn1 < 1e-14) next
    n0 <- n - n1
    if (n1 == 0L || n0 == 0L) next          # constant focal column: failed fit
    if (separation == "fail") {
      x1 <- seq_len(n1 - 1L)                # interior counts only (else failed)
      x0 <- seq_len(n0 - 1L)
    } else {
      x1 <- 0:n1
      x0 <- 0:n0
    }
    if (!length(x1) || !length(x0)) next
    p1h <- x1 / n1
    p0h <- x0 / n0
    est <- outer(p1h, p0h, "-")
    se <- sqrt(outer(p1h * (1 - p1h) / n1, p0h * (1 - p0h) / n0, "+"))
    dec <- switch(claim,
                  detect = est - z * se > 0,
                  minimum = est - z * se > sesoi,
                  equivalence = (est + z * se < sesoi) & (est - z * se > -sesoi))
    if (separation == "keep")                 # constant y: failed fit
      dec <- dec & !(outer(x1 == 0, x0 == 0, "&") | outer(x1 == n1, x0 == n0, "&"))
    w <- outer(dbinom(x1, n1, p1), dbinom(x0, n0, p0))
    total <- total + pn1 * sum(w[dec])
  }
  total
}
