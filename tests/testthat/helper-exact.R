# Exact finite-sample power of the engine's decision rule in the saturated
# case (binary focal, no covariates), by triple binomial enumeration.
# Shared by test-exact-power.R and test-ape-n-confirm.R; mirrored in
# validation/run-validation.R (V1).
exact_power_sat <- function(n, pd, p0, p1, conf, claim, sesoi = NA) {
  z <- qnorm(1 - (1 - conf) / 2)
  total <- 0
  for (n1 in 0:n) {
    pn1 <- dbinom(n1, n, pd)
    if (pn1 < 1e-14) next
    n0 <- n - n1
    if (n1 == 0L || n0 == 0L) next          # constant focal column: failed fit
    x1 <- 1:(n1 - 1)                        # interior counts only (else failed)
    x0 <- 1:(n0 - 1)
    p1h <- x1 / n1
    p0h <- x0 / n0
    est <- outer(p1h, p0h, "-")
    se <- sqrt(outer(p1h * (1 - p1h) / n1, p0h * (1 - p0h) / n0, "+"))
    dec <- switch(claim,
                  detect = est - z * se > 0,
                  minimum = est - z * se > sesoi,
                  equivalence = (est + z * se < sesoi) & (est - z * se > -sesoi))
    w <- outer(dbinom(x1, n1, p1), dbinom(x0, n0, p0))
    total <- total + pn1 * sum(w[dec])
  }
  total
}
