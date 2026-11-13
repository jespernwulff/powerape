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

# Exact power of the difference-in-differences engine for repeated
# cross-sections without covariates (test-did-known.R). The design is
# saturated in the four group-by-period cells, so given the cell sizes every
# study is one of prod(n_c + 1) outcomes with a known binomial probability,
# and the estimate and its standard error are closed forms in the four cell
# means:
#   index analysis:  p11 - G(G^-1(p10) + G^-1(p01) - G^-1(p00)), delta method
#   linear analysis: p11 - p10 - p01 + p00, HC1
# A study with a cell without events or without non-events is failed under
# the index analysis (separation = "fail"); a study whose outcomes are all
# equal, or with an empty cell, is failed under both. Cell sizes are
# multinomial: the exact conditional power is averaged over `draws` draws of
# them (the caller's RNG state is left alone), and `size_se` is the
# standard error that averaging leaves. `pc` and `share` are ordered
# (control before, control after, treated before, treated after).
exact_did_xs <- function(pc, n, link = "logit", analysis = c("index", "linear"),
                         share = rep(0.25, 4), conf = 0.95, draws = 300L, seed = 1L) {
  analysis <- match.arg(analysis)
  G <- if (link == "probit") pnorm else plogis
  Gi <- if (link == "probit") qnorm else qlogis
  g <- if (link == "probit") dnorm else dlogis
  z <- qnorm(1 - (1 - conf) / 2)
  ws <- getFromNamespace("with_seed", "powerape")
  sizes <- ws(seed, t(rmultinom(draws, n, share)))
  out <- fail <- numeric(draws)
  for (r in seq_len(draws)) {
    nc <- sizes[r, ]
    if (any(nc == 0)) { fail[r] <- 1; next }
    dims <- nc + 1L
    ## expand a vector over cell i to the full four-way array
    E <- function(x, i) aperm(array(x, dim = c(dims[i], dims[-i])), order(c(i, setdiff(1:4, i))))
    ph <- lapply(nc, function(m) (0:m) / m)
    W <- E(dbinom(0:nc[1], nc[1], pc[1]), 1) * E(dbinom(0:nc[2], nc[2], pc[2]), 2) *
      E(dbinom(0:nc[3], nc[3], pc[3]), 3) * E(dbinom(0:nc[4], nc[4], pc[4]), 4)
    P <- lapply(1:4, function(i) E(ph[[i]], i))
    if (analysis == "linear") {
      est <- P[[4]] - P[[3]] - P[[2]] + P[[1]]
      v <- lapply(1:4, function(i) E(ph[[i]] * (1 - ph[[i]]) / nc[i], i))
      se <- sqrt(n / (n - 4) * (v[[1]] + v[[2]] + v[[3]] + v[[4]]))
      allsame <- (P[[1]] == P[[2]]) & (P[[2]] == P[[3]]) & (P[[3]] == P[[4]]) &
        (P[[1]] == 0 | P[[1]] == 1)
      ok <- se > 0 & !allsame
    } else {
      inner <- lapply(ph, function(x) x > 0 & x < 1)
      ok <- E(inner[[1]], 1) & E(inner[[2]], 2) & E(inner[[3]], 3) & E(inner[[4]], 4)
      li <- lapply(seq_along(ph), function(i) ifelse(inner[[i]], Gi(pmin(pmax(ph[[i]], 1e-12), 1 - 1e-12)), NA))
      L <- E(li[[3]], 3) + E(li[[2]], 2) - E(li[[1]], 1)
      est <- P[[4]] - G(L)
      vl <- lapply(1:3, function(i) E(ph[[i]] * (1 - ph[[i]]) / (nc[i] * g(li[[i]])^2), i))
      se <- sqrt(E(ph[[4]] * (1 - ph[[4]]) / nc[4], 4) + g(L)^2 * (vl[[1]] + vl[[2]] + vl[[3]]))
    }
    det <- ok & (est - z * se > 0)
    det[is.na(det)] <- FALSE
    out[r] <- sum(W[det])
    fail[r] <- sum(W[!ok])
  }
  list(power = mean(out), size_se = sd(out) / sqrt(draws), failed = mean(fail),
       failed_se = sd(fail) / sqrt(draws))
}
