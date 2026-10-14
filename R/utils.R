# Internal helpers ------------------------------------------------------------

# Two-sided critical value for a conf-level CI.
zcrit <- function(conf) qnorm(1 - (1 - conf) / 2)

`%||%` <- function(a, b) if (is.null(a)) b else a

# Distance to the claim boundary, in standard errors, at which a claim
# reaches power `p` under the normal approximation: the inverse power map
# the searches of ape_n() and ape_mde() step with. A one-sided claim
# (detect, minimum; r = Inf) needs z + qnorm(p). The equivalence claim needs
# BOTH interval bounds inside the margin; with the truth `delta` SEs from
# the near bound and `r * delta` from the far one,
#   P = pnorm(delta - z) + pnorm(r * delta - z) - 1,
# so at a truth of 0 (r = 1) delta = z + qnorm((1 + p) / 2). Inverting the
# one-sided map instead put the first equivalence proposal at about 72% of
# the requirement at a truth of 0, and the search oscillated (through
# 1.11.1). The bracket is safe: f <= 0 at max(0, z + qnorm(p)) and
# f >= 0 at z + qnorm((1 + p) / 2) whenever r >= 1.
req_delta <- function(p, z, r = Inf) {
  if (!is.finite(r)) return(z + qnorm(p))
  hi <- z + qnorm((1 + p) / 2)
  f <- function(d) pnorm(d - z) + pnorm(r * d - z) - 1 - p
  if (f(hi) <= 1e-12) return(hi)          # r = 1: the bound is the root
  uniroot(f, c(max(0, z + qnorm(p)), hi), tol = 1e-12)$root
}

# Equivalence margin at which a design whose estimate has standard error
# `s` reaches power `p` when the truth sits at |t|: the inverse in m of the
# TOST power pnorm((m - t) / s - z) + pnorm((m + t) / s - z) - 1 (the
# margin search of ape_mde(); the standard error does not move with m).
req_margin <- function(p, s, t, z) {
  f <- function(m) pnorm((m - t) / s - z) + pnorm((m + t) / s - z) - 1 - p
  hi <- t + s * (z + qnorm((1 + p) / 2))
  if (f(hi) <= 1e-12) return(hi)          # t = 0: the bound is the root
  uniroot(f, c(t + s * max(0, z + qnorm(p)), hi), tol = 1e-12)$root
}

# Evaluate `code` under `seed` without disturbing the caller's RNG state.
# seed = NULL means: use the ambient RNG stream as-is.
with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old <- if (had) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had) {
      assign(".Random.seed", old, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  })
  set.seed(seed)
  force(code)
}
