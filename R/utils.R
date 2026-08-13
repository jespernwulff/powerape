# Internal helpers ------------------------------------------------------------

# Two-sided critical value for a conf-level CI.
zcrit <- function(conf) qnorm(1 - (1 - conf) / 2)

`%||%` <- function(a, b) if (is.null(a)) b else a

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
