# The 1,200 pilot rows of "Difference-in-differences designs"
# (vignettes/did-designs.Rmd), written to did-pilot.csv for the Stata run
# did.do: Stata cannot draw R's random numbers, so it reads the rows the
# vignette's R code draws (same seed, same order of draws; on the website the
# vignette checks that the file matches). All four columns are integers. Run
# from the package root:
#
#     Rscript precompute/stata/did-pilot.R

set.seed(1)
u <- rnorm(600)
treated <- rbinom(600, 1, 0.5)
pilot <- data.frame(
  firm = rep(1:600, 2), year = rep(1:2, each = 600), treated = rep(treated, 2),
  y = rbinom(1200, 1, pnorm(-0.7 + 0.3 * rep(treated, 2) + 0.8 * rep(u, 2))))

write.csv(pilot, file.path("precompute", "stata", "did-pilot.csv"),
          row.names = FALSE, quote = FALSE)
