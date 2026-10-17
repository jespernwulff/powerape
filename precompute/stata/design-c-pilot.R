# The 900 pilot rows of Design C in the worked example
# (vignettes/articles/worked-example.Rmd.orig), written to design-c-pilot.csv
# for the Stata run worked-example-c.do: Stata cannot draw R's random
# numbers, so it reads the rows the article's R code draws (same seed, same
# order of draws; the article checks that the file matches). Values carry 17
# significant digits, so Stata's -import delimited, asdouble- reads back the
# doubles R holds. Run from the package root:
#
#     Rscript precompute/stata/design-c-pilot.R

set.seed(31)
z        <- rnorm(900)
program  <- rbinom(900, 1, 0.5)
survived <- rbinom(900, 1, pnorm(-0.55 + 0.35 * program + 0.45 * z))

out <- data.frame(survived = survived, program = program,
                  z = sprintf("%.17g", z))
write.csv(out, file.path("precompute", "stata", "design-c-pilot.csv"),
          row.names = FALSE, quote = FALSE)
