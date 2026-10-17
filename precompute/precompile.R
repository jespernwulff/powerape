# Precompute the home page of the website (and any other *.Rmd.orig source).
#
# pkgdown builds its home page from a Markdown file, so the R output on it
# cannot be computed when the site is built. pkgdown/index.Rmd.orig is
# knitted once, here, into pkgdown/index.md, which already contains the code,
# its output and the Stata output; pkgdown then only formats that file. The
# worked example (vignettes/articles/worked-example.Rmd.orig) is knitted here
# too, into worked-example.Rmd, because its R simulations take about an hour.
# The other articles and vignettes need no precompiling: pkgdown knits them,
# and their Stata output is read from precompute/stata/output/ at that time.
#
# Usage, from the package root (the installed powerape is the one used):
#
#     Rscript precompute/precompile.R              # every *.Rmd.orig
#     Rscript precompute/precompile.R index        # one or more, by name
#
# What it needs is listed in precompute/README.md.

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(".", winslash = "/")
if (!file.exists(file.path(root, "DESCRIPTION")))
  stop("Run this script from the package root.")
options(powerape.precompute.dir = file.path(root, "precompute"))

sources <- list.files(file.path(root, "vignettes"), pattern = "\\.Rmd\\.orig$",
                      recursive = TRUE, full.names = TRUE)
names(sources) <- sub("\\.Rmd\\.orig$", "", basename(sources))
# The home page of the site: pkgdown/index.Rmd.orig -> pkgdown/index.md
home <- file.path(root, "pkgdown", "index.Rmd.orig")
if (file.exists(home)) sources <- c(sources, index = home)
if (length(args)) {
  unknown <- setdiff(args, names(sources))
  if (length(unknown)) stop("No such source: ", paste(unknown, collapse = ", "))
  sources <- sources[args]
}

knit_one <- function(name, orig) {
  dir <- dirname(orig)
  out <- file.path(dir, paste0(name, if (name == "index") ".md" else ".Rmd"))
  message("== ", name)

  knitr::opts_knit$restore()
  knitr::opts_chunk$restore()
  knitr::opts_knit$set(base.dir = paste0(dir, "/"), upload.fun = knitr::image_uri)
  knitr::opts_chunk$set(
    collapse   = TRUE,
    comment    = "#>",
    fig.path   = paste0("figures/", name, "-"),
    cache.path = paste0(file.path(root, "precompute", "cache", name), "/"),
    fig.width  = 7,
    fig.height = 4.5,
    dpi        = 144,
    dev        = "png"
  )

  knitr::knit(orig, out, envir = new.env(parent = globalenv()), quiet = TRUE)

  # Say where the file comes from: right after the YAML header, or at the top
  # if there is none.
  x   <- readLines(out, warn = FALSE, encoding = "UTF-8")
  end <- which(x == "---")[2]
  note <- sprintf("<!-- Generated from %s.Rmd.orig by precompute/precompile.R. Edit the .orig file, not this one. -->", name)
  x <- if (is.na(end)) c(note, "", x) else append(x, c("", note), after = end)
  writeLines(enc2utf8(x), out, useBytes = TRUE)
  message("   wrote ", sub(paste0(root, "/"), "", out, fixed = TRUE))
}

for (n in names(sources)) knit_one(n, sources[[n]])
