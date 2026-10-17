# The website: Stata runs and the precomputed home page

The website (<https://jespernwulff.github.io/powerape/>) documents the R
package and the Stata command side by side. Nothing in this folder is part
of the package (`.Rbuildignore`).

| Source | What it is | Stata runs it quotes |
|---|---|---|
| `pkgdown/index.Rmd.orig` | the home page; knitted into `pkgdown/index.md` | `home` |
| `vignettes/articles/worked-example.Rmd.orig` | Worked example: five designs, R and Stata tabs; knitted into `worked-example.Rmd` | `worked-example-a` to `worked-example-e` |
| `vignettes/articles/stata.Rmd` | Get started in Stata | `getting-started` |
| `vignettes/articles/stata-help.Rmd` | the Stata help file as text | `stata-help` |
| `vignettes/articles/r-and-stata.Rmd`, `validation.Rmd`, `stata-news.Rmd` | written by hand | none |
| `vignettes/classic-sesoi-detection.Rmd` | workflow, R and Stata tabs | `classic-sesoi` |
| `vignettes/equivalence.Rmd` | workflow, R and Stata tabs | `equivalence` |
| `vignettes/interaction-effects.Rmd` | workflow, R and Stata tabs | `interaction` |
| `vignettes/panel-designs.Rmd` | workflow, R and Stata tabs | `panel` |
| `vignettes/iv-designs.Rmd` | workflow, R and Stata tabs | `iv` |
| `vignettes/justified-alpha.Rmd` | workflow, R and Stata tabs | `justified-alpha` |
| `vignettes/minimum-effect.Rmd` | Get started in R (R only) | none |

The six workflow vignettes are also the package's vignettes. Their Stata
tabs exist only on the website: their setup chunk sources `helpers.R` when
pkgdown builds them (`IN_PKGDOWN` is `"true"`) and defines empty stand-ins
otherwise, so `R CMD build` and `R CMD check` render the same R-only
vignettes as before.

## Rebuilding

1. **Stata.** Run each file from `precompute/stata/`, with
   `POWERAPE_STATA_DIR` pointing to the `ado/` folder of the version to
   document (a clone of <https://github.com/jespernwulff/powerape-stata>),
   or with that version installed:

   ```
   cd precompute/stata
   stata -e do getting-started.do
   stata -e do home.do
   ...
   stata -e do stata-help.do
   ```

   Everything an article quotes is written to `precompute/stata/output/`:
   the log (`<run>.txt`), key results (`<run>-results.txt`), timings and
   graphs. Stata's own batch log, `<run>.log`, is not kept. Most runs take a
   few minutes, `panel.do` included: version 0.1.3 of the Stata command
   calibrates a panel world in seconds (under 0.1.2, `panel.do` took about
   20 minutes). The worked example takes longest: `worked-example-e.do`
   about half an hour, most of it the instrumented sample-size search;
   `worked-example-d.do` about 20 minutes and `worked-example-a.do` about
   10, most of it the `probit` fits of every simulated study in the single
   power calls. These times are for four runs at a time on a four-core
   machine. **Do not edit a do-file while it runs**: Stata reads do-files as
   it goes.

   `worked-example-c.do` reads the 900 pilot rows that the worked example's
   R code draws from `design-c-pilot.csv` (Stata cannot draw R's random
   numbers). The file is under version control; `design-c-pilot.R` writes
   it again (`Rscript precompute/stata/design-c-pilot.R` from the package
   root), and the article stops if the two ever differ.

2. **Home page and worked example.** With the version of the R package to
   document installed, from the package root:

   ```
   Rscript precompute/precompile.R
   ```

   This knits `pkgdown/index.Rmd.orig` into `pkgdown/index.md` (and draws
   `man/figures/home-curves.png`) and
   `vignettes/articles/worked-example.Rmd.orig` into
   `vignettes/articles/worked-example.Rmd`; `precompile.R index` or
   `precompile.R worked-example` knits one of them. Long chunks are cached
   in `precompute/cache/` (not under version control); delete the cache
   after changing what the package prints. The worked example's R
   simulations take about an hour. Run with the environment variable
   `WE_SKIP_STATA=1`, the knit skips the Stata blocks, so the R cache can
   be filled while the Stata runs are still going; knit once more without
   it when they are done. The last chunk of the worked example checks the
   numbers its text is written around and stops the knit when a version of
   the package moves one of them.

3. **Site.** `pkgdown::build_site()` writes the site to `docs/`, which is
   what GitHub Pages serves. Build outside Dropbox (a synced folder locks
   pkgdown's temporary files) and copy `docs/` back.

## How Stata output gets into the articles

The run files mark every block that an article shows:

```stata
* >>> power
powerape power d, n(1500) claim(minimum) sesoi(0.05) nsim(500) seed(1)
* <<<
```

`stata_show("getting-started", "power")` in an article (see `helpers.R`)
prints the commands of that block, read from the do-file, and below them
what Stata printed, read from the log. No output is typed or pasted by hand,
and the numbers quoted in the running text are read from the same results.
The scheme follows the bootmakr site
(<https://jespernwulff.github.io/bootmakr/>).
