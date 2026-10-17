# Helpers for the website: the R | Stata tabs and the Stata output that the
# articles quote. Sourced by the articles and vignettes when pkgdown builds
# the site, and by precompute/precompile.R when it knits the home page.
# Nothing here is part of the package.

.pc_dir <- function() {
  d <- getOption("powerape.precompute.dir")
  if (!is.null(d)) return(d)
  # Knitted by pkgdown: walk up from the article's folder to the package root.
  root <- normalizePath(getwd(), winslash = "/")
  while (!file.exists(file.path(root, "DESCRIPTION"))) {
    up <- dirname(root)
    if (identical(up, root)) stop("precompute/helpers.R: package root not found")
    root <- up
  }
  file.path(root, "precompute")
}

# ---- Stata ------------------------------------------------------------------
# The do-files in precompute/stata/ mark every block an article shows:
#
#     * >>> name
#     ...commands...
#     * <<<
#
# stata_code()   returns the commands of a block, read from the do-file;
# stata_output() returns what Stata printed for it, read from the log that the
#                run wrote, without the echoed commands;
# stata_show()   writes both as fenced code blocks (chunk option results = "asis").

# The do-file is precompute/stata/<run>.do; the log it wrote is
# precompute/stata/output/<run>.txt.
.stata_file <- function(run, what = c("do", "log")) {
  what <- match.arg(what)
  f <- if (what == "do") file.path(.pc_dir(), "stata", paste0(run, ".do"))
       else file.path(.pc_dir(), "stata", "output", paste0(run, ".txt"))
  if (!file.exists(f)) stop("Missing ", f, ": run precompute/stata/", run, ".do first")
  readLines(f, warn = FALSE, encoding = "UTF-8")
}

.block <- function(lines, start, end, what) {
  i <- grep(start, lines)
  if (length(i) != 1) stop(sprintf("%s: found %d start markers, expected 1", what, length(i)))
  j <- grep(end, lines)
  j <- j[j > i]
  if (!length(j)) stop(sprintf("%s: no end marker", what))
  if (j[1] == i + 1) character(0) else lines[(i + 1):(j[1] - 1)]
}

.trim_blank <- function(x) {
  keep <- which(nzchar(trimws(x)))
  if (!length(keep)) return(character(0))
  x <- x[min(keep):max(keep)]
  # at most one blank line in a row
  blank <- !nzchar(trimws(x))
  x[!(blank & c(FALSE, head(blank, -1)))]
}

stata_code <- function(run, name) {
  x <- .stata_file(run, "do")
  .trim_blank(.block(x, sprintf("^\\* >>> %s\\s*$", name), "^\\* <<<\\s*$",
                     sprintf("%s.do, block '%s'", run, name)))
}

# Stata wraps lines longer than its linesize (120 in these runs) with "> " at
# the start of the continuation. After a command those lines belong to the
# command and are dropped with it; after output they are joined back to the
# line they continue. Lines still longer than `width` are then wrapped at
# word boundaries, keeping their indentation.
.unwrap_log <- function(b, width = 100) {
  out <- character(0)
  state <- "out"
  for (line in b) {
    if (grepl("^\\. ", line) || identical(line, ".")) {
      state <- "cmd"
    } else if (grepl("^>", line)) {
      if (state != "cmd" && length(out)) {
        out[length(out)] <- paste0(out[length(out)], substring(line, 3))
      }
    } else if (state == "cmd" && grepl("^\\s+[0-9]+\\.( |$)", line)) {
      # the echoed body of a loop or program definition
    } else {
      out <- c(out, line)
      state <- "out"
    }
  }
  out <- sub("\\s+$", "", out)
  long <- nchar(out) > width + 10
  if (any(long)) {
    out <- unlist(lapply(seq_along(out), function(i) {
      if (!long[i]) return(out[i])
      ind <- nchar(sub("^( *).*$", "\\1", out[i]))
      paste0(strrep(" ", ind),
             strwrap(trimws(out[i]), width = width - ind, exdent = 2))
    }))
  }
  out
}

stata_output <- function(run, name) {
  x <- .stata_file(run, "log")
  b <- .block(x, sprintf("^\\. \\* >>> %s\\s*$", name), "^\\. \\* <<<\\s*$",
              sprintf("%s.txt, block '%s'", run, name))
  .trim_blank(.unwrap_log(b))
}

# stata_show(run, names): the commands of one or more blocks in one code
# block, then what Stata printed. `output` is TRUE (all blocks), FALSE (none)
# or the names of the blocks whose output to show.
stata_show <- function(run, names, code = TRUE, output = TRUE) {
  if (code) {
    cmd <- unlist(lapply(names, function(n) c(stata_code(run, n), "")))
    cat("```stata", .trim_blank(cmd), "```", "", sep = "\n")
  }
  show <- if (isTRUE(output)) names else if (isFALSE(output)) character(0) else output
  out <- unlist(lapply(show, function(n) c(stata_output(run, n), "")))
  out <- .trim_blank(out)
  if (length(out)) cat("```{.stata-output}", out, "```", "", sep = "\n")
  invisible(NULL)
}

# Key results that a Stata run stored in <run>-results.txt (see _prelude.do),
# one row per label.
stata_results <- function(run) {
  f <- file.path(.pc_dir(), "stata", "output", paste0(run, "-results.txt"))
  x <- utils::read.delim(f, na.strings = ".", strip.white = TRUE,
                         stringsAsFactors = FALSE)
  rownames(x) <- x$label
  x
}

# Full path of a graph exported by a Stata run.
stata_graph <- function(file) {
  src <- file.path(.pc_dir(), "stata", "output", file)
  if (!file.exists(src)) stop("Missing ", src)
  src
}

# A graph exported by a Stata run, embedded in the page as a data URI so that
# the article does not depend on a separate image file (chunk option
# results = "asis").
stata_figure <- function(file, alt, caption = NULL) {
  cat(sprintf('<figure class="stata-figure"><img src="%s" alt="%s" />%s</figure>',
              knitr::image_uri(stata_graph(file)), alt,
              if (is.null(caption)) "" else sprintf("<figcaption>%s</figcaption>", caption)),
      "", sep = "\n")
  invisible(NULL)
}

# Seconds that the blocks of a run took (<run>-timings.txt), by label.
stata_seconds <- function(run) {
  f <- file.path(.pc_dir(), "stata", "output", paste0(run, "-timings.txt"))
  x <- utils::read.delim(f, strip.white = TRUE, stringsAsFactors = FALSE)
  stats::setNames(x$seconds, x$label)
}

# ---- Numbers in running text ---------------------------------------------------
# `lead = FALSE` drops the leading zero of quantities that cannot exceed 1.
num <- function(x, digits = 2, lead = TRUE) {
  s <- formatC(x, format = "f", digits = digits, big.mark = ",")
  minus <- intToUtf8(8722)                          # a proper minus sign
  s <- sub("^-", minus, s)
  if (!lead) s <- sub(paste0("^(", minus, "?)0\\."), "\\1.", s)
  s
}

# ---- R | Stata tabs -----------------------------------------------------------
# Used as inline R code, each call on a line of its own with a blank line
# before and after:
#
#     `r tabs_r("id")`
#     ...R chunk(s)...
#     `r tabs_stata("id")`
#     ...Stata block(s)...
#     `r tabs_end()`
#
# `id` must be unique on the page. The markup is plain Bootstrap 5; the script
# in pkgdown/extra.js keeps all tabs on a page on the same language.
.tab_button <- function(id, lang, label, active) sprintf(
  paste0('<li class="nav-item" role="presentation"><button class="nav-link%s" ',
         'id="%s-%s-tab" data-bs-toggle="tab" data-bs-target="#%s-%s" type="button" ',
         'role="tab" aria-controls="%s-%s" aria-selected="%s" data-lang="%s">%s</button></li>'),
  if (active) " active" else "", id, lang, id, lang, id, lang,
  if (active) "true" else "false", lang, label)

tabs_r <- function(id) paste(
  '<div class="lang-tabs">',
  '<ul class="nav nav-tabs" role="tablist">',
  .tab_button(id, "r", "R", TRUE),
  .tab_button(id, "stata", "Stata", FALSE),
  '</ul>',
  '<div class="tab-content">',
  sprintf('<div class="tab-pane active" id="%s-r" role="tabpanel" aria-labelledby="%s-r-tab">', id, id),
  sep = "\n")

tabs_stata <- function(id) paste(
  '</div>',
  sprintf('<div class="tab-pane" id="%s-stata" role="tabpanel" aria-labelledby="%s-stata-tab">', id, id),
  sep = "\n")

tabs_end <- function() paste('</div>', '</div>', '</div>', sep = "\n")
