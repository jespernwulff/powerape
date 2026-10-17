* Shared set-up for the Stata runs behind the website articles. Not shown on
* the site. Each run file calls it as
*     do "_prelude.do" <name of the run>
* and is started from this folder (precompute/stata).
*
* Environment variable:
*   POWERAPE_STATA_DIR  folder holding powerape.ado and lpowerape.mlib (the
*                       ado/ folder of a checkout of
*                       https://github.com/jespernwulff/powerape-stata), if the
*                       version to document is not the installed one
args run

clear all
set more off
set linesize 120

* Everything a run produces goes to output/: the log the articles quote from
* (<run>.txt), key results, timings and graphs. Stata's own batch log,
* <run>.log in this folder, starts with the licence banner and is not kept.
capture mkdir "output"
global pc_out "`c(pwd)'/output"

local adodir : environment POWERAPE_STATA_DIR
if "`adodir'" != "" {
    adopath ++ "`adodir'"
}
mata: mata mlib index
capture program drop powerape

* Key results of each powerape call go to <run>-results.txt so that the
* article text can quote them without retyping. A result that a command does
* not store is written as a missing value.
capture file close res
file open res using "$pc_out/`run'-results.txt", write text replace
file write res "label" _tab "power" _tab "mcse" _tab "n" _tab "mde" _tab "alpha" ///
    _tab "nmax" _tab "nsim" _n

capture program drop keep_r
program define keep_r
    * usage: keep_r label   (right after a powerape command)
    args label
    file write res "`label'"
    foreach s in power mcse n mde alpha nmax nsim {
        local v = r(`s')
        file write res _tab %15.0g (`v')
    }
    file write res _n
end

* A robustness sweep stores its scenarios, not one power: keep_range writes
* the lowest and the highest scenario power as rows <label>_lo and <label>_hi
* (power column; the other columns missing).
capture program drop keep_range
program define keep_range
    * usage: keep_range label   (right after powerape robust)
    args label
    tempname S
    matrix `S' = r(scenarios)
    local j = colnumb(`S', "power")
    local lo = .
    local hi = .
    forvalues i = 1/`=rowsof(`S')' {
        local lo = min(`lo', `S'[`i', `j'])
        local hi = max(`hi', `S'[`i', `j'])
    }
    foreach k in lo hi {
        file write res "`label'_`k'" _tab %15.0g (``k'')
        forvalues c = 1/6 {
            file write res _tab "."
        }
        file write res _n
    }
end

* Timing: tic ... toc label  (writes seconds to <run>-timings.txt). toc runs
* -timer list-, which replaces r(): call keep_r first.
capture file close tim
file open tim using "$pc_out/`run'-timings.txt", write text replace
file write tim "label" _tab "seconds" _n
capture program drop tic
program define tic
    timer clear 99
    timer on 99
end
capture program drop toc
program define toc
    args label
    timer off 99
    quietly timer list 99
    file write tim "`label'" _tab %9.1f (r(t99)) _n
end

* Graphs in the site's colours.
capture program drop pc_graph
program define pc_graph
    * usage: pc_graph name   (exports the current graph as output/<name>.svg)
    args name
    quietly graph export "$pc_out/`name'.svg", replace
end

* The log that the articles quote from. It is opened last: after Stata's
* start-up banner (no licence details) and after the set-up above (no local
* paths).
capture log close pc
log using "$pc_out/`run'.txt", text replace nomsg name(pc)
