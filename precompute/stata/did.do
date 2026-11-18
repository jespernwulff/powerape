* Runs behind the Stata tabs of "Difference-in-differences designs". The
* single power calls use the default engine, which fits every simulated study
* with -logit- and -vce(cluster)- on the unit identifier; the sample-size
* search and the sweep use the Mata engine. The pilot rows are the ones the
* vignette's R code draws (did-pilot.csv, written by did-pilot.R). Start from
* this folder (see _prelude.do).
do "_prelude.do" did

* >>> design
powerape dgp d, did model(logit) baseline(0.20) gap(0.20) trend(0.10)      ///
    design(panel) rho(0.5)
* <<<

* >>> pin
powerape setatt d, target(0.067)
* <<<

tic
* >>> power
powerape power d, n(3000) claim(detect) nsim(300) seed(1)
* <<<
keep_r power
toc power

tic
* >>> n
powerape n d, power(0.80) claim(detect) nsim(1000) seed(2) engine(mata)
* <<<
keep_r n
toc n

* >>> mismatch-world
powerape dgp d_mis, did model(logit) baseline(0.20) gap(0.20) trend(0.10)  ///
    rho(0.5) trends(linear) analysis(index)
powerape setatt d_mis, target(0.10)
* <<<

tic
* >>> mismatch
powerape power d_mis, n(1000) claim(detect) nsim(300) seed(3)
* <<<
keep_r mismatch
toc mismatch

tic
* >>> robust
powerape robust d, n(3000) claim(detect) nsim(200) seed(4)                 ///
    vary(rho(0.3 0.7) trends(index linear)) gridpoints(2) nonmax engine(mata)
* <<<
keep_range robust
toc robust

* >>> periods-world
powerape dgp d4, did model(logit) baseline(0.20) gap(0.20) trend(0.10)     ///
    rho(0.5) npre(2) npost(2)
powerape setatt d4, target(0.067)
* <<<

tic
* >>> periods
powerape power d4, n(1500) claim(detect) nsim(300) seed(5)
* <<<
keep_r periods
toc periods

* >>> subgroup
powerape dgp ds, did model(logit) baseline(0.30 0.20) gap(0.10 0.15)       ///
    trend(0.10 0.05) rho(0.5) subgroup(female, p(0.4))
* <<<

* >>> subgroup-pin
powerape setatt ds, target(0.10 0.03)
* <<<

* >>> pilot
* The pilot rows are the ones the R code draws (Stata cannot draw R's random
* numbers): precompute/stata/did-pilot.csv in the R package's repository.
import delimited using "did-pilot.csv", encoding(utf-8) clear
powerape dgp dp, didpilot outcome(y) group(treated) id(firm) period(year)  ///
    trend(0.05) model(logit)
* <<<

* >>> pilot-pin
powerape setatt dp, target(0.08)
* <<<

tic
* >>> statement-power
powerape power d, n(3000) claim(minimum) sesoi(0.02) nsim(300) seed(6)
* <<<
keep_r statement
toc statement

* >>> statement
powerape statement
* <<<

file close res
file close tim
log close pc
exit, clear
