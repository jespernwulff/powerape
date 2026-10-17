* Runs behind the Stata tabs of Design A in "Worked example: five designs":
* a support program expected to raise firm survival from 30% by ten points,
* worth its cost from five points upward. The world and the two searches
* are those of the home page (home.do), with the same seeds. Single power
* calls use the default engine (probit fits every study), searches and the
* sweep the Mata engine. Start from this folder (see _prelude.do).
do "_prelude.do" worked-example-a

* >>> world
powerape dgp a, focal(program, binary p(0.5)) covariates((z, normal))      ///
    baseline(0.30) signal(0.10)
powerape setape a, target(0.10)
* <<<

tic
* >>> power
powerape power a, n(1500) claim(minimum) sesoi(0.05) nsim(4000) seed(501)
* <<<
keep_r power
toc power

tic
* >>> n-minimum
powerape n a, power(0.80) claim(minimum) sesoi(0.05) nsim(1000) seed(502)  ///
    engine(mata)
local n_min = r(n)
* <<<
keep_r n_minimum
toc n_minimum

tic
* >>> n-detect
powerape n a, power(0.80) claim(detect) nsim(1000) seed(503) engine(mata)
* <<<
keep_r n_detect
toc n_detect

tic
* >>> robust
powerape robust a, n(`n_min') claim(minimum) sesoi(0.05)                   ///
    vary(baseline(0.20 0.40) signal(0 0.30)) gridpoints(3)                 ///
    nsim(600) seed(504) nonmax engine(mata)
* <<<
keep_range robust
toc robust

tic
* >>> statement-power
powerape power a, n(`n_min') claim(minimum) sesoi(0.05) nsim(4000)         ///
    seed(20260821)
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
