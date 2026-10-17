* Runs behind the Stata tabs of Design C in "Worked example: five designs":
* last year's 900-row pilot informs next year's design. Stata reads the
* pilot rows that the article's R code draws (design-c-pilot.csv, written by
* design-c-pilot.R), fits the same probit, and builds the world from the
* fit. Both searches use the Mata engine. Start from this folder (see
* _prelude.do).
do "_prelude.do" worked-example-c

* >>> pilot
import delimited using "design-c-pilot.csv", asdouble encoding(utf-8) clear
probit survived program z, nolog
* <<<

* >>> world
powerape dgp c, fromfit focal(program)
powerape setape c, target(0.08)
* <<<

* >>> world-25
powerape dgp c25, fromfit focal(program) baseline(0.25)
powerape setape c25, target(0.08)
* <<<

tic
* >>> n
powerape n c, power(0.80) claim(minimum) sesoi(0.03) nsim(1000) seed(511)  ///
    engine(mata)
* <<<
keep_r n
toc n

tic
* >>> n-25
powerape n c25, power(0.80) claim(minimum) sesoi(0.03) nsim(1000)          ///
    seed(512) engine(mata)
* <<<
keep_r n_25
toc n_25

file close res
file close tim
log close pc
exit, clear
