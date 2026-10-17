* Runs behind the Stata tabs of Design B in "Worked example: five designs":
* does the program work better for family firms? The same eight-point
* effect as a main effect and as an average interaction effect (AIE), with
* no covariates. Both searches use the Mata engine. Start from this folder
* (see _prelude.do).
do "_prelude.do" worked-example-b

* >>> world-ape
powerape dgp b, focal(program, binary p(0.5)) baseline(0.30)
powerape setape b, target(0.08)
* <<<

* >>> world-aie
powerape dgp bi, focal(program, binary p(0.5))                             ///
    moderator(family, binary p(0.5)) baseline(0.30)
powerape setaie bi, target(0.08) mainfocal(0.10) mainmoderator(0.03)
* <<<

tic
* >>> n-ape
powerape n b, power(0.80) claim(detect) nsim(1000) seed(521) engine(mata)
* <<<
keep_r n_ape
toc n_ape

tic
* >>> n-aie
powerape n bi, power(0.80) claim(detect) nsim(1000) seed(522) engine(mata)
* <<<
keep_r n_aie
toc n_aie

file close res
file close tim
log close pc
exit, clear
