* Runs behind the home page: a support program that raises firm survival
* from 30% by ten points and is worth its cost from five points upward.
* Start from this folder (see _prelude.do).
do "_prelude.do" home

* >>> world
powerape dgp program, focal(program, binary p(0.5)) covariates((z, normal)) ///
    baseline(0.30) signal(0.10)
powerape setape program, target(0.10)
* <<<

tic
* >>> n-minimum
powerape n program, power(0.80) claim(minimum) sesoi(0.05)               ///
    nsim(1000) seed(502) engine(mata)
* <<<
keep_r n_minimum
toc n_minimum

tic
* >>> n-detect
powerape n program, power(0.80) claim(detect)                            ///
    nsim(1000) seed(503) engine(mata)
* <<<
keep_r n_detect
toc n_detect

file close res
file close tim
log close pc
exit, clear
