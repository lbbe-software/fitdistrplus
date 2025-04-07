require("fitdistrplus")
require("actuar")

set.seed(1234)

n <- 1e3
x <- rfpareto(n, 1, 2, 2, 2, scale=2)

fitdistrplus:::startarg_fellerpareto_family(x, "pareto4")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto3")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto2")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto1")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto")

fitdistrplus:::startarg_fellerpareto_family(x, "llogis")

eps <- 1e-3
fitdist(x, "pareto4", control=list(trace=1), lower=eps, optim="L-BFGS-B")
fitdist(x, "pareto3", lower=eps, optim="L-BFGS-B")
#fitdist(x, "pareto2", lower=eps, optim="L-BFGS-B")
#fitdist(x, "pareto1", lower=eps, optim="L-BFGS-B")
#fitdist(x, "pareto", lower=eps, optim="L-BFGS-B")

fitdist(x, "llogis", lower=eps, optim="L-BFGS-B")


fitdistrplus:::startarg_fellerpareto_family(x, "fpareto")

fitdist(x, "fpareto", lower=eps, optim="L-BFGS-B")

fitdistrplus:::startarg_fellerpareto_family(x-1, "trbeta")

fitdist(x-1, "trbeta", lower=eps, optim="L-BFGS-B")


fitdistrplus:::startarg_fellerpareto_family(x-1, "genpareto")

fitdist(x-1, "genpareto", lower=eps, optim="L-BFGS-B")

fitdistrplus:::startarg_fellerpareto_family(x-1, "paralogis")

fitdist(x-1, "paralogis", lower=eps, optim="L-BFGS-B")


x <- rfpareto(n, 0, 1, 2, 2, scale=2)

fitdistrplus:::startarg_fellerpareto_family(x, "invburr")
fitdistrplus:::startarg_fellerpareto_family(x, "invpareto")
fitdistrplus:::startarg_fellerpareto_family(x, "invparalogis")

