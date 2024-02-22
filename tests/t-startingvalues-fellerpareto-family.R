library(fitdistrplus)
library(actuar)

set.seed(1234)

n <- 1e3
x <- rfpareto(n, 1, 2, 2, 2, scale=2)

fitdistrplus:::startarg_fellerpareto_family(x, "pareto4")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto3")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto2")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto1")

fitdistrplus:::startarg_fellerpareto_family(x, "pareto")

fitdistrplus:::startarg_fellerpareto_family(x, "llogis")

fitdist(x, "pareto4")
fitdist(x, "pareto3")
fitdist(x, "pareto2", lower=0)
fitdist(x, "pareto1", lower=0)
fitdist(x, "pareto", lower=0)

fitdist(x, "llogis", lower=0)


fitdistrplus:::startarg_fellerpareto_family(x, "fpareto")

fitdist(x, "fpareto", lower=0)

fitdistrplus:::startarg_fellerpareto_family(x-1, "trbeta")

fitdist(x-1, "trbeta", lower=0)


fitdistrplus:::startarg_fellerpareto_family(x-1, "genpareto")

fitdist(x-1, "genpareto", lower=0)

fitdistrplus:::startarg_fellerpareto_family(x-1, "paralogis")

fitdist(x-1, "paralogis", lower=0)


x <- rfpareto(n, 0, 1, 2, 2, scale=2)

fitdistrplus:::startarg_fellerpareto_family(x, "invburr")
fitdistrplus:::startarg_fellerpareto_family(x, "invpareto")
fitdistrplus:::startarg_fellerpareto_family(x, "invparalogis")

