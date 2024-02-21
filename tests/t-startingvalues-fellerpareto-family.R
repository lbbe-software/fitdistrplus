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
