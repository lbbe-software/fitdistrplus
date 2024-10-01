require("fitdistrplus")
require("actuar")

data(groundbeef)
obs <- groundbeef$serving

gammastart <- fitdistrplus:::startargdefault(obs, "gamma")
ll <- fitdistrplus:::manageparam(gammastart, list(rate = 0.5), obs, "gamma")

fitdist(obs, "gamma", fix.arg = list(rate = 0.5))
fitdist(obs, "gamma", fix.arg = list(scale = 2))


pareto2start <- fitdistrplus:::startargdefault(obs, "pareto2")
ll <- fitdistrplus:::manageparam(pareto2start, list(rate = 0.5), obs, "pareto2")

fitdist(obs, "pareto2", fix.arg = list(rate = 0.5), lower=0)
fitdist(obs, "pareto2", fix.arg = list(scale = 2), lower=0)
