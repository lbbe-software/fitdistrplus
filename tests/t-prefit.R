require(fitdistrplus)

# (1) gamma
x <- rgamma(1e3, 5/2, 7/2)

prefit(x, "gamma", "mle", list(shape=3, scale=3), lower=-Inf, upper=Inf, silent=TRUE, control=list(trace=1, REPORT=1))
prefit(x, "gamma", "mle", list(shape=1, scale=1), lower=-Inf, upper=Inf, silent=TRUE)

prefit(x, "gamma", "mle", list(shape=3), fix.arg=list(scale=7/2), lower=-Inf, upper=Inf, silent=TRUE)

prefit(x, "gamma", "qme", list(shape=1, scale=1), probs=1:2/3, lower=-Inf, upper=Inf, silent=TRUE)

# (2) geometric
x <- rgeom(1e3, 1/7)
prefit(x, "geom", "mle", list(prob=1/2), lower=-Inf, upper=Inf, silent=TRUE)
tbx <- table(x)
prefit(as.numeric(names(tbx)), "geom", "mle", list(prob=1/2), lower=-Inf, upper=Inf, silent=TRUE, weights=tbx)
prefit(x, "geom", "qme", list(prob=1/2), probs=1/2, lower=-Inf, upper=Inf)

# (3) Pareto
require(actuar)
x  <-  rpareto(1000, 6, 2)

prefit(x, "pareto", "mme", list(shape=10, scale=10), order=1:2, memp=function(x, order) mean(x^order), lower=-Inf, upper=Inf)

