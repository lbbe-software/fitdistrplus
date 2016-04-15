# exemples à intégrer dans l'aide de llplot

# 1 parametre
x <- rexp(50)
fite <- fitdist(x, "exp")
llplot(fite)
llplot(fite, col = "red")

# 2 parametres
data(groundbeef)
serving <- groundbeef$serving
fitg <- fitdist(serving, "gamma")
llplot(fitg)
llplot(fitg, expansion = 2)
llplot(fitg, pal.col = heat.colors(100))
llplot(fitg, back.col = FALSE, nlev = 25)

# 2 parametres dont un fixé
fitg2 <- fitdist(serving, "gamma", fix.arg = list(rate = 0.5))
llplot(fitg2)

# 3 parametres
data(endosulfan)
ATV <-endosulfan$ATV
library("actuar")
fBurr <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1))
llplot(fBurr)
llplot(fBurr, back.col = FALSE)
llplot(fBurr, nlev = 0, pal.col = rainbow(100), lseq = 100)

# 3 parametres dont un fixé
fBurr2 <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1),
                  fix.arg = list(rate = 1.5))
llplot(fBurr2)

fBurr3 <- fitdist(ATV, "burr", start = list(shape1 = 0.3, rate = 1),
                  fix.arg = list(shape2 = 1.5))
llsurface(data = ATV, distr = "burr", plot.arg = c("shape1", "rate"),
          min.arg = c(0.1, 1), max.arg = c(0.3, 2), fix.arg = list(shape2 = 1.5))
llplot(fBurr3)
