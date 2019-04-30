library(fitdistrplus)

set.seed(123)


#--------------------------------------------------------
# exponential sample

x1 <- rexp(1e3)

mseKL_exp <- msedist(x1, "exp")
mseJ_exp <- msedist(x1, "exp", phidiv = "J")
mseR2_exp <- msedist(x1, "exp", phidiv = "R", power.phidiv=2)
mseR1o2_exp <- msedist(x1, "exp", phidiv = "R", power.phidiv=1/2)
mseH3o2_exp <- msedist(x1, "exp", phidiv = "H", power.phidiv=3/2)
mseV3o2_exp <- msedist(x1, "exp", phidiv = "V", power.phidiv=3/2)

c(true=1, mseKL_exp$estimate, mseJ_exp$estimate, mseR2_exp$estimate, 
  mseR1o2_exp$estimate, mseH3o2_exp$estimate, mseV3o2_exp$estimate)



mseKL_exp <- fitdist(x1, "exp", method="mse", phidiv="KL")
mseJ_exp <- fitdist(x1, "exp", method="mse", phidiv = "J")
mseR2_exp <- fitdist(x1, "exp", method="mse", phidiv = "R", power.phidiv=2)
mseR1o2_exp <- fitdist(x1, "exp", method="mse", phidiv = "R", power.phidiv=1/2)
mseH3o2_exp <- fitdist(x1, "exp", method="mse", phidiv = "H", power.phidiv=3/2)
mseV3o2_exp <- fitdist(x1, "exp", method="mse", phidiv = "V", power.phidiv=3/2)



gofstat(list(mseKL_exp, mseJ_exp, mseR2_exp, mseR1o2_exp, mseH3o2_exp, mseV3o2_exp))


cdfcomp(list(mseKL_exp, mseJ_exp, mseR2_exp, mseR1o2_exp, mseH3o2_exp, mseV3o2_exp), do.points=FALSE,
        legendtext = c("Kullback-Leibler", "Jeffreys", "Renyi alpha=2", "Renyi alpha=1/2", "Hellinger p=3/2",
                       "Vajda beta=3/2"))

qqcomp(list(mseKL_exp, mseJ_exp, mseR2_exp, mseR1o2_exp, mseH3o2_exp, mseV3o2_exp), 
        legendtext = c("Kullback-Leibler", "Jeffreys", "Renyi alpha=2", "Renyi alpha=1/2", "Hellinger p=3/2",
                       "Vajda beta=3/2"))

denscomp(list(mseKL_exp, mseJ_exp, mseR2_exp, mseR1o2_exp, mseH3o2_exp, mseV3o2_exp), demp = TRUE,
       legendtext = c("Kullback-Leibler", "Jeffreys", "Renyi alpha=2", "Renyi alpha=1/2", "Hellinger p=3/2",
                      "Vajda beta=3/2"))


#defensive test 

try(msedist(x1, "exp", phidiv="ABC"))
try(msedist(x1, "exp", phidiv="K", power.phidiv = "a"))
try(msedist(x1, "exp", phidiv="J", power.phidiv = "a"))
try(msedist(x1, "exp", phidiv="R", power.phidiv = 0))
try(msedist(x1, "exp", phidiv="R", power.phidiv = 1:10))
try(msedist(x1, "exp", phidiv="H", power.phidiv = 0))
try(msedist(x1, "exp", phidiv="H", power.phidiv = 1:10))



#--------------------------------------------------------
# Poisson sample

x1 <- rpois(1e3, lambda=15)

#no weight
mseKL_pois1 <- fitdist(x1, "pois", method="mse", phidiv="KL")
mseJ_pois1 <- fitdist(x1, "pois", method="mse", phidiv = "J")
mseR2_pois1 <- fitdist(x1, "pois", method="mse", phidiv = "R", power.phidiv=2)
mseR1o2_pois1 <- fitdist(x1, "pois", method="mse", phidiv = "R", power.phidiv=1/2)
mseH3o2_pois1 <- fitdist(x1, "pois", method="mse", phidiv = "H", power.phidiv=3/2)
mseV3o2_pois1 <- fitdist(x1, "pois",method="mse",  phidiv = "V", power.phidiv=3/2)


#with weight
mseKL_pois2 <- fitdist(unique(sort(x1)), "pois", method="mse", phidiv="KL", weights=as.numeric(table(x1)))
mseJ_pois2 <- fitdist(unique(sort(x1)), "pois", method="mse", phidiv = "J", weights=as.numeric(table(x1)))
mseR2_pois2 <- fitdist(unique(sort(x1)), "pois", method="mse", phidiv = "R", power.phidiv=2, weights=as.numeric(table(x1)))
mseR1o2_pois2 <- fitdist(unique(sort(x1)), "pois", method="mse", phidiv = "R", power.phidiv=1/2, weights=as.numeric(table(x1)))
mseH3o2_pois2 <- fitdist(unique(sort(x1)), "pois", method="mse", phidiv = "H", power.phidiv=3/2, weights=as.numeric(table(x1)))
mseV3o2_pois2 <- fitdist(unique(sort(x1)), "pois",method="mse",  phidiv = "V", power.phidiv=3/2, weights=as.numeric(table(x1)))

#taking into account partially unbias the estimation of the true cdf
cdfcomp(list(mseKL_pois1, mseJ_pois1, mseR2_pois1), do.points=FALSE, fitlty=1,
        legendtext = c("Kullback-Leibler", "Jeffreys", "Renyi alpha=2"))

cdfcomp(list(mseKL_pois2, mseJ_pois2, mseR2_pois2), do.points=FALSE, add=TRUE, 
        fitlty=2, addlegend = FALSE, datacol="white")


cdfcomp(list(mseR1o2_pois1, mseH3o2_pois1, mseV3o2_pois1), do.points=FALSE, fitlty=1,
        legendtext = c("Renyi alpha=1/2", "Hellinger p=3/2", "Vajda beta=3/2"))

cdfcomp(list(mseR1o2_pois2, mseH3o2_pois2, mseV3o2_pois2), do.points=FALSE, add=TRUE, 
        fitlty=2, addlegend = FALSE, datacol="white")


