
require(fitdistrplus)
require(actuar)

data("danishuni")
try(fitdist(danishuni$Loss, "burr")) #Error code 1
try(fitdist(danishuni$Loss, "burr", control=list(maxit=1000))) #Error system is computationally singular: reciprocal condition number = 5.84972e-19
#Error non-finite finite-difference value
try(fitdist(danishuni$Loss, "burr", lower=0)) #Error code 7
try(fitBurr_cvg1 <- fitdist(danishuni$Loss, "burr", upper=100, control=list(trace=1))) #converged, loglik = -3373.199
try(fitdist(danishuni$Loss, "burr", upper=1000)) #Error code 7
try(fitBurr_cvg2 <- fitdist(danishuni$Loss, "burr", lower=.Machine$double.eps,
                            optim.method="L-BFGS-B", control=list(trace=1))) #converged, loglik = -3369.521
