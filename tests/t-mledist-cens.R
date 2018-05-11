library(fitdistrplus)




# (1) fit a continuous distribution (Gumbel) to censored data.
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
# definition of the Gumbel distribution
dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))

mledist(log10EC50, "gumbel", start=list(a=0,b=2), optim.method="Nelder-Mead")
mledist(log10EC50, "gumbel", start=list(a=0,b=2)) #default NM



# (2) test optimization arguments to censored data MLE.
#

mledist(log10EC50, "lnorm", optim.method="BFGS")
mledist(log10EC50, "lnorm", optim.method="Nelder")

#optim() is used
mledist(log10EC50, "lnorm", optim.method="L-BFGS-B", lower=0)
mledist(log10EC50, "lnorm", optim.method="BFGS", lower=0) #a simple warning

#constrOptim() is used
mledist(log10EC50, "lnorm", optim.method="Nelder", lower=0)
mledist(log10EC50, "lnorm", lower=0) #error
