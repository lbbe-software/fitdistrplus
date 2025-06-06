require("fitdistrplus")

# (1) fit a user-defined continuous distribution (Gumbel) to censored data.
#

data(fluazinam)
log10EC50 <-log10(fluazinam)
# definition of the Gumbel distribution
dgumbel  <-  function(x,a,b) 1/b*exp((a-x)/b)*exp(-exp((a-x)/b))
pgumbel  <-  function(q,a,b) exp(-exp((a-q)/b))
qgumbel  <-  function(p,a,b) a-b*log(-log(p))

mledist(log10EC50, "gumbel", start=list(a=0,b=2), optim.method="Nelder-Mead")
f1 <- mledist(log10EC50, "gumbel", start=list(a=0,b=2)) #default NM
c(f1$estimate, f1$loglik)
#fitted coef is 1.63 1.14 , fitted loglik is -20.3

try(mledist(rbind(log10EC50, c("a", "b")), "gamma"))
try(mledist(rbind(log10EC50, c(NA, NA)), "gamma"))
try(mledist(rbind(log10EC50, c(3, Inf)), "gamma"))
try(mledist(rbind(log10EC50, c(3, NaN)), "gamma"))


# (2) test optimization arguments to censored data MLE.
#

mledist(log10EC50, "lnorm", optim.method="BFGS")
f1 <- mledist(log10EC50, "lnorm", optim.method="Nelder")
c(f1$estimate, f1$loglik)
#fitted coef is 0.623   0.928, fitted loglik is -21.1

#optim() is used
mledist(log10EC50, "lnorm", optim.method="L-BFGS-B", lower=0)
mledist(log10EC50, "lnorm", optim.method="BFGS", lower=0) #a simple warning

#constrOptim() is used
mledist(log10EC50, "lnorm", optim.method="Nelder", lower=0)
mledist(log10EC50, "lnorm", lower=0) #error

# (3) weighted MLE
#
xleft <- c(-1.8, -0.6, -0.1, 0.07, 0.14, 1, 1.2, 1.2, 1.2)
xright <- c(-1.8, -0.6, -0.1, 0.07, 0.14, 1, NA, NA, NA)
d <- data.frame(left = xleft, right = xright)
f1 <- mledist(d, "norm")
c(f1$estimate, f1$loglik)
#fitted coef is 0.545 1.342, fitted loglik is -12.9

dbis <- data.frame(left = c(-1.8, -0.6, -0.1, 0.07, 0.14, 1, 1.2), 
                   right = c(-1.8, -0.6, -0.1, 0.07, 0.14, 1, NA))
f2 <- mledist(dbis, "norm", weights = c(rep(1,6),3))

# f1 and f2 must give quite the same results (only starting values differ)
f1$estimate
f2$estimate

# (4) test the definition of fix.arg/start.arg as functions
#     if defined as functions, start.arg and fix.arg must be
#     functions of pseudo data (output pseudo of cens2pseudo())

mledist(d, "norm", start = function(x) list(mean = 0, sd = 1))$estimate
mledist(d, "norm", start = function(x) list(mean = mean(x), sd = 1))$estimate
mledist(d, "norm", fix.arg = function(x) list(mean = 0))$estimate
mledist(d, "norm", fix.arg = function(x) list(mean = 0.544))$estimate
mledist(d, "norm", fix.arg = function(x) list(mean = mean(x)))$estimate
