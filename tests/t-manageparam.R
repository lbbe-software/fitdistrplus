library(fitdistrplus)


manageparam <- fitdistrplus:::manageparam

obs1 <- rnorm(10)
s1 <- NULL
s2 <- list("mean"=2, "sd"=3)
s3 <- function(x)
  list("mean"=1.01*mean(x)) 
s4 <- list("mean"=1)
  
f1 <- NULL  
f2 <- list("sd"=3)
f3 <- function(x) 
  list("sd"=1.01*sd(x))  
f4 <- list("toto"=2)

#no error

manageparam(s1, f1, obs1, "norm")
manageparam(s2, f1, obs1, "norm")
manageparam(s3, f1, obs1, "norm")

manageparam(s1, f2, obs1, "norm")
manageparam(s1, f3, obs1, "norm")

#raise error

try(manageparam(matrix(3), f1, obs1, "norm"))
try(manageparam(function(x) c("a"=33), f1, obs1, "norm"))
try(manageparam(function(x) list(33), f1, obs1, "norm"))
try(manageparam(NULL, list(mean=1, sd=1), obs1, "norm"))


#no error

checkparamlist <- fitdistrplus:::checkparamlist

myformal <- names(formals("dnorm"))

res <- manageparam(s1, f1, obs1, "norm")
checkparamlist(res$start.arg, res$fix.arg, myformal)

res <- manageparam(s1, f2, obs1, "norm")
checkparamlist(res$start.arg, res$fix.arg, myformal)

res <- manageparam(s1, f3, obs1, "norm")
checkparamlist(res$start.arg, res$fix.arg, myformal)

res <- manageparam(s2, f1, obs1, "norm")
checkparamlist(res$start.arg, res$fix.arg, myformal)

#raise errors

res <- manageparam(s1, f4, obs1, "norm")
try(checkparamlist(res$start.arg, res$fix.arg, myformal))

res <- manageparam(s2, f2, obs1, "norm")
try(checkparamlist(res$start.arg, res$fix.arg, myformal))

res <- manageparam(s2, f3, obs1, "norm")
try(checkparamlist(res$start.arg, res$fix.arg, myformal))



#no error
fitdist(obs1, "norm", start=NULL, fix.arg=NULL)
fitdist(obs1, "norm", start=NULL, fix.arg=f3)


#raise error
try(fitdist(obs1, "norm", start=NULL, fix.arg=f4))
try(fitdist(obs1, "norm", start=s2, fix.arg=f2))
try(fitdist(obs1, "norm", start=s2, fix.arg=f3))
