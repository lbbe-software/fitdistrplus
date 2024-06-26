library(fitdistrplus)
nbboot <- 100
nbboot <- 10

nsample <- 100
nsample <- 10

visualize <- FALSE # TRUE for manual tests with visualization of results


#### sanity check -- data ####

try(fitdist(c(serving, "a"), "gamma"))
try(fitdist(c(serving, NA), "gamma"))
try(fitdist(c(serving, Inf), "gamma"))
try(fitdist(c(serving, -Inf), "gamma"))
try(fitdist(c(serving, NaN), "gamma"))

#### sanity check -- distr ####

try(fitdist(serving, "toto"))

#### sanity check -- method ####

try(fitdist(serving, "gamma", method="toto"))
try(fitdist(serving, "gamma", method=1))

#### sanity check -- start ####

try(fitdist(serving, "gamma", start=list("a"=1, b=2)))

#### sanity check -- fix.arg ####

try(fitdist(serving, "gamma", fix.arg=list("a"=1, b=2)))
try(fitdist(serving, "gamma", fix.arg=list("shape"=1, rate=2)))


#### sanity check -- discrete ####

try(fitdist(serving, "gamma", discrete=3))

#### sanity check -- keepdata ####

try(fitdist(serving, "gamma", keepdata=3))
try(fitdist(serving, "gamma", keepdata=TRUE, keepdata.nb = 1))


#### sanity check -- calcvcov ####

try(fitdist(serving, "gamma", calcvcov=3))


#### check the warning messages when using weights in the fit followed by functions ####
# that do not yet take weights into account
# with an example to be used later to see if weights are well taken into account
#
if(visualize)
{
  x3 <- rnorm(100) # this sample size must be fixed here (see next lines, 50+50)
  x3 <- sort(x3)
  (f <- fitdist(x3, "norm", method="mle", weights= c(rep(1, 50), rep(2, 50))))
  try(plot(f))
  try(cdfcomp(f))
  (f2 <- fitdist(x3, "logis", method="mle", weights= c(rep(1, 50), rep(2, 50))))
  try(cdfcomp(list(f,f2)))
  try(denscomp(f))
  try(denscomp(list(f,f2)))
  try(ppcomp(f))
  try(ppcomp(list(f,f2)))
  try(qqcomp(f))
  try(qqcomp(list(f,f2)))
  try(gofstat(f))
  try(gofstat(list(f,f2)))
  try(bootdist(f))
}