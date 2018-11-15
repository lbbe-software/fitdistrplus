# INPUTS 
# argdistname : argument names of the distribution from names(formals())

# OUTPUTS 
# parameter names (as a vector) of the distribution (excluding non parameter argument)

computegetparam <- function(argdistname)
{
  #remove first argument, that should be "x", "p", "q", or "n", see ?dgamma, pgamma, qgamma
  argdistname <- argdistname[-1]
  nonparaminR <- c("x", "p", "q", "n") #defensive programming
  #remove other arguments, see ?dgamma, pgamma, qgamma, dbeta
  nonparaminR <- c(nonparaminR, "log", "log.p", "lower.tail", "ncp")
  nonparaminActuar <- c("limit", "order", "t")
  nonparaminGamlssdist <- "fast"
  nonparamspecial <- c("...", "..1", "..2")
  #see ?dnig, dhyperb, dskewlap, dgig,...
  nonparaminGenHyperbolic <- c("param", "KOmega", "ibfTol", "nmax", "method", "intTol",
                               "valueOnly", "nInterpol", "uniTol", "subdivisions", "logPars")
  #see ?dsn
  nonparamsn <- "dp"
  
  plist <- setdiff(argdistname, nonparaminR)
  plist <- setdiff(plist, nonparaminActuar)
  plist <- setdiff(plist, nonparaminGamlssdist)
  plist <- setdiff(plist, nonparamspecial)
  plist <- setdiff(plist, nonparaminGenHyperbolic)
  plist <- setdiff(plist, nonparamsn)
  
  plist
}

