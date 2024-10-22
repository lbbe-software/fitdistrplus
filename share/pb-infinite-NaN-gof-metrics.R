
mygof <- function(gof, echo=FALSE)
{
  if (gof == "CvM")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
    }
  }else if (gof == "KS")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam) 
    {
      n <- length(obs)
      s <- sort(obs)
      obspu <- seq(1,n)/n
      obspl <- seq(0,n-1)/n
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      max(pmax(abs(theop-obspu), abs(theop-obspl)))
    }
  }else if (gof == "AD")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      if(echo)
        print(cbind(2 * 1:n - 1, s, theop, 1-rev(theop), log(theop) + log(1 - rev(theop)), (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ),
              digit=22)
      ilogpi <- log(theop * (1 - rev(theop))) * (2 * 1:n - 1)
      idx <- is.finite(ilogpi)
      - sum(idx) - mean( ilogpi[idx] ) 
    }
  }else if (gof == "ADR")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      if(echo)
        print(cbind(2 * 1:n - 1, s, 1-rev(theop), log(1 - rev(theop)), (2 * 1:n - 1) * log(1 - rev(theop)) ))
      
      ilogpi <- log(1 - rev(theop)) * (2 * 1:n - 1)
      idx <- is.finite(ilogpi)
      sum(idx)/2 - 2 * sum(theop[idx]) - mean ( ilogpi[idx] )
    }
  }else if (gof == "ADL")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      ilogpi <- (2 * 1:n - 1) * log(theop)
      idx <- is.finite(ilogpi)
      -3*sum(idx)/2 + 2 * sum(theop[idx]) - mean ( ilogpi[idx] )
    }
  }else if (gof == "AD2R")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      
      logpi <- log(1 - theop)
      i1pi2 <- (2 * 1:n - 1) / (1 - rev(theop))
      idx <- is.finite(logpi) & is.finite(i1pi2)
      
      2 * sum(logpi[idx]) + mean( i1pi2[idx] )
    }
  }else if (gof == "AD2L")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      logpi <- log(theop)
      i1pi <- (2 * 1:n - 1) / theop 
      idx <- is.finite(logpi) & is.finite(i1pi)
      2 * sum( logpi[idx] ) + mean ( i1pi[idx] )
    }
  }else if (gof == "AD2")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      if(echo)
        print(cbind(theop, log(theop), log(1-theop), log(theop) + log(1 - theop), ((2 * 1:n - 1) / (1 - rev(theop)))))
      
      logpi <- log(theop * (1 - theop)) 
      i1pi <- (2 * 1:n - 1) / theop
      i1pi2 <- (2 * 1:n - 1) / (1 - rev(theop))
      idx <- is.finite(logpi) & is.finite(i1pi) & is.finite(i1pi2)
      
      2 * sum( logpi[idx] ) + mean ( i1pi[idx]  + i1pi2[idx] )
    }
  }else
    fnobj <- NULL
  
  fnobj
}

pdistnam <- pweibull
x <- c(rep(1, 14), 10)
par <- fitdistrplus:::startarg_transgamma_family(x, "weibull")
fix.arg <- NULL

eps <- sqrt(.Machine$double.eps)

mygof("AD", TRUE)(par, fix.arg, x, pdistnam)
mygof("ADR", TRUE)(par, fix.arg, x, pdistnam)
mygof("AD2", TRUE)(par, fix.arg, x, pdistnam)

optim(par, function(theta) mygof("AD", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))

optim(par, function(theta) mygof("AD2", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))

optim(par, function(theta) mygof("ADR", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))

optim(par, function(theta) mygof("AD2R", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))

optim(par, function(theta) mygof("AD2L", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))



