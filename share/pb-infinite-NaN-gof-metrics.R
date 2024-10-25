

mygofold <- function(gof, echo=FALSE)
{
  if (gof == "CvM")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
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
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      max(pmax(abs(theop-obspu),abs(theop-obspl)))
    }
  }else if (gof == "AD")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ) 
    }
  }else if (gof == "ADR")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      n/2 - 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(1 - rev(theop)) )
    }
  }else if (gof == "ADL")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      -3*n/2 + 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(theop) )
    }
  }else if (gof == "AD2R")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      2 * sum(log(1 - theop)) + mean ( (2 * 1:n - 1) / (1 - rev(theop)) )
    }
  }else if (gof == "AD2L")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      2 * sum(log(theop)) + mean ( (2 * 1:n - 1) / theop )
    }
  }else if (gof == "AD2")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      2 * sum(log(theop) + log(1 - theop) ) + 
        mean ( ((2 * 1:n - 1) / theop) + ((2 * 1:n - 1) / (1 - rev(theop))) )
    }
  }else
    fnobj <- NULL
  
  fnobj
}

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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
      2 * sum( logpi[idx] ) + mean ( i1pi[idx]  + i1pi2[idx] )
    }
  }else
    fnobj <- NULL
  
  fnobj
}


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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
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
      print(sum(idx))
      2 * sum( logpi[idx] ) + mean ( i1pi[idx]  + i1pi2[idx] )
    }
  }else
    fnobj <- NULL
  
  fnobj
}

pdistnam <- pweibull
obs <- c(rep(1, 14), 10)
par <- fitdistrplus:::startarg_transgamma_family(x, "weibull")
fix.arg <- NULL

eps <- sqrt(.Machine$double.eps)

mygof("AD", TRUE)(par, fix.arg, x, pdistnam)
mygof("ADR", TRUE)(par, fix.arg, x, pdistnam)
mygof("AD2", TRUE)(par, fix.arg, x, pdistnam)


#### AD ####

optim(par, function(theta) mygof("AD", FALSE)(theta, fix.arg, obs, pdistnam), control=list(trace=1, REPORT=1))
optim(par, function(theta) mygofold("AD", FALSE)(theta, fix.arg, obs, pdistnam), control=list(trace=1, REPORT=1))

mygof1 <- function(shape)
{
  sapply(shape, function(y) mygof("AD", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofold1 <- function(shape)
{
  sapply(shape, function(y) mygofold("AD", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}

curve(mygof1(x), from=.25, to =2, main="gof=AD")
curve(mygofold1(x), add=TRUE, col="blue")

#### AD2 ####

optim(par, function(theta) mygof("AD2", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))



mygof1 <- function(shape)
{
  sapply(shape, function(y) mygof("AD2", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofold1 <- function(shape)
{
  sapply(shape, function(y) mygofold("AD2", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}

mygof1(1.7)
mygofold(1.7)

curve(mygof1(x), from=.25, to =1.6, main="gof=AD2", log="y")
curve(mygofold1(x), add=TRUE, col="blue")


#### ADR ####

optim(par, function(theta) mygof("ADR", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))


mygof1 <- function(shape)
{
  sapply(shape, function(y) mygof("ADR", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofold1 <- function(shape)
{
  sapply(shape, function(y) mygofold("ADR", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}

mygof1(2)
mygofold(2)

curve(mygof1(x), from=.25, to =2, main="gof=ADR")
curve(mygofold1(x), add=TRUE, col="blue")


#### AD2R ####

optim(par, function(theta) mygof("AD2R", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))


mygof1 <- function(shape)
{
  sapply(shape, function(y) mygof("AD2R", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofold1 <- function(shape)
{
  sapply(shape, function(y) mygofold("AD2R", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}

mygof1(2)
mygofold(2)

curve(mygof1(x), from=.25, to =2, main="gof=AD2R", log="y")
curve(mygofold1(x), add=TRUE, col="blue")


#### AD2L ####


optim(par, function(theta) mygof("AD2L", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))

mygof1 <- function(shape)
{
  sapply(shape, function(y) mygof("AD2L", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofold1 <- function(shape)
{
  sapply(shape, function(y) mygofold("AD2L", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}

mygof1(2)
mygofold(2)

curve(mygof1(x), from=.25, to =2, main="gof=AD2L")
curve(mygofold1(x), add=TRUE, col="blue")


#### test large sample size ####


mygofscale <- function(gof, echo=FALSE)
{
  if (gof == "CvM")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
      1/(12*n^2) + mean( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
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
      print(sum(idx))
      - sum(idx)/n - mean( ilogpi[idx] ) /n
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
      print(sum(idx))
      sum(idx)/2/n - 2 * sum(theop[idx]) /n - mean ( ilogpi[idx] ) /n
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
      print(sum(idx))
      -3*sum(idx)/2/n + 2 * sum(theop[idx]) /n - mean ( ilogpi[idx] ) /n
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
      print(sum(idx))
      2 * sum(logpi[idx])/n + mean( i1pi2[idx] )/n
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
      print(sum(idx))
      2 * sum( logpi[idx] )/n + mean ( i1pi[idx] )/n
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
      print(sum(idx))
      2 * sum( logpi[idx] )/n + mean ( i1pi[idx]  + i1pi2[idx] )/n
    }
  }else
    fnobj <- NULL
  
  fnobj
}


obs <- rlnorm(1e6, 3, 2)

mygof1 <- function(shape)
{
  sapply(shape, function(y) mygof("AD", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofold1 <- function(shape)
{
  sapply(shape, function(y) mygofold("AD", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}
mygofscale1 <- function(shape)
{
  sapply(shape, function(y) mygofscale("AD", FALSE)(c(y, 1), fix.arg, obs, pdistnam))
}

curve(mygof1(x), from=.25, to =2, main="gof=AD")
curve(mygofold1(x), add=TRUE, col="blue")

curve(mygofscale1(x), from=.25, to =2, main="gof=AD")

fitdistrplus::mgedist(obs, "lnorm", gof="AD")

