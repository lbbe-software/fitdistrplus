
mygof <- function(gof, echo=FALSE)
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
      if(echo)
        print(cbind(2 * 1:n - 1, s, theop, 1-rev(theop), log(theop) + log(1 - rev(theop)), (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ),
              digit=22)
      - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ) 
    }
  }else if (gof == "ADR")
  {
    fnobj <- function(par, fix.arg, obs, pdistnam)
    { 
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam,c(list(s),as.list(par),as.list(fix.arg)))
      if(echo)
        print(cbind(2 * 1:n - 1, s, 1-rev(theop), log(1 - rev(theop)), (2 * 1:n - 1) * log(1 - rev(theop)) ))
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
      if(echo)
        print(cbind(theop, log(theop), log(1-theop), log(theop) + log(1 - theop), ((2 * 1:n - 1) / (1 - rev(theop)))))
      2 * sum(log(theop) + log(1 - theop) ) + 
        mean ( ((2 * 1:n - 1) / theop) + ((2 * 1:n - 1) / (1 - rev(theop))) )
    }
  }else
    fnobj <- NULL
  
  fnobj
}

pdistnam <- pweibull
x <- c(rep(1, 13), 10)
par <- fitdistrplus:::startarg_transgamma_family(x, "weibull")
fix.arg <- NULL

eps <- sqrt(.Machine$double.eps)

mygof("AD", TRUE)(par, fix.arg, x, pdistnam)
mygof("ADR", TRUE)(par, fix.arg, x, pdistnam)
mygof("AD2", TRUE)(par, fix.arg, x, pdistnam)


optim(par, function(theta) mygof("AD2", FALSE)(theta, fix.arg, x, pdistnam), control=list(trace=1, REPORT=1))

abs(mygof("AD")(as.list(unlist(par)+eps), fix.arg, x, pdistnam) - mygof("AD")(par, fix.arg, x, pdistnam))/eps
