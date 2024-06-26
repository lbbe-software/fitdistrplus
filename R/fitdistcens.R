#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### fit parametric distributions for censored data
###
###         R functions
### 

fitdistcens <- function (censdata, distr, start=NULL, fix.arg=NULL, 
                         keepdata = TRUE, keepdata.nb=100, calcvcov = TRUE, ...) 
{
  if (missing(censdata) || !is.data.frame(censdata) ||
      !(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[, 1])>1))
    stop("censdata must be a dataframe with two columns named left 
            and right and more than one line")
  
  checkCensoredDataFrameNAInfNan(censdata)
  
  if (!is.character(distr)) 
    distname <- substring(as.character(match.call()$distr), 2)
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  pdistname <- paste("p", distname, sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", pdistname, " function must be defined"))
  if(!is.logical(keepdata) || !is.numeric(keepdata.nb) || keepdata.nb < 3)
    stop("wrong arguments 'keepdata' and 'keepdata.nb'.")
  if(!is.logical(calcvcov) || length(calcvcov) != 1)
    stop("wrong argument 'calcvcov'.")
  #encapsulate three dots arguments
  my3dots <- list(...)    
  if (length(my3dots) == 0) 
    my3dots <- NULL
  
  #format data for calculation of starting values
  pseudodata <- cens2pseudo(censdata)$pseudo
  
  # manage starting/fixed values: may raise errors or return two named list
  arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=pseudodata, 
                              distname=distname)
  
  #check inconsistent parameters
  argddistname <- names(formals(ddistname))
  hasnodefaultval <- sapply(formals(ddistname), is.name)
  arg_startfix <- checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg, 
                                 argddistname, hasnodefaultval)
  #arg_startfix contains two names list (no longer NULL nor function)
  
  #store fix.arg.fun if supplied by the user
  if(is.function(fix.arg))
    fix.arg.fun <- fix.arg
  else
    fix.arg.fun <- NULL
  
  # check d, p, q, functions of distname
  dpq2test <- c("d", "p")
  resdpq <- testdpqfun(distname, dpq2test, start.arg=arg_startfix$start.arg, 
                       fix.arg=arg_startfix$fix.arg, discrete=FALSE)
  if(any(!resdpq$ok))
  {
    for(x in resdpq[!resdpq$ok, "txt"])
      warning(x)
  }
  
  
  # MLE fit with mledist 
  mle <- mledist(censdata, distname, start=arg_startfix$start.arg, 
                 fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE,
                 calcvcov=calcvcov, ...)
  if (mle$convergence>0) 
    stop("the function mle failed to estimate the parameters, 
        with the error code ", mle$convergence) 
  estimate <- mle$estimate
  varcovar <- mle$vcov
  if(!is.null(varcovar))
  {
    correl <- cov2cor(varcovar)
    sd <- sqrt(diag(varcovar))
  }else
    correl <- sd <- NULL
  
  loglik <- mle$loglik
  n <- nrow(censdata)
  npar <- length(estimate)
  aic <- -2*loglik+2*npar
  bic <- -2*loglik+log(n)*npar
  
  fix.arg <- mle$fix.arg
  weights <- mle$weights
  
  #needed for bootstrap
  if (!is.null(fix.arg)) 
    fix.arg <- as.list(fix.arg)
  
  if(keepdata)
  {
    reslist <- list(estimate = estimate, method="mle", sd = sd, cor = correl, 
                    vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n=n, censdata=censdata, 
                    distname=distname, fix.arg=fix.arg, fix.arg.fun = fix.arg.fun, 
                    dots=my3dots, convergence=mle$convergence, discrete=FALSE,
                    weights = weights)
  }else
  {
    n2keep <- min(keepdata.nb, n)-4
    imin <- unique(apply(censdata, 2, which.min))
    imax <- unique(apply(censdata, 2, which.max))
    subdata <- censdata[sample((1:n)[-c(imin, imax)], size=n2keep, replace=FALSE), ]
    subdata <- rbind.data.frame(subdata, censdata[c(imin, imax), ])
    
    reslist <- list(estimate = estimate, method="mle", sd = sd, cor = correl, 
                    vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n=n, censdata=subdata, 
                    distname=distname, fix.arg=fix.arg, fix.arg.fun = fix.arg.fun, 
                    dots=my3dots, convergence=mle$convergence, discrete=FALSE,
                    weights = weights)
  }
  
  return(structure(reslist, class = "fitdistcens"))
  
}

print.fitdistcens <- function(x, ...){
  if (!inherits(x, "fitdistcens"))
    stop("Use only with 'fitdistcens' objects")
  cat("Fitting of the distribution '", x$distname, "' on censored data by maximum likelihood \n")
  cat("Parameters:\n")
  print(cbind.data.frame("estimate" = x$estimate), ...)
  
  if(!is.null(x$fix.arg))
  {
    if(is.null(x$fix.arg.fun))
    {
      cat("Fixed parameters:\n")
    }else
    {
      cat("Fixed parameters (computed by a user-supplied function):\n")
    }
    print(cbind.data.frame("value" = unlist(x$fix.arg)), ...)
  }
  
}

plot.fitdistcens <- function(x, ...){
  if (!inherits(x, "fitdistcens"))
    stop("Use only with 'fitdistcens' objects")
  if(!is.null(x$weights))
    stop("The plot of the fit is not yet available when using weights")
  plotdistcens(censdata=x$censdata, distr=x$distname, 
               para=c(as.list(x$estimate), as.list(x$fix.arg)), ...)
  
}

summary.fitdistcens <- function(object, ...){
  if (!inherits(object, "fitdistcens"))
    stop("Use only with 'fitdistcens' objects")
  object$ddistname <- paste("d", object$distname, sep="")
  object$pdistname <- paste("p", object$distname, sep="")
  
  class(object) <- c("summary.fitdistcens", class(object))    
  object
}

print.summary.fitdistcens <- function(x, ...){
  if (!inherits(x, "summary.fitdistcens"))
    stop("Use only with 'fitdistcens' objects")
  ddistname <- x$ddistname
  pdistname <- x$pdistname
  
  cat("Fitting of the distribution '", x$distname, 
      "' By maximum likelihood on censored data \n")
  cat("Parameters\n")
  print(cbind.data.frame("estimate" = x$estimate, "Std. Error" = x$sd), ...)
  
  if(!is.null(x$fix.arg))
  {
    if(is.null(x$fix.arg.fun))
    {
      cat("Fixed parameters:\n")
    }else
    {
      cat("Fixed parameters (computed by a user-supplied function):\n")
    }
    print(cbind.data.frame("value" = unlist(x$fix.arg)), ...)
  }
  
  cat("Loglikelihood: ", x$loglik, "  ")
  cat("AIC: ", x$aic, "  ")
  cat("BIC: ", x$bic, "\n")
  if (length(x$estimate) > 1) {
    cat("Correlation matrix:\n")
    print(x$cor)
    cat("\n")
  }
  invisible(x)
}

#see quantiles.R for quantile.fitdistcens
#see logLik.R for loglik.fitdistcens
#see vcov.R for vcov.fitdistcens
#see coef.R for coef.fitdistcens
#see AIC.R for AIC.fitdistcens
#see BIC.R for BIC.fitdistcens

