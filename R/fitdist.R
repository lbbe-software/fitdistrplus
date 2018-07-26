#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis, Christophe Dutang                                                                                                  
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
### fit parametric distributions for non-censored data
###
###         R functions
### 

fitdist <- function (data, distr, method = c("mle", "mme", "qme", "mge"), start=NULL, 
                     fix.arg=NULL, discrete, keepdata = TRUE, keepdata.nb=100, ...) 
{
    #check argument distr
    if (!is.character(distr)) 
        distname <- substring(as.character(match.call()$distr), 2)
    else 
        distname <- distr
    ddistname <- paste("d", distname, sep="")
    if (!exists(ddistname, mode="function"))
        stop(paste("The ", ddistname, " function must be defined"))
    
    #pdistname <- paste("p", distname, sep="")
    #if (!exists(pdistname, mode="function"))
    #    stop(paste("The ", pdistname, " function must be defined"))
    #check argument discrete
    if(missing(discrete))
    {
      if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
        discrete <- TRUE
      else
        discrete <- FALSE
    }
    if(!is.logical(discrete))
      stop("wrong argument 'discrete'.")
    if(!is.logical(keepdata) || !is.numeric(keepdata.nb) || keepdata.nb < 2)
      stop("wrong arguments 'keepdata' and 'keepdata.nb'")
    #check argument method
    if(any(method == "mom"))
        warning("the name \"mom\" for matching moments is NO MORE used and is replaced by \"mme\"")
    
    method <- match.arg(method, c("mle", "mme", "qme", "mge"))
    if(method %in% c("mle", "mme", "mge"))
      dpq2test <- c("d", "p")
    else
      dpq2test <- c("d", "p", "q")
    #check argument data
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")
 
    #encapsulate three dots arguments
    my3dots <- list(...)    
    if (length(my3dots) == 0) 
      my3dots <- NULL
    n <- length(data)
    
    # manage starting/fixed values: may raise errors or return two named list
    arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=data, 
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
    resdpq <- testdpqfun(distname, dpq2test, start.arg=arg_startfix$start.arg, 
               fix.arg=arg_startfix$fix.arg, discrete=discrete)
    if(any(!resdpq$ok))
    {
      for(x in resdpq[!resdpq$ok, "txt"])
        warning(x)
    }
    
    
    # Fit with mledist, qmedist, mgedist or mmedist
    if (method == "mme")
    {
        if (!is.element(distname, c("norm", "lnorm", "pois", "exp", "gamma", 
                "nbinom", "geom", "beta", "unif", "logis")))
            if (!"order" %in% names(my3dots))
                stop("moment matching estimation needs an 'order' argument")   
        
        mme <- mmedist(data, distname, start=arg_startfix$start.arg, 
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)
                
        sd <- NULL
        correl <- varcovar <- NULL
        
        estimate <- mme$estimate
        loglik <- mme$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        convergence <- mme$convergence
        fix.arg <- mme$fix.arg
        weights <- mme$weights
    }else if (method == "mle")
    {
        mle <- mledist(data, distname, start=arg_startfix$start.arg, 
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)
        if (mle$convergence>0) 
           stop("the function mle failed to estimate the parameters, 
                with the error code ", mle$convergence, "\n") 
        estimate <- mle$estimate
        if(!is.null(mle$hessian)){
            #check for NA values and invertible Hessian
            if(all(!is.na(mle$hessian)) && qr(mle$hessian)$rank == NCOL(mle$hessian)){
                varcovar <- solve(mle$hessian)
                sd <- sqrt(diag(varcovar))
                correl <- cov2cor(varcovar)
            }else{
                varcovar <- NA
                sd <- NA
                correl <- NA                            
            }
        }else{
            varcovar <- NA
            sd <- NA
            correl <- NA            
        }
        loglik <- mle$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        convergence <- mle$convergence
        fix.arg <- mle$fix.arg
        weights <- mle$weights
    }else if (method == "qme")
    {
        if (!"probs" %in% names(my3dots))
            stop("quantile matching estimation needs an 'probs' argument") 
                
        qme <- qmedist(data, distname, start=arg_startfix$start.arg, 
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)

        estimate <- qme$estimate
        sd <- NULL
        loglik <- qme$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        correl <- varcovar <- NULL
        
        convergence <- qme$convergence   
        fix.arg <- qme$fix.arg
        weights <- qme$weights
    }else if (method == "mge")
    {
        if (!"gof" %in% names(my3dots))
            warning("maximum GOF estimation has a default 'gof' argument set to 'CvM'")    

        mge <- mgedist(data, distname, start=arg_startfix$start.arg, 
                       fix.arg=arg_startfix$fix.arg, checkstartfix=TRUE, ...)

        estimate <- mge$estimate
        sd <- NULL
        loglik <- mge$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        correl <- varcovar <- NULL
        
        convergence <- mge$convergence
        fix.arg <- mge$fix.arg
        weights <- NULL
    }else
    {
        stop("match.arg() does not work correctly")
    }
    
    #needed for bootstrap
    if (!is.null(fix.arg)) 
      fix.arg <- as.list(fix.arg)
    
    if(keepdata)
    {
      reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl, 
                  vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=data,
                  distname = distname, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                  dots = my3dots, convergence = convergence, discrete = discrete, 
                  weights = weights)
    }else #just keep a sample set of all observations
    {
      n2keep <- min(keepdata.nb, n)-2
      imin <- which.min(data)
      imax <- which.max(data)
      subdata <- data[sample((1:n)[-c(imin, imax)], size=n2keep, replace=FALSE)]
      subdata <- c(subdata, data[c(imin, imax)])
      
      reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl, 
                  vcov = varcovar, loglik = loglik, aic=aic, bic=bic, n = n, data=subdata,
                  distname = distname, fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, 
                  dots = my3dots, convergence = convergence, discrete = discrete, 
                  weights = weights)  
    }
    
    
    return(structure(reslist, class = "fitdist"))

}

print.fitdist <- function(x, ...)
{
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    if (x$method=="mme") 
        cat("Fitting of the distribution '", x$distname, "' by matching moments \n")
    else if (x$method=="mle") 
       cat("Fitting of the distribution '", x$distname, "' by maximum likelihood \n")
    else if (x$method=="qme") 
        cat("Fitting of the distribution '", x$distname, "' by matching quantiles \n")
    else if (x$method=="mge") 
        cat("Fitting of the distribution '", x$distname, "' by maximum goodness-of-fit \n")
    
    cat("Parameters:\n")
    if (x$method=="mle") 
        print(cbind.data.frame("estimate" = x$estimate, "Std. Error" = x$sd), ...)
     else 
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

plot.fitdist <- function(x, breaks="default", ...)
{
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    if(!is.null(x$weights))
      stop("The plot of the fit is not yet available when using weights")
  if(!is.null(x$data))
      plotdist(data=x$data, distr=x$distname, 
        para=c(as.list(x$estimate), as.list(x$fix.arg)), breaks=breaks, 
        discrete = x$discrete, ...)
    if(!is.null(x$weights))
      stop("The plot of the fit is not yet available when using weights")
}

summary.fitdist <- function(object, ...)
{
    if (!inherits(object, "fitdist"))
        stop("Use only with 'fitdist' objects")
    object$ddistname <- paste("d", object$distname, sep="")
    object$pdistname <- paste("p", object$distname, sep="")
    object$qdistname <- paste("q", object$distname, sep="")
    
    class(object) <- c("summary.fitdist", class(object))    
    object
}

print.summary.fitdist <- function(x, ...)
{
    if (!inherits(x, "summary.fitdist"))
        stop("Use only with 'summary.fitdist' objects")

    ddistname <- x$ddistname
    pdistname <- x$pdistname
    
    if (x$method=="mme") 
        cat("Fitting of the distribution '", x$distname, "' by matching moments \n")
    else if (x$method=="mle") 
       cat("Fitting of the distribution '", x$distname, "' by maximum likelihood \n")
    else if (x$method=="qme") 
        cat("Fitting of the distribution '", x$distname, "' by matching quantiles \n")
    else if (x$method=="mge") 
        cat("Fitting of the distribution '", x$distname, "' by maximum goodness-of-fit \n")
    
    cat("Parameters : \n")
    if (x$method == "mle")
        print(cbind.data.frame("estimate" = x$estimate, "Std. Error" = x$sd), ...)
    else
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

    cat("Loglikelihood: ", x$loglik, "  ")
    cat("AIC: ", x$aic, "  ")
    cat("BIC: ", x$bic, "\n")
    
    if (x$method=="mle") {
        if (length(x$estimate) > 1) {
            cat("Correlation matrix:\n")
            print(x$cor)
            cat("\n")
        }
    }
    
    invisible(x)
}

#see quantiles.R for quantile.fitdist
#see logLik.R for loglik.fitdist
#see vcov.R for vcov.fitdist
#see coef.R for coef.fitdist

