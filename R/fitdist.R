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

fitdist <- function (data, distr, method=c("mle", "mme", "qme"), start=NULL, fix.arg=NULL, ...) 
{
    if (!is.character(distr)) 
        distname <- substring(as.character(match.call()$distr),2)
    else 
        distname <- distr
    ddistname <- paste("d",distname,sep="")
    if (!exists(ddistname, mode="function"))
        stop(paste("The ",ddistname," function must be defined"))
    
    pdistname <- paste("p",distname,sep="")
    if (!exists(pdistname, mode="function"))
        stop(paste("The ",pdistname," function must be defined"))
        
    if(any(method == "mom"))
        warning("the name \"mom\" for matching moments is NO MORE used and is replaced by \"mme\".")
    
    method <- match.arg(method, c("mle", "mme", "qme"))
    dots <- list(...)
    if (length(dots)==0) 
        dots <- NULL
    
#   if (!is.null(start) & method=="mme")
#        warnings("Starting values for parameters will not be used with matching moments")  
#    if (!is.null(fix.arg) & method=="mme")
#       stop("Matching moments cannot be used when some distribution parameters are fixed")  
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")
        
    
    n <- length(data)
    # Fit with mledist, qmedist or mmedist
    if (method == "mme")
    {
        mme <- mmedist(data, distname, start=start, fix.arg=fix.arg, ...)
                
        sd <- NULL
        correl <- NULL
        
        estimate <- mme$estimate
        loglik <- mme$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        convergence <- mme$convergence
    }else if (method == "mle")
    {
        mle <- mledist(data, distname, start, fix.arg, ...)
        if (mle$convergence>0) 
           stop("the function mle failed to estimate the parameters, 
                with the error code ",mle$convergence, "\n") 
        estimate <- mle$estimate
        if(!is.null(mle$hessian)){
            if(all(!is.na(mle$hessian))){
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
    }else if (method == "qme")
    {
        qme <- qmedist(data, distname, start, fix.arg, ...)

        estimate <- qme$estimate
        sd <- NULL
        loglik <- qme$loglik
        npar <- length(estimate)
        aic <- -2*loglik+2*npar
        bic <- -2*loglik+log(n)*npar
        correl <- NULL
        
        convergence <- qme$convergence      
    }else
    {
        stop("match.arg does not work correctly.")
    }
    
    if (!is.null(fix.arg)) fix.arg <- as.list(fix.arg)
    
    reslist <- list(estimate = estimate, method = method, sd = sd, cor = correl, 
                    loglik = loglik, aic=aic, bic=bic, n = n, data=data, 
                    distname = distname, fix.arg = fix.arg, dots = dots, 
                    convergence = convergence)

    
    if (method == "qme")
        reslist <- c(reslist, list(probs=qme$probs))
    if (method == "mme")
        reslist <- c(reslist, list(order=mme$order, memp=mme$memp))

    
    return(structure(reslist, class = "fitdist"))

}

print.fitdist <- function(x, ...){
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    if (x$method=="mme") 
        cat("Fitting of the distribution '",x$distname,"' by matching moments \n")
    else if (x$method=="mle") 
       cat("Fitting of the distribution '",x$distname,"' by maximum likelihood \n")
    else if (x$method=="qme") 
        cat("Fitting of the distribution '",x$distname,"' by matching quantiles \n")
    
    cat("Parameters:\n")
    if (x$method=="mle") 
        print(cbind.data.frame("estimate" = x$estimate, "Std. Error" = x$sd), ...)
     else 
        print(cbind.data.frame("estimate" = x$estimate))

}

plot.fitdist <- function(x, breaks="default", ...){
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    plotdist(data=x$data,distr=x$distname,
    para=c(as.list(x$estimate),as.list(x$fix.arg)),breaks=breaks,...)
}

summary.fitdist <- function(object, ...){
    if (!inherits(object, "fitdist"))
        stop("Use only with 'fitdist' objects")
    object$ddistname <- paste("d", object$distname,sep="")
    object$pdistname <- paste("p", object$distname,sep="")
    object$qdistname <- paste("q", object$distname,sep="")
    
    class(object) <- c("summary.fitdist", class(object))    
    object
}

print.summary.fitdist <- function(x, ...){
    if (!inherits(x, "summary.fitdist"))
        stop("Use only with 'fitdist' objects")

    ddistname <- x$ddistname
    pdistname <- x$pdistname
    
    if (x$method=="mme") 
        cat("Fitting of the distribution '",x$distname,"' by matching moments \n")
    else if (x$method=="mle") 
       cat("Fitting of the distribution '",x$distname,"' by maximum likelihood \n")
    else if (x$method=="qme") 
        cat("Fitting of the distribution '",x$distname,"' by matching quantiles \n")
    
    cat("Parameters : \n")
    if (x$method == "mle")
        print(cbind.data.frame("estimate" = x$estimate, "Std. Error" = x$sd), ...)
    else
        print(cbind.data.frame("estimate" = x$estimate), ...)

    cat("Loglikelihood: ",x$loglik,"  ")
    cat("AIC: ",x$aic,"  ")
    cat("BIC: ",x$bic,"\n")
    
    if (x$method=="mle") {
        if (length(x$estimate) > 1) {
            cat("Correlation matrix:\n")
            print(x$cor)
            cat("\n")
        }
    }
    
    invisible(x)
}
