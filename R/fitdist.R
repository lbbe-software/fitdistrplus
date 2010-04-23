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

fitdist <- function (data, distr, method=c("mle", "mme"), start, ...) 
{
    if (!is.character(distr)) 
        distname <- substring(as.character(match.call()$distr),2)
    else 
        distname <- distr
    ddistname <- paste("d",distname,sep="")
    
    if (!exists(ddistname,mode="function"))
        stop(paste("The ",ddistname," function must be defined"))
    pdistname <- paste("p",distname,sep="")
    if (!exists(pdistname,mode="function"))
        stop(paste("The ",pdistname," function must be defined"))
        
    if(any(method == "mom"))
        warning("the name \"mom\" for matching moments is NO MORE used and is replaced by \"mme\".")
    
    method <- match.arg(method)

    
    if (!missing(start) & method=="mme")
        warnings("Starting values for parameters will not be used with matching moments")  
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")
    n <- length(data)
    # MLE fit with mledist or matching moments fit with mmedist
    if (method=="mme")
    {
        estimate <- mmedist(data, distname)
        sd <- NULL
        loglik <- NULL
        aic <- NULL
        bic <- NULL
        correl <- NULL
    }
    else
    {
        if (missing(start))
            mle <- mledist(data, distname, ...) 
        else 
            mle <- mledist(data, distname, start, ...)
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
    }     
    
    return(structure(list(estimate = estimate, method = method, sd = sd, 
    cor = correl, loglik = loglik, aic=aic, bic=bic,
    n = n, data=data, distname=distname), class = "fitdist"))
        
}

print.fitdist <- function(x, ...){
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    if (x$method=="mme") 
        cat("Fitting of the distribution '",x$distname,"' by matching moments \n")
    else
       cat("Fitting of the distribution '",x$distname,"' by maximum likelihood \n")
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
    para=as.list(x$estimate),breaks=breaks,...)
}

summary.fitdist <- function(object, ...){
    if (!inherits(object, "fitdist"))
        stop("Use only with 'fitdist' objects")
    object$ddistname <- paste("d", object$distname,sep="")
    object$pdistname <- paste("p", object$distname,sep="")
    
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
    else
       cat("Fitting of the distribution '",x$distname,"' by maximum likelihood \n")
    cat("Parameters : \n")
    
    if (x$method=="mle") {
        print(cbind.data.frame("estimate" = x$estimate, "Std. Error" = x$sd), ...)
        cat("Loglikelihood: ",x$loglik,"  ")
        cat("AIC: ",x$aic,"  ")
        cat("BIC: ",x$bic,"\n")
        if (length(x$estimate) > 1) {
            cat("Correlation matrix:\n")
            print(x$cor)
            cat("\n")
        }
    }
    else {
        print(cbind.data.frame("estimate" = x$estimate))
    }
    
    invisible(x)
}
