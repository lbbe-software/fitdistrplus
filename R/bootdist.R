#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis                                                                                                  
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
### bootstrap in fitdistrplus
###
###         R functions
### 


bootdist <- function (f, bootmethod="param", niter=1001)
{ 
    if (niter<10) 
        stop("niter must be an integer above 10")
        bootmethod <- match.arg(bootmethod, c("param", "nonparam"))
    
    if (!inherits(f, "fitdist"))
        stop("Use only with 'fitdist' objects")
        
    #simulate bootstrap data
    if (bootmethod == "param") { # parametric bootstrap
        rdistname <- paste("r", f$distname, sep="")
        if (!exists(rdistname, mode="function"))
            stop(paste("The ", rdistname, " function must be defined"))
        rdata <- do.call(rdistname, c(list(n=niter*f$n), as.list(f$estimate), f$fix.arg))
        dim(rdata) <- c(f$n, niter)
    }
    else { # non parametric bootstrap
        rdata <- sample(f$data, size=niter*f$n, replace=TRUE)
        dim(rdata) <- c(f$n, niter)
    }
    
    #compute bootstrap estimates
    foncestim <- switch(f$method, "mle"=mledist, "qme"=qmedist, "mme"=mmedist, "mge"=mgedist)
    start <- f$estimate
    if (is.null(f$dots))
        func <- function(iter) {
            res <- do.call(foncestim, list(data=rdata[, iter], distr=f$distname, start=start, fix.arg=f$fix.arg))
            return(c(res$estimate, res$convergence))
        }
    else
        func <- function(iter) {
            res <- do.call(foncestim, c(list(data=rdata[, iter], distr=f$distname, start=start, fix.arg=f$fix.arg), f$dots))
            return(c(res$estimate, res$convergence))
        }
    resboot <- sapply(1:niter, func)
    rownames(resboot) <- c(names(start), "convergence")
    if (length(resboot[, 1])>2) {
        estim <- data.frame(t(resboot)[, -length(resboot[, 1])])
        bootCI <- cbind(apply(resboot[-length(resboot[, 1]), ], 1, median, na.rm=TRUE), 
        apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.025, na.rm=TRUE), 
        apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.975, na.rm=TRUE))
        colnames(bootCI) <- c("Median", "2.5%", "97.5%")
    }
    else {
        estim <- as.data.frame(t(resboot)[, -length(resboot[, 1])])
        names(estim) <- names(f$estimate)
        bootCI <- c(median(resboot[-length(resboot[, 1]), ], na.rm=TRUE), 
        quantile(resboot[-length(resboot[, 1]), ], 0.025, na.rm=TRUE), 
        quantile(resboot[-length(resboot[, 1]), ], 0.975, na.rm=TRUE)) 
        names(bootCI) <- c("Median", "2.5%", "97.5%") 
    } 
    
    # code of convergence of the optimization function for each iteration
    converg <- t(resboot)[, length(resboot[, 1])]

    res <- structure(list(estim=estim, converg=converg, 
                          method=bootmethod, nbboot=niter, CI=bootCI, fitpart=f), 
                     class="bootdist")
    res    
}

print.bootdist <- function(x, ...){
    if (!inherits(x, "bootdist"))
        stop("Use only with 'bootdist' objects")
    if (x$method=="param") 
        cat("Parameter values obtained with parametric bootstrap \n")
    else
       cat("Parameter values obtained with nonparametric bootstrap \n")
    print(x$estim, ...)    
    nconverg <- length(x$converg[x$converg==0])
    if (nconverg < length(x$converg))
    {
        cat("\n")
        cat("The estimation method converged only for ", nconverg, " among ", 
                length(x$converg), " iterations \n")
    }

}

plot.bootdist <- function(x, ...){
    if (!inherits(x, "bootdist"))
        stop("Use only with 'bootdist' objects")
    if (dim(x$estim)[2]==1) {
        stripchart(x$estim, method="jitter", 
            xlab="Bootstrapped values of the parameter", ...)
    }
    else {
        if (dim(x$estim)[2]==2)
            plot(x$estim, 
            main="Bootstrapped values of parameters", ...)
        else 
            plot(x$estim, 
            main="Bootstrapped values of parameters", ...)
    }
}

summary.bootdist <- function(object, ...){
    if (!inherits(object, "bootdist"))
        stop("Use only with 'bootdist' objects")
    
    class(object) <- c("summary.bootdist", class(object))  
    object
}

print.summary.bootdist <- function(x, ...){
    
    if (!inherits(x, "summary.bootdist"))
        stop("Use only with 'summary.bootdist' objects")

    if (x$method=="param") 
        cat("Parametric bootstrap medians and 95% percentile CI \n")
    else
        cat("Nonparametric bootstrap medians and 95% percentile CI \n")
    print(x$CI)
    
    nconverg <- length(x$converg[x$converg==0])
    if (nconverg < length(x$converg))
    {
        cat("\n")
        cat("The estimation method converged only for ", nconverg, " among ", 
            length(x$converg), " iterations \n")
    }
}
