#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### bootstrap in fitdistrplus with censored data
###
###         R functions
### 

bootdistcens <- function (f, niter=1001)
{ 
    if (niter<10) 
        stop("niter must be an integer above 10")
    if (!inherits(f, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    # non parametric bootstrap
    n <- length(f$censdata[, 1])
    numrow <- seq(1, n)
    rnumrow <- sample(numrow, size=niter*n, replace=TRUE)
    dim(rnumrow) <- c(n, niter)
    start <- as.list(f$estimate) #a named vector is no longer is accepted as starting values.
    if(is.function(f$fix.arg.fun))
      fix.arg <- f$fix.arg.fun
    else 
      fix.arg <- f$fix.arg
    if (is.null(f$dots))
        funcmle <- function(iter) {
        mle <- do.call(mledist, list(data=data.frame(left=f$censdata[rnumrow[, iter], ]$left, 
            right=f$censdata[rnumrow[, iter], ]$right), distr=f$distname, start=start, 
            fix.arg=fix.arg))
        return(c(mle$estimate, mle$convergence))
        }
    else
        funcmle <- function(iter) {
        mle <- do.call(mledist, c(list(data=data.frame(left=f$censdata[rnumrow[, iter], ]$left, 
            right=f$censdata[rnumrow[, iter], ]$right), distr=f$distname, start=start), 
            fix.arg=fix.arg, f$dots))
        return(c(mle$estimate, mle$convergence))
        }
    resboot <- sapply(1:niter, funcmle)
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
    res <- structure(list(estim=estim, converg=t(resboot)[, length(resboot[, 1])], 
                          method="nonparam", nbboot=niter, CI=bootCI, fitpart=f), 
        class="bootdistcens")
    res
}

print.bootdistcens <- function(x, ...){
    if (!inherits(x, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    cat("Parameter values obtained with nonparametric bootstrap \n")
    print(x$estim, ...)    
    nconverg <- length(x$converg[x$converg==0])
    if (nconverg < length(x$converg))
    {
        cat("\n")
        cat("The estimation method converged only for", nconverg, "among", 
                length(x$converg), "iterations \n")
    }

}

plot.bootdistcens <- function(x, ...){
    if (!inherits(x, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    if (dim(x$estim)[2]==1) {
        stripchart(x$estim, method="jitter", 
        xlab="Bootstrapped values of the parameter", ...)
    }
    else {
        if (dim(x$estim)[2]==2)
            plot(x$estim, 
            main="Bootstrapped values of the two parameters", ...)
        else 
            plot(x$estim, 
            main="Bootstrapped values of parameters", ...)
    }
}

summary.bootdistcens <- function(object, ...){
    if (!inherits(object, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    
    class(object) <- c("summary.bootdistcens", class(object))  
    object
    
}


print.summary.bootdistcens <- function(x, ...){
    if (!inherits(x, "summary.bootdistcens"))
        stop("Use only with 'summary.bootdistcens' objects")
    cat("Nonparametric bootstrap medians and 95% percentile CI \n")
    print(x$CI)
    
    nconverg <- length(x$converg[x$converg==0])
    if (nconverg < length(x$converg))
    {
        cat("\n")
        cat("The estimation method converged only for", nconverg, "among", 
                length(x$converg), "iterations \n")
    }
}
