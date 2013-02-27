#############################################################################
#   Copyright (c) 2012 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### calculation of theoretical quantiles from a parametric distribution
### fitted on censored or non-censored data
### 
###         R functions
### 

#quantile function for fitdist objects
quantile.fitdist <- function(x, probs = seq(0.1, 0.9, by=0.1), ...)
{
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    myquantiles.fitdist(f = x, probs = probs, cens = FALSE)
}

#quantile function for fitdistcens objects
quantile.fitdistcens <- function(x, probs = seq(0.1, 0.9, by=0.1), ...)
{
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    myquantiles.fitdist(f = x, probs = probs, cens = TRUE)
}

#internal quantile function for fitdist
myquantiles.fitdist <- function(f, probs, cens)
{    
    qdistname<-paste("q", f$distname, sep="")
    if (!exists(qdistname, mode="function"))
        stop(paste("The ", qdistname, " function must be defined")) 
     
    # computation and print of quantiles using estimations of parameters   
    para=c(as.list(f$estimate), as.list(f$fix.arg))
    quantiles <- do.call(qdistname, c(list(p=probs), as.list(para)))
    if (length(probs)>1)
        quantiles <- as.data.frame(t(quantiles))
    else
        quantiles <- as.data.frame(quantiles)
    colnames(quantiles) <- paste("p=", probs, sep="")
    rownames(quantiles) <- "estimate"
    
    reslist <- list(quantiles = quantiles, probs = probs)
    if(!cens)
        class(reslist) <- "quantile.fitdist"    
    else
        class(reslist) <- "quantile.fitdistcens"    

    reslist
}

print.quantile.fitdist <- function(x, ...)
{
    if (!inherits(x, "quantile.fitdist"))
        stop("Use only with 'quantile.fitdist' objects")
    typedata <- "(non-censored data)"
    
    cat("Estimated quantiles for each specified probability ", typedata,"\n", sep="")
    print(x$quantiles)  
    invisible(x)    
}

print.quantile.fitdistcens <- function(x, ...)
{
    if (!inherits(x, "quantile.fitdistcens"))
    stop("Use only with 'quantile.fitdistcens' objects")
    typedata <- "(censored data)"
    
    cat("Estimated quantiles for each specified probability ", typedata,"\n", sep="")
    print(x$quantiles)  
    invisible(x)    
}




#############################################################################
### calculation of theoretical quantiles from a parametric distribution
### fitted on censored or non-censored data
### and associated bootstrap confidence intervals
###
###         R functions
### 


#quantile function for bootdist objects
quantile.bootdist <- function(x, probs = seq(0.1, 0.9, by=0.1), 
    CI.type = "two.sided", CI.level = 0.95, ...)
{
    if (!inherits(x, "bootdist"))
        stop("Use only with 'bootdist' objects")
    myquantiles.bootdist(b = x, probs = probs, CI.type = CI.type, 
                         CI.level = CI.level, cens = FALSE)
}

#quantile function for bootdistcens objects
quantile.bootdistcens <- function(x, probs = seq(0.1, 0.9, by=0.1), 
    CI.type = "two.sided", CI.level = 0.95, ...)
{
    if (!inherits(x, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    myquantiles.bootdist(b = x, probs = probs, CI.type = CI.type, 
                         CI.level = CI.level, cens = TRUE)
}

#internal quantile function for bootdist
myquantiles.bootdist <- function(b, probs, CI.type, CI.level, cens)
{    
    CI.type <- match.arg(CI.type, c("two.sided", "less", "greater"))
    if(!is.logical(cens))
        stop("wrong argument cens.")
    CI.level <- CI.level[1]
        
    # 1/ computation of quantiles using quantile.fitdist
    basequant <- quantile(b$fitpart, probs=probs)
        
    # 2/ computation of bootstraped quantiles and alpha-percent CI of quantiles     
    qdistname <- paste("q", b$fitpart$distname, sep="")
    calcquant <- function(i)
    {
        parai <- c(as.list(b$estim[i, ]), as.list(b$fitpart$fix.arg))
        do.call(qdistname, c(list(p=probs), as.list(parai)))
    }
    
    bootquant <- sapply(1:b$nbboot, calcquant)
    if (length(probs)>1)
        bootquant <- as.data.frame(t(bootquant))
    else
        bootquant <- as.data.frame(bootquant)
	colnames(bootquant) <- paste("p=", probs, sep="")
	
    quantmedian <- rbind(apply(bootquant, MARGIN=2, median, na.rm=TRUE))
	colnames(quantmedian) <- paste("p=", probs, sep="")
	rownames(quantmedian) <- "estimate"
    
    if (CI.type == "two.sided")
    {
        alpha <- (1-CI.level)/2
        quantCI <- rbind(
                         apply(bootquant, MARGIN=2, quantile, alpha, na.rm=TRUE), 
                         apply(bootquant, MARGIN=2, quantile, 1-alpha, na.rm=TRUE))
        rownames(quantCI) <- format.perc(c(alpha, 1-alpha), 3)
    }else if (CI.type == "less")
    {
        quantCI <- t(apply(bootquant, MARGIN=2, quantile, CI.level, na.rm=TRUE))
        rownames(quantCI) <- format.perc(CI.level, 3)
    }else
    {
        quantCI <- t(apply(bootquant, MARGIN=2, quantile, 1-CI.level, na.rm=TRUE))
        rownames(quantCI) <- format.perc(1-CI.level, 3)
    }
    
    
    # message when lack of convergence
    nbconverg <- length(b$converg[b$converg == 0])
    
    reslist <- list(quantiles = basequant$quantiles, probs=probs, bootquant = bootquant, 
                    quantCI = as.data.frame(quantCI), quantmedian = quantmedian, 
                    CI.type =  CI.type, CI.level = CI.level, 
                    nbboot = b$nbboot, nbconverg = nbconverg)
    if(!cens)
        class(reslist) <- "quantile.bootdist"
    else
        class(reslist) <- "quantile.bootdistcens"
    reslist
}


print.quantile.bootdist <- function(x, ...)
{
    if (!inherits(x, "quantile.bootdist"))
        stop("Use only with 'quantile.bootdist' objects")
    typedata <- "(non-censored data)"
    
    #base quantiles
    cat("(original) estimated quantiles for each specified probability ", typedata,"\n", sep="")
    print(x$quantiles)      
    cat("Median of bootstrap estimates\n")
	print(x$quantmedian)
    
    #confidence intervals
    cat("\n")

    if (x$CI.type == "two.sided")
    {
        cat("two-sided ", format.perc(x$CI.level, 3)," CI of each quantile\n", sep="")
        print(x$quantCI)
    }else if (x$CI.type == "less")
    {
        cat("right bound of one-sided ", format.perc(x$CI.level, 3)," CI of each quantile\n")
        print(x$quantCI)
    }else
    {
        cat("left bound of one-sided ", format.perc(x$CI.level, 3)," CI of each quantile\n")
        print(x$quantCI)
    }   
    
    if (x$nbconverg < x$nbboot)
    {
        cat("\n")
        cat("The estimation method converged only for ", x$nbconverg, " among ", 
            x$nbboot, " bootstrap iterations.\n")
    }
    invisible(x)
}

print.quantile.bootdistcens <- function(x, ...)
{
    if (!inherits(x, "quantile.bootdistcens"))
        stop("Use only with 'quantile.bootdistcens' objects")
    typedata <- "(censored data)"
    
    #base quantiles
    cat("(original) estimated quantiles for each specified probability ", typedata,"\n", sep="")
    print(x$quantiles)      
    cat("Median of bootstrap estimates\n")
	print(x$quantmedian)
    
    #confidence intervals
    cat("\n")
    
    if (x$CI.type == "two.sided")
    {
        cat("two-sided ", format.perc(x$CI.level, 3)," CI of each quantile\n", sep="")
        print(x$quantCI)
    }else if (x$CI.type == "less")
    {
        cat("right bound of one-sided ", format.perc(x$CI.level, 3)," CI of each quantile\n")
        print(x$quantCI)
    }else
    {
        cat("left bound of one-sided ", format.perc(x$CI.level, 3)," CI of each quantile\n")
        print(x$quantCI)
    }   
    
    if (x$nbconverg < x$nbboot)
    {
        cat("\n")
        cat("The estimation method converged only for ", x$nbconverg, " among ", 
            x$nbboot, " bootstrap iterations.\n")
    }
    invisible(x)
}

#from the stat package (not exported in fitdistrplus)
format.perc <- function(probs, digits)
    ## Not yet exported, maybe useful in other contexts:
    ## quantile.default() sometimes uses a version of it
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits), "%")
