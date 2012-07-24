#############################################################################
#   Copyright (c) 2012 Marie Laure Delignette-Muller                                                                                                  
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
### and associated bootstrap confidence intervals
###
###         R functions
### 


quantile.fitdist <- function(x, probs = seq(0.1,0.9,0.1), bootstrap = TRUE, CI.type = "two.sided",
                bootstrap.arg = list(bootmethod="param", niter=1001), ...)
{
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    cens <- FALSE
    quantiles(f = x,probs = probs,bootstrap = bootstrap,CI.type = CI.type,bootstrap.arg = bootstrap.arg,cens)
}

quantile.fitdistcens <- function(x, probs = seq(0.1,0.9,0.1), bootstrap = TRUE, CI.type = "two.sided",
                bootstrap.arg = list(niter=1001), ...)
{
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    cens <- TRUE
    quantiles(f = x,probs = probs,bootstrap = bootstrap,CI.type = CI.type,bootstrap.arg = bootstrap.arg,cens)
}

quantiles <- function(f, probs, bootstrap, CI.type, bootstrap.arg, cens )
{    
    CI.type <- match.arg(CI.type, c("two.sided","less","greater"))
    
    qdistname<-paste("q",f$distname,sep="")
    if (!exists(qdistname,mode="function"))
        stop(paste("The ",qdistname," function must be defined")) 
     
    # 1/ calculation and print of quantiles using estimations of parameters   
    para=c(as.list(f$estimate),as.list(f$fix.arg))
    quantiles <- do.call(qdistname,c(list(p=probs),as.list(para)))
    if (length(probs)>1)
        quantiles <- as.data.frame(t(quantiles))
    else
        quantiles <- as.data.frame(quantiles)
    colnames(quantiles) <- paste("prob=",probs,sep="")
        
    cat("Estimated quantiles for each specified probability \n")
    print(quantiles)
        
    # 2/ calculation of bootstraped quantiles and 95 percent CI      
    if (bootstrap)
    {
        if (cens)
        resbootdist <- do.call(bootdistcens,c(list(f),as.list(bootstrap.arg)))
        else
        resbootdist <- do.call(bootdist,c(list(f),as.list(bootstrap.arg)))
        
        calcquant <- function(i)
        {
            parai=c(as.list(resbootdist$estim[i,]),as.list(f$fix.arg))
            do.call(qdistname,c(list(p=probs),as.list(parai)))
        }
        
      
        bootquant <- sapply(1:nrow(resbootdist$estim),calcquant)
        if (length(probs)>1)
            bootquant <- as.data.frame(t(bootquant))
        else
            bootquant <- as.data.frame(bootquant)
        colnames(bootquant) <- paste("prob=",probs,sep="")
        
        if (CI.type == "two.sided")
        {
            quantCI <- rbind(
            apply(bootquant,MARGIN=2,quantile,0.025,na.rm=TRUE),
            apply(bootquant,MARGIN=2,quantile,0.975,na.rm=TRUE))
            rownames(quantCI) <- c("2.5%","97.5%")
            cat("\n")
            cat("two-sided 95% CI of each quantile\n")
            print(quantCI)
        }
        else
        {
            if (CI.type == "less")
            {
                quantCI <- 
                t(apply(bootquant,MARGIN=2,quantile,0.95,na.rm=TRUE))
                rownames(quantCI) <- c("95%")
                cat("\n")
                cat("right bound of one-sided 95% CI of each quantile\n")
                print(quantCI)
            }
            else
            {
                quantCI <- 
                t(apply(bootquant,MARGIN=2,quantile,0.05,na.rm=TRUE))
                rownames(quantCI) <- c("5%")
                cat("\n")
                cat("left bound of one-sided 95% CI (bootstrap) of each quantile\n")
                print(quantCI)
            }
        }
        
        # message when lack of convergence
        nconverg<-length(resbootdist$converg[resbootdist$converg==0])
        if (nconverg < length(resbootdist$converg))
        {
            cat("\n")
            cat("The estimation method converged only for ",nconverg," among ",
                    length(resbootdist$converg)," bootstrap iterations \n")
        }

        reslist <- list(quantiles = quantiles, resbootdist = resbootdist, bootquant = bootquant, quantCI = as.data.frame(quantCI))

    }
    else # if (bootstrap)
    {
        reslist <- list(quantiles = quantiles, resbootdist = NULL, bootquant = NULL, quantCI = NULL)
    }
        
    invisible(reslist)
}
