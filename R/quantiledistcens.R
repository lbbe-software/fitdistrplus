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
### fitted on non-censored data
### and associated bootstrap confidence intervals
###
###         R functions
### 
quantiledistcens <- function(f, probs = seq(0.1,0.9,0.1), bootstrap = TRUE, 
                bootstrap.arg = list(niter=1001))
{
    if (!inherits(f, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    
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
        resbootdist <- do.call(bootdistcens,c(list(f),as.list(bootstrap.arg)))
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
        
        quantCI <- rbind(
        apply(bootquant,MARGIN=2,quantile,0.025,na.rm=TRUE),
        apply(bootquant,MARGIN=2,quantile,0.975,na.rm=TRUE))
        rownames(quantCI) <- c("2.5%","97.5%")
        
        cat("\n")
        cat("95% percentile CI (bootstrap) of each quantile\n")
        print(quantCI)

        reslist <- list(quantiles = quantiles, resbootdistcens = resbootdist, bootquant = bootquant, quantCI = as.data.frame(quantCI))

    }
    else # if (bootstrap)
    {
        reslist <- list(quantiles = quantiles, resbootdistcens = NULL, bootquant = NULL, quantCI = NULL)
    }
        
    invisible(reslist)
}
