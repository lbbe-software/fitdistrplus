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
### fit parametric distributions for censored data
###
###			R functions
### 

fitdistcens<-function (censdata, distr, start) 
{
    if (!is.character(distr)) distname<-substring(as.character(match.call()$distr),2)
    else distname<-distr
    ddistname<-paste("d",distname,sep="")
    if (!exists(ddistname,mode="function"))
        stop(paste("The ",ddistname," function must be defined"))
    pdistname<-paste("p",distname,sep="")
    if (!exists(pdistname,mode="function"))
        stop(paste("The ",pdistname," function must be defined"))
    # MLE fit with mledistcens 
    if (missing(start))
        mle<-mledistcens(censdata,distname) 
    else 
    mle<-mledistcens(censdata,distname,start)
    if (mle$convergence>0) 
        stop("the function mle failed to estimate the parameters, 
        with the error code ",mle$convergence) 
    estimate<-mle$estimate
    varcovar<-solve(mle$hessian)
    sd<-sqrt(diag(varcovar))
    correl<-cov2cor(varcovar)
    loglik<-mle$loglik
         
    return(structure(list(estimate = estimate, sd = sd, cor = correl, loglik = loglik, 
        censdata=censdata, distname=distname), class = "fitdistcens"))
        
}

print.fitdistcens <- function(x,...){
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    cat("Fitting of the distribution '",x$distname,"' on censored data by maximum likelihood \n")
    cat("Parameters:\n")
    op<-options()
    options(digits=3)
    print(data.frame("estimate" = x$estimate),...)
    options(op)

}

plot.fitdistcens <- function(x,...){
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    plotdistcens(censdata=x$censdata,distr=x$distname,
    para=as.list(x$estimate),...)
}

summary.fitdistcens <- function(object,...){
    if (!inherits(object, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    ddistname<-paste("d",object$distname,sep="")
    pdistname<-paste("p",object$distname,sep="")
    
    op<-options()
    options(digits=3)
    cat("FITTING OF THE DISTRIBUTION '",object$distname,
    "' BY MAXIMUM LIKELIHOOD ON CENSORED DATA \n")
    cat("PARAMETERS\n")
    print(cbind.data.frame("estimate" = object$estimate, "Std. Error" = object$sd))
    cat("Loglikelihood: ",object$loglik,"\n")
    if (length(object$estimate) > 1) {
        cat("Correlation matrix:\n")
        print(object$cor)
        cat("\n")
    }
    options(op)
}
