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
### fit parametric distributions for non-censored data
###
###         R functions
### 

fitdist<-function (data, distr, method=c("mle", "mme"), start, chisqbreaks, meancount, ...) 
{
    if (!is.character(distr)) 
		distname<-substring(as.character(match.call()$distr),2)
    else 
		distname<-distr
    ddistname<-paste("d",distname,sep="")
	
    if (!exists(ddistname,mode="function"))
        stop(paste("The ",ddistname," function must be defined"))
    pdistname<-paste("p",distname,sep="")
    if (!exists(pdistname,mode="function"))
        stop(paste("The ",pdistname," function must be defined"))
    
	method <- match.arg(method)
#	if (!is.element(method,c("mle","mom")))
#        stop("method must be affected to 'mle' or 'mom'") 
    
	if (!missing(start) & method=="mme")
        warnings("Starting values for parameters will not be used with matching moments")  
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")
    n<-length(data)
    # MLE fit with mledist or matching moments fit with momdist
    if (method=="mme")
    {
        estimate <- momdist(data,distname)
        sd <- NULL
        loglik <- NULL
        aic <- NULL
        bic <- NULL
        correl <- NULL
    }
    else
    {
        if (missing(start))
            mle<-mledist(data,distname,...) 
        else 
            mle<-mledist(data,distname,start,...)
        if (mle$convergence>0) 
           stop("the function mle failed to estimate the parameters, 
            with the error code ",mle$convergence) 
        estimate<-mle$estimate
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
    # Goodness of fit statistics
    if (is.element(distname,c("binom","nbinom","geom","hyper","pois"))) 
        discrete<-TRUE
    else 
        discrete<-FALSE
    
    # chi-squared statistic
    if (missing(chisqbreaks)) { 
        if (missing(meancount))
            meancount<-round(n/((4*n)^(2/5)))
        sdata<-sort(data)
        if (length(sdata)>ceiling(1.5*meancount)) {
            limit<-sdata[meancount]
            sdata<-sdata[sdata>limit]
            chisqbreaks<-limit
        }
        else {
            warnings("The sample is too small to automatically define chisqbreaks")
            chisq<-NULL
            chisqbreaks<-NULL
            chisqpvalue<-NULL
            chisqtable<-NULL
            chisqdf<-NULL
            
        }
        while (length(sdata)>ceiling(1.5*meancount)) {
            limit<-sdata[meancount]
            sdata<-sdata[sdata>limit]
            chisqbreaks<-c(chisqbreaks,limit)
        } 
    }
                
    if (!is.null(chisqbreaks)) {
        if(!is.numeric(chisqbreaks))
            stop("chisqbreaks must be a numeric vector defining the cell boundaries")
        nbreaks<-length(chisqbreaks)  
        pbreaks<-do.call(pdistname,c(list(q=chisqbreaks),as.list(estimate)))
        Fobsbreaks<-ecdf(data)(chisqbreaks)
        Fobsunder<-c(0,Fobsbreaks[1:nbreaks-1]) 
        punder<-c(0,pbreaks[1:nbreaks-1])
        if (pbreaks[nbreaks]==1 & Fobsbreaks[nbreaks]==1) {
            p<-pbreaks-punder
            Fobs<-Fobsbreaks-Fobsunder
        }
        else {
            p<-c(pbreaks-punder,1-pbreaks[nbreaks])
            Fobs<-c(Fobsbreaks-Fobsunder,1-Fobsbreaks[nbreaks])            
        }
        obscounts<-round(Fobs*n)
        theocounts<-p*n
        chisq<-sum(((obscounts-theocounts)^2)/theocounts)
        chisqdf<-length(obscounts)-1-length(estimate)
        if (chisqdf>0) {
            chisqpvalue<-pchisq(chisq,df=chisqdf,lower.tail=FALSE)
        }
        else
            chisqpvalue<-NULL
        chisqtable<-as.table(cbind(obscounts,theocounts))
        for (i in 1:length(obscounts)-1)
            rownames(chisqtable)[i]<-paste("<=",signif(chisqbreaks[i],digits=4))
        rownames(chisqtable)[length(obscounts)]<-paste(">",signif(chisqbreaks[i],digits=4))
   }
        
    if (!discrete) {
        # Kolmogorov-Smirnov statistic
        s<-sort(data)
        obspu<-seq(1,n)/n
        obspl<-seq(0,n-1)/n
        theop<-do.call(pdistname,c(list(q=s),as.list(estimate)))
        ks<-max(pmax(abs(theop-obspu),abs(theop-obspl)))
        Dmod<-ks*(sqrt(n)+0.12+0.11/sqrt(n))
        if (n>=30)
            kstest<-ifelse(Dmod>1.358,"rejected","not rejected")
        else
            kstest<-NULL
        
        # Anderson-Darling statistic
        ad<- -n-sum((2*(1:n)-1)*log(theop) + (2*n+1-2*(1:n))*log(1-theop))/n 
        if ((distname == "norm" | distname == "lnorm") & n>=5) {
            a2mod<-ad*(1+0.75/n+2.25/n^2)
            adtest<-ifelse(a2mod>0.752,"rejected","not rejected")
        } 
        else
            if (distname == "exp" & n>=5) {
                a2mod<-ad*(1+0.6/n)
                adtest<-ifelse(a2mod>1.321,"rejected","not rejected")
            }
            else
                if (distname == "gamma" & n>=5) {
                    m<-as.list(estimate)$shape
                    interp<-approxfun(c(1,2,3,4,5,6,8,10,12,15,20),
                    c(0.786,0.768,0.762,0.759,0.758,0.757,0.755,0.754,0.754,0.754,0.753),
                    yright=0.752)
                    adtest<-ifelse(ad>interp(m),"rejected","not rejected")
                }
                else
                    if (distname == "weibull" & n>=5) {
                        a2mod<-ad*(1+0.2/sqrt(n))
                        adtest<-ifelse(a2mod>0.757,"rejected","not rejected")
                    }
                    else
                        if (distname == "logis" & n>=5) {
                            a2mod<-ad*(1+0.25/n)
                            adtest<-ifelse(a2mod>0.66,"rejected","not rejected")
                        }
                        else
                            if (distname == "cauchy" & n>=5) {
                                interp<-approxfun(c(5,8,10,12,15,20,25,30,40,50,60,100),
                                c(1.77,3.2,3.77,4.14,4.25,4.05,3.57,3.09,2.48,2.14,1.92,1.52),
                                yright=1.225)
                                adtest<-ifelse(ad>interp(n),"rejected","not rejected")
                            }
                            else adtest<-NULL
        
        if (length(table(data))!=length(data))
        warnings("Kolmogorov-Smirnov and Anderson-Darling statistics may not be correct with ties")
    }
    else { # so if discrete
        ks<-NULL
        kstest<-NULL
        ad<-NULL
        adtest<-NULL
    }
    
    return(structure(list(estimate = estimate, method = method, sd = sd, 
    cor = correl, loglik = loglik, aic=aic, bic=bic,
    n = n, data=data, distname=distname,chisq = chisq, chisqbreaks=chisqbreaks,
    chisqpvalue=chisqpvalue,
    chisqdf=chisqdf,chisqtable=chisqtable, 
    ad = ad,adtest=adtest,ks = ks,kstest=kstest), class = "fitdist"))
        
}

print.fitdist <- function(x,...){
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    if (x$method=="mme") 
        cat("Fitting of the distribution '",x$distname,"' by matching moments \n")
    else
       cat("Fitting of the distribution '",x$distname,"' by maximum likelihood \n")
    cat("Parameters:\n")
    op<-options()
    options(digits=3)
    print(data.frame("estimate" = x$estimate),...)
    options(op)

}

plot.fitdist <- function(x,breaks="default",...){
    if (!inherits(x, "fitdist"))
        stop("Use only with 'fitdist' objects")
    plotdist(data=x$data,distr=x$distname,
    para=as.list(x$estimate),breaks=breaks,...)
}

summary.fitdist <- function(object,...){
    if (!inherits(object, "fitdist"))
        stop("Use only with 'fitdist' objects")
    ddistname<-paste("d",object$distname,sep="")
    pdistname<-paste("p",object$distname,sep="")
    
    op<-options()
    options(digits=3)
    if (object$method=="mme") 
        cat("FITTING OF THE DISTRIBUTION '",object$distname,"' BY MATCHING MOMENTS \n")
    else
       cat("FITTING OF THE DISTRIBUTION '",object$distname,"' BY MAXIMUM LIKELIHOOD \n")
    cat("PARAMETERS\n")
    if (object$method=="mle") {
        print(cbind.data.frame("estimate" = object$estimate, "Std. Error" = object$sd))
        cat("Loglikelihood: ",object$loglik,"  ")
        cat("AIC: ",object$aic,"  ")
        cat("BIC: ",object$bic,"\n")
        if (length(object$estimate) > 1) {
            cat("Correlation matrix:\n")
            print(object$cor)
            cat("\n")
        }
    }
    else {
        print(cbind.data.frame("estimate" = object$estimate))
    }
    cat("------\n")
    cat("GOODNESS-OF-FIT STATISTICS \n")
    cat("\n")
    cat("_____________ Chi-squared_____________\n")
    if(!is.null(object$chisq)) {
        cat("Chi-squared statistic: ",object$chisq,"\n")
        cat("Degree of freedom of the Chi-squared distribution: ",object$chisqdf,"\n")
        if (object$chisqdf<=0) {
            cat("!!! The degree of freedom of the chi-squared distribution is less than 1 !!! \n") 
            cat("The number of cells is insufficient to calculate the p-value.  \n") 
        }
        else
        { 
            cat("Chi-squared p-value: ",object$chisqpvalue,"\n")
            if (any(object$chisqtable[,2]<5)) cat("!!! the p-value may be wrong 
            with some theoretical counts < 5 !!! \n")
        }
    }
    else cat("The sample is too small to automatically define cells for Chi-squared test \n")
    if (!is.null(object$ks) & !is.null(object$chisq)) {
        cat("\n")
        cat("!!! For continuous distributions, Kolmogorov-Smirnov and  \n")
        cat("      Anderson-Darling statistics should be prefered !!! \n")
    }
    if(!is.null(object$ks)) {
        cat("\n")
        cat("_____________ Kolmogorov-Smirnov_____________\n")
        cat("Kolmogorov-Smirnov statistic: ",object$ks,"\n")
        if (!is.null(object$kstest)) {
            cat("Kolmogorov-Smirnov test: ",object$kstest,"\n")
            cat("!!! The result of this test may be too conservative as it  \n")
            cat("     assumes that the distribution parameters are known !!! \n")
        }
        else
            cat("Kolmogorov-Smirnov test: not calculated \n")
    }
    if(!is.null(object$ad)) {
        cat("\n")
        cat("_____________ Anderson-Darling_____________\n")
        cat("Anderson-Darling statistic: ",object$ad,"\n")
        if (!is.null(object$adtest)) 
            cat("Anderson-Darling test: ",object$adtest,"\n")
        else
            cat("Anderson-Darling test: not calculated \n")
    }
    options(op)
}
