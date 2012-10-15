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
### calculate goodness-of-fit statistics for 
### fit of a parametric distribution on non-censored data
###
###         R functions
### 

gofstat <- function (f, chisqbreaks, meancount, print.test = FALSE) 
{

    if (!inherits(f, "fitdist"))
        stop("Use only with 'fitdist' objects")
        
    data <- f$data
    n <- f$n
    distname <- f$distname
    pdistname <- paste("p", distname,sep="")
    estimate <- f$estimate
    fix.arg <- f$fix.arg

    # Goodness of fit statistics
    if (is.element(distname,c("binom","nbinom","geom","hyper","pois"))) 
        discrete <- TRUE
    else 
        discrete <- FALSE
    
    # chi-squared statistic
    if (missing(chisqbreaks)) { 
        if (missing(meancount))
            meancount <- round(n/((4*n)^(2/5)))
        sdata <- sort(data)
        if (length(sdata)>ceiling(1.5*meancount)) {
            limit <- sdata[meancount]
            sdata <- sdata[sdata>limit]
            chisqbreaks <- limit
        }
        else {
            warnings("The sample is too small to automatically define chisqbreaks")
            chisq <- NULL
            chisqbreaks <- NULL
            chisqpvalue <- NULL
            chisqtable <- NULL
            chisqdf <- NULL
            
        }
        while (length(sdata)>ceiling(1.5*meancount)) {
            limit <- sdata[meancount]
            sdata <- sdata[sdata>limit]
            chisqbreaks <- c(chisqbreaks,limit)
        } 
    }
                
    if (!is.null(chisqbreaks)) {
        if(!is.numeric(chisqbreaks))
            stop("chisqbreaks must be a numeric vector defining the cell boundaries")
        nbreaks <- length(chisqbreaks)  
        pbreaks <- do.call(pdistname,c(list(q=chisqbreaks),as.list(estimate),fix.arg))
        Fobsbreaks <- ecdf(data)(chisqbreaks)
        Fobsunder <- c(0,Fobsbreaks[1:nbreaks-1]) 
        punder <- c(0,pbreaks[1:nbreaks-1])
        if (pbreaks[nbreaks]==1 & Fobsbreaks[nbreaks]==1) {
            p <- pbreaks-punder
            Fobs <- Fobsbreaks-Fobsunder
        }
        else {
            p <- c(pbreaks-punder,1-pbreaks[nbreaks])
            Fobs <- c(Fobsbreaks-Fobsunder,1-Fobsbreaks[nbreaks])            
        }
        obscounts <- round(Fobs*n)
        theocounts <- p*n
        chisq <- sum(((obscounts-theocounts)^2)/theocounts)
        chisqdf <- length(obscounts)-1-length(estimate)
        if (chisqdf>0) {
            chisqpvalue <- pchisq(chisq,df=chisqdf,lower.tail=FALSE)
        }
        else
            chisqpvalue <- NULL
        chisqtable <- as.table(cbind(obscounts,theocounts))
        for (i in 1:length(obscounts)-1)
            rownames(chisqtable)[i] <- paste("<=",signif(chisqbreaks[i],digits=4))
        rownames(chisqtable)[length(obscounts)] <- paste(">",signif(chisqbreaks[i],digits=4))
   }
        
    if (!discrete) {
        s <- sort(data)
        obspu <- seq(1,n)/n
        obspl <- seq(0,n-1)/n
        theop <- do.call(pdistname,c(list(q=s),as.list(estimate),fix.arg))
        # Kolmogorov-Smirnov statistic
        ks <- max(pmax(abs(theop-obspu),abs(theop-obspl)))
        Dmod <- ks*(sqrt(n)+0.12+0.11/sqrt(n))
        if (n>=30)
            kstest <- ifelse(Dmod>1.358,"rejected","not rejected")
        else
            kstest <- NULL
        
        # Anderson-Darling statistic
        ad <- - n - mean( (2 * seq(1:n) - 1) * (log(theop) + log(1 - rev(theop))) ) 
        # ad <-  -n-sum((2*(1:n)-1)*log(theop) + (2*n+1-2*(1:n))*log(1-theop))/n 
        if (is.null(fix.arg) & f$method == "mle")
        {
            # the following test does not correspond to MLE estimate but to unbiased 
            # estimate of the variance
          #if ((distname == "norm" | distname == "lnorm") & n>=5) {
          #  a2mod <- ad*(1+0.75/n+2.25/n^2)
          #  adtest <- ifelse(a2mod>0.752,"rejected","not rejected")
          #} 
          #else
            if (distname == "exp" & n>=5) {
                a2mod <- ad*(1+0.6/n)
                adtest <- ifelse(a2mod>1.321,"rejected","not rejected")
            }
            else
                if (distname == "gamma" & n>=5) {
                    m <- as.list(estimate)$shape
                    interp <- approxfun(c(1,2,3,4,5,6,8,10,12,15,20),
                    c(0.786,0.768,0.762,0.759,0.758,0.757,0.755,0.754,0.754,0.754,0.753),
                    yright=0.752)
                    adtest <- ifelse(ad>interp(m),"rejected","not rejected")
                }
                else
                    if (distname == "weibull" & n>=5) {
                        a2mod <- ad*(1+0.2/sqrt(n))
                        adtest <- ifelse(a2mod>0.757,"rejected","not rejected")
                    }
                    else
                        # the following test does not correspond to MLE estimate  
                        # if (distname == "logis" & n>=5) {
                        #    a2mod <- ad*(1+0.25/n)
                        #    adtest <- ifelse(a2mod>0.66,"rejected","not rejected")
                        #}
                        # else
                            if (distname == "cauchy" & n>=5) {
                                interp <- approxfun(c(5,8,10,12,15,20,25,30,40,50,60,100),
                                c(1.77,3.2,3.77,4.14,4.25,4.05,3.57,3.09,2.48,2.14,1.92,1.52),
                                yright=1.225)
                                adtest <- ifelse(ad>interp(n),"rejected","not rejected")
                            }
                            else adtest <- NULL
        }  # if (is.null(fix.arg)...)
        else 
            adtest <- NULL
            
                    # Cramer-von Mises statistic
        cvm <- 1/(12*n) + sum( ( theop - (2 * seq(1:n) - 1)/(2 * n) )^2 )
        
        if (is.null(fix.arg) & f$method == "mle")
        {
            # the following test does not correspond to MLE estimate but to unbiased 
            # estimate of the variance
          # if ((distname == "norm" | distname == "lnorm") & n>=5) {
          #  w2mod <- cvm*(1+0.5/n)
          #  cvmtest <- ifelse(w2mod>0.126,"rejected","not rejected")
          # } 
          # else
            if (distname == "exp" & n>=5) {
                w2mod <- cvm*(1+0.16/n)
                cvmtest <- ifelse(w2mod>0.222,"rejected","not rejected")
            }
            else
                if (distname == "gamma" & n>=5) {
                    m <- as.list(estimate)$shape
                    interp <- approxfun(c(1,2,3,4,5,6,8,10,12,15,20),
                    c(0.136,0.131,0.129,0.128,0.128,0.128,0.127,0.127,0.127,0.127,0.126),
                    yright=0.126)
                    cvmtest <- ifelse(cvm>interp(m),"rejected","not rejected")
                }
                else
                    if (distname == "weibull" & n>=5) {
                        w2mod <- cvm*(1+0.2/sqrt(n))
                        cvmtest <- ifelse(w2mod>0.124,"rejected","not rejected")
                    }
                    else
                        # the following test does not correspond to MLE estimate 
                        # if (distname == "logis" & n>=5) {
                        #     w2mod <- (n*cvm - 0.08)/(n - 1)
                        #     cvmtest <- ifelse(w2mod>0.098,"rejected","not rejected")
                        # }
                        # else
                            if (distname == "cauchy" & n>=5) {
                                interp <- approxfun(c(5,8,10,12,15,20,25,30,40,50,60,100),
                                c(0.393,0.703,0.833,0.896,0.904,0.835,0.726,0.615,0.460,0.381,0.330,0.2378),
                                yright=0.170)
                                cvmtest <- ifelse(cvm>interp(n),"rejected","not rejected")
                            }
                            else cvmtest <- NULL
        }  # if (is.null(fix.arg))
        else 
            cvmtest <- NULL
            
        if (length(table(data))!=length(data))
        warnings("Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling statistics may not be correct with ties")
    }
    else { # so if discrete
        ks <- NULL
        kstest <- NULL
        cvm <- NULL
        cvmtest <- NULL
        ad <- NULL
        adtest <- NULL
    }

    res<-list(chisq = chisq, chisqbreaks=chisqbreaks,
    chisqpvalue = chisqpvalue,
    chisqdf = chisqdf,chisqtable = chisqtable, 
    cvm = cvm,cvmtest = cvmtest,ad = ad,adtest = adtest,ks = ks,kstest=kstest)
    
    if (discrete) 
    {
        if(!is.null(chisq)) 
        {
            cat("Chi-squared statistic: ",chisq,"\n")
            if (print.test) 
            {
                cat("Degree of freedom of the Chi-squared distribution: ",chisqdf,"\n")
                if (chisqdf<=0) 
                {
                    cat("  The degree of freedom of the chi-squared distribution is less than 1  \n") 
                    cat("  The number of cells is insufficient to calculate the p-value.  \n") 
                }
                else
                { 
                    cat("Chi-squared p-value: ",chisqpvalue,"\n")
                    if (any(chisqtable[,2]<5)) 
                    cat("   the p-value may be wrong with some theoretical counts < 5  \n")
                }
            }
        }
        else cat("The sample is too small to automatically define cells for Chi-squared test \n")
    }
    else 
    { # if (discrete)
        cat("Kolmogorov-Smirnov statistic: ",ks,"\n")
        if (print.test) 
        {
            if (!is.null(kstest)) 
            {
                cat("Kolmogorov-Smirnov test: ",kstest,"\n")
                cat("   The result of this test may be too conservative as it  \n")
                cat("   assumes that the distribution parameters are known\n")
            }
            else
                cat("Kolmogorov-Smirnov test: not calculated \n")
        }
        cat("Cramer-von Mises statistic: ",cvm,"\n")
        if (print.test) 
        {
            if (!is.null(cvmtest)) 
                cat("Cramer-von Mises test: ",cvmtest,"\n")
            else
                cat("Crame-von Mises test: not calculated \n")
        }
        cat("Anderson-Darling statistic: ",ad,"\n")
        if (print.test) 
        {
            if (!is.null(adtest)) 
                cat("Anderson-Darling test: ",adtest,"\n")
            else
                cat("Anderson-Darling test: not calculated \n")
        }
    }

    invisible(res)        
}
