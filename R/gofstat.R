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
### calculate goodness-of-fit statistics for 
### fit of a parametric distribution on non-censored data
###
###         R functions
### 

gofstat <- function (f, chisqbreaks, meancount, discrete, 
	fitnames=NULL) 
{
	if(inherits(f, "fitdist"))
    {
        f <- list(f)
    }else if(length(f) == 1)
    {
        if(!inherits(f, "fitdist"))
		stop("argument f must a 'fitdist' object or a list of 'fitdist' objects.")
    }else if(!is.list(f))
    {
        stop("argument f must be a list of 'fitdist' objects")
    }else
    {
        if(any(sapply(f, function(x) !inherits(x, "fitdist"))))        
		stop("argument f must be a list of 'fitdist' objects")
    }
	        
	odata <- f[[1]]$data
	sdata <- sort(odata)
	n <- f[[1]]$n
	distname <- f[[1]]$distname
	pdistname <- paste("p", distname, sep="")
	estimate <- f[[1]]$estimate
	fix.arg <- f[[1]]$fix.arg
	
    verif.ftidata <- function(fti)
    {
        if (any(fti$data != odata))
			stop("All compared fits must have been obtained with the same dataset")
        invisible()
    }
    lapply(f, verif.ftidata)	
	

    # initiate discrete if not given 
	if(missing(discrete))
	{
      discrete <- f[[1]]$discrete
	}
	if(!is.logical(discrete))
		stop("wrong argument 'discrete'.")
        
	#define chisqbreaks if not defined
	if (missing(chisqbreaks)) { 
		if (missing(meancount))
		meancount <- round( n / ((4*n)^(2/5)) )
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
		sdata <- sort(odata)
	}
	
	nbfit <- length(f)

	
	if(is.null(fitnames))
		fitnames <- paste(1:nbfit, sapply(f, function(x) x$method), 
						  sapply(f, function(x) x$distname), sep="-")
	

	
	Chi2 <- compute.gofstat.Chi2(sdata, n, distname, pdistname, estimate, fix.arg, 
								 chisqbreaks)
	if(length(f) > 1)
	{
		#renaming
		names(Chi2$chisq) <- names(Chi2$chisqpvalue) <- names(Chi2$chisqdf) <- fitnames[1]
		colnames(Chi2$chisqtable)[2] <- paste("theo", fitnames[1], sep=" ")
		
		#computation and storing
		for(i in 2:nbfit)
		{
			Chi2temp <- compute.gofstat.Chi2(sdata, n, f[[i]]$distname, paste("p", f[[i]]$distname, sep=""), 
											 f[[i]]$estimate, f[[i]]$fix.arg, chisqbreaks)
			
			names(Chi2temp$chisq) <- names(Chi2temp$chisqpvalue) <- names(Chi2temp$chisqdf) <- fitnames[i]
			Chi2$chisq <- c(Chi2$chisq, Chi2temp$chisq)
			Chi2$chisqpvalue <- c(Chi2$chisqpvalue, Chi2temp$chisqpvalue)
			Chi2$chisqdf <- c(Chi2$chisqdf, Chi2temp$chisqdf)
			
			Chi2$chisqtable <- cbind(Chi2$chisqtable, Chi2temp$chisqtable[,2])
			colnames(Chi2$chisqtable)[NCOL(Chi2$chisqtable)] <- paste("theo", fitnames[i], sep=" ")
		}
	}
 
	if(discrete)
	{	
		addres <- Chi2
		
	}else
	{
		KSCvMAD <- compute.gofstat.KSCvMAD(sdata, n, distname, pdistname, estimate, 
										   fix.arg, f[[1]]$method)
		#renaming
		names(KSCvMAD$cvm) <- names(KSCvMAD$ad) <- names(KSCvMAD$ks) <- fitnames[1]
		
		if(!is.null(KSCvMAD$cvmtest))
		names(KSCvMAD$cvmtest) <- names(KSCvMAD$adtest) <- names(KSCvMAD$kstest) <- fitnames[1]
		
		if(length(f) > 1)
		{			
			#computation and storing
			for(i in 2:nbfit)
			{
				KSCvMADtemp <- compute.gofstat.KSCvMAD(sdata, n, f[[i]]$distname, paste("p", f[[i]]$distname, sep=""), 
													   f[[i]]$estimate, f[[i]]$fix.arg, f[[i]]$method)

				names(KSCvMADtemp$cvm) <- names(KSCvMADtemp$ad) <- names(KSCvMADtemp$ks) <- fitnames[i]

				if(!is.null(KSCvMADtemp$cvmtest))
					names(KSCvMADtemp$cvmtest) <- names(KSCvMADtemp$adtest) <- names(KSCvMADtemp$kstest) <- fitnames[i]

				KSCvMAD$cvm <- c(KSCvMAD$cvm, KSCvMADtemp$cvm)
				KSCvMAD$cvmtest <- c(KSCvMAD$cvmtest, KSCvMADtemp$cvmtest)
				KSCvMAD$ad <- c(KSCvMAD$ad, KSCvMADtemp$ad)
				KSCvMAD$adtest <- c(KSCvMAD$adtest, KSCvMADtemp$adtest)
				KSCvMAD$ks <- c(KSCvMAD$ks, KSCvMADtemp$ks)
				KSCvMAD$kstest <- c(KSCvMAD$kstest, KSCvMADtemp$kstest)

			}
		}
		
		addres <- c(Chi2, KSCvMAD)
	}
	
	aics <- sapply(f, function(x) x$aic)
	names(aics) <- fitnames
	bics <- sapply(f, function(x) x$bic)
	names(bics) <- fitnames

    res <- c(addres, aic=list(aics), bic=list(bics), discrete=discrete, nbfit=nbfit)
	class(res) <- c("gofstat.fitdist", "fitdist")
	res
}


#----------------------------------------------------------------------
#KS, CvM, AD statistics : only for continuous distributions
compute.gofstat.KSCvMAD <- function(sdata, n, distname, pdistname, estimate, 
	fix.arg, method)
{
	
	obspu <- seq(1,n)/n
    obspl <- seq(0,n-1)/n
    theop <- do.call(pdistname, c(list(q=sdata), as.list(estimate), fix.arg))
	
	# Kolmogorov-Smirnov statistic
	ks <- max(pmax(abs(theop-obspu), abs(theop-obspl)))
	Dmod <- ks*(sqrt(n)+0.12+0.11/sqrt(n))
	# Kolmogorov-Smirnov test
	if (n>=30)
		kstest <- ifelse(Dmod>1.358,"rejected","not rejected")
	else
		kstest <- "not computed"
        
	# Anderson-Darling statistic
	ad <- - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ) 
	
		# ad <-  -n-sum((2*(1:n)-1)*log(theop) + (2*n+1-2*(1:n))*log(1-theop))/n 

	# Anderson-Darling test        
	if (is.null(fix.arg) & method == "mle")
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
			adtest <- ifelse(a2mod>1.321, "rejected", "not rejected")
		}
		else
		if (distname == "gamma" & n>=5) {
			m <- as.list(estimate)$shape
			interp <- approxfun(c(1,2,3,4,5,6,8,10,12,15,20),
								c(0.786,0.768,0.762,0.759,0.758,0.757,0.755,0.754,0.754,0.754,0.753),
								yright=0.752)
			adtest <- ifelse(ad>interp(m), "rejected", "not rejected")
		}
		else
		if (distname == "weibull" & n>=5) {
			a2mod <- ad*(1+0.2/sqrt(n))
			adtest <- ifelse(a2mod>0.757, "rejected", "not rejected")
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
			adtest <- ifelse(ad>interp(n), "rejected", "not rejected")
		}
		else adtest <- "not computed"
	}  # if (is.null(fix.arg)...)
	else 
		adtest <- "not computed"
	
	# Cramer-von Mises statistic
	cvm <- 1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )

	# Cramer-von Mises test
	if (is.null(fix.arg) & method == "mle")
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
		else cvmtest <- "not computed"
	}  # if (is.null(fix.arg))
	else 
		cvmtest <- "not computed"
	
	if (length(table(sdata)) != length(sdata))
		warnings("Kolmogorov-Smirnov, Cramer-von Mises and Anderson-Darling statistics may not be correct with ties")

	list(cvm = cvm, cvmtest = cvmtest, ad = ad,adtest = adtest, ks = ks, kstest=kstest)
}


#----------------------------------------------------------------------
#chi-squared statistic : both for continuous and discrete distributions
compute.gofstat.Chi2 <- function(sdata, n, distname, pdistname, estimate, fix.arg, chisqbreaks)
{
	# chi-squared statistic and pvalues
	if (!is.null(chisqbreaks)) {
		if(!is.numeric(chisqbreaks))
			stop("chisqbreaks must be a numeric vector defining the cell boundaries")

		nbreaks <- length(chisqbreaks)  
		pbreaks <- do.call(pdistname, c(list(q=chisqbreaks), as.list(estimate), fix.arg))
		Fobsbreaks <- ecdf(sdata)(chisqbreaks)
		
		Fobsunder <- c(0, Fobsbreaks[1:nbreaks-1]) 
		punder <- c(0, pbreaks[1:nbreaks-1])
		if (pbreaks[nbreaks]==1 & Fobsbreaks[nbreaks]==1) 
		{
			p <- pbreaks-punder
			Fobs <- Fobsbreaks-Fobsunder
		}else 
		{
			p <- c(pbreaks-punder, 1-pbreaks[nbreaks])
			Fobs <- c(Fobsbreaks-Fobsunder, 1-Fobsbreaks[nbreaks])            
		}
		
		obscounts <- round(Fobs*n)
		theocounts <- p*n
		
		chisq <- sum(((obscounts-theocounts)^2)/theocounts)
		chisqdf <- length(obscounts)-1-length(estimate)
		
		chisqpvalue <- ifelse(chisqdf>0, pchisq(chisq, df=chisqdf, lower.tail=FALSE), NULL)
	
		chisqtable <- as.table(cbind(obscounts, theocounts))
		for (i in 1:length(obscounts)-1)
			rownames(chisqtable)[i] <- paste("<=", signif(chisqbreaks[i], digits=4))
		rownames(chisqtable)[length(obscounts)] <- paste(">", signif(chisqbreaks[i], digits=4))
		
		return( list(chisq = chisq, chisqbreaks=chisqbreaks, chisqpvalue = chisqpvalue,
					 chisqdf = chisqdf, chisqtable = chisqtable) )
	}else
		return(NULL)
}


print.gofstat.fitdist <- function(x, ...)
{
    if (!inherits(x, "gofstat.fitdist"))
		stop("Use only with 'gofstat.fitdist' objects")
	
	if (x$discrete) #discrete distribution
    {
        if(!is.null(x$chisq)) 
        {
            cat("Chi-squared statistic: ",x$chisq,"\n")
			cat("Degree of freedom of the Chi-squared distribution: ",x$chisqdf,"\n")
			if (any(x$chisqdf <= 0))
			{
				cat("  The degree of freedom of the chi-squared distribution is less than 1  \n") 
				cat("  The number of cells is insufficient to calculate the p-value.  \n") 
			}
			else
			{ 
				cat("Chi-squared p-value: ",x$chisqpvalue,"\n")
				if (any(x$chisqtable[,-1] < 5)) 
					cat("   the p-value may be wrong with some theoretical counts < 5  \n")
			}

			cat("Chi-squared table:\n")
			print(x$chisqtable)	

			cat("\nGoodness-of-fit criteria\n")
			mm <- rbind(AIC=x$aic, BIC=x$bic)
			rownames(mm) <- c("Aikake's Information Criterion",
							  "Bayesian Information Criterion")
			print(mm)
        }
        else 
			cat("The sample is too small to automatically define cells for Chi-squared test \n")
    }else # continuous distribution
    { 
#		if(x$nbfit == 1)
#		{
#			cat("Kolmogorov-Smirnov statistic: ", x$ks,"\n")
#			if (!is.null(x$kstest)) 
#			{
#				cat("Kolmogorov-Smirnov test: ",x$kstest,"\n")
#				cat("   The result of this test may be too conservative as it  \n")
#				cat("   assumes that the distribution parameters are known\n")
#			}
#			else
#				cat("Kolmogorov-Smirnov test: not calculated \n")
#			
#			cat("Cramer-von Mises statistic: ",x$cvm,"\n")
#			if (!is.null(x$cvmtest)) 
#				cat("Cramer-von Mises test: ",x$cvmtest,"\n")
#			else
#				cat("Crame-von Mises test: not calculated \n")
#			
#			cat("Anderson-Darling statistic: ",x$ad,"\n")
#			if (!is.null(x$adtest)) 
#				cat("Anderson-Darling test: ",x$adtest,"\n")
#			else
#				cat("Anderson-Darling test: not calculated \n")
#			
#			cat("Aikake's Information Criterion: ", x$aic, "\n")
#			cat("Bayesian Information Criterion: ", x$bic, "\n")
#		}else
#		{
			cat("Goodness-of-fit statistics\n")
			mm <- rbind(KS=x$ks, CvM=x$cvm, AD=x$ad)
			rownames(mm) <- c("Kolmogorov-Smirnov statistic", "Cramer-von Mises statistic",
							  "Anderson-Darling statistic")
			print(mm)

			cat("\nGoodness-of-fit criteria\n")
			mm <- rbind(AIC=x$aic, BIC=x$bic)
			rownames(mm) <- c("Aikake's Information Criterion",
							  "Bayesian Information Criterion")
			print(mm)
#		}

    }
	
	invisible(x)
}