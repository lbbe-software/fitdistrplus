#############################################################################
#   Copyright (c) 2012 Christophe Dutang
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
### plot density functions for various fits
### of continuous distribution(s) (fitdist results)
### on a same dataset
###
###         R functions
###


denscomp <- function(ft, addlegend=TRUE, legendtext, datapch, datacol, probability,
	fitcol, fitlty, xlab, ylab, xlim, ylim, main, xlegend = "topright", ylegend = NULL,
	..., demp=FALSE, dempcol="grey")
{
	if(inherits(ft, "fitdist"))
	{
		ft <- list(ft)
	}else if(length(ft) == 1)
	{
		if(!inherits(ft, "fitdist"))
		stop("argument ft must a 'fitdist' object or a list of 'fitdist' objects.")
	}else if(!is.list(ft))
	{
		stop("argument ft must be a list of 'fitdist' objects")
	}else
	{
		if(any(sapply(ft, function(x) !inherits(x, "fitdist"))))		
		stop("argument ft must be a list of 'fitdist' objects")
	}
	
	
    nft <- length(ft)
    if (missing(datapch)) datapch <- 16
    if (missing(datacol)) datacol <- NULL    
	if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (length(fitcol)!=nft)
		stop("if specified, fitcol must be a vector of length
		 the number of fitted distributions to represent")
    if (missing(fitlty)) fitlty <- 1:nft
    if (length(fitlty) != nft)
		stop("if specified, fitlty must be a vector of length
		 the number of fitted distributions to represent")
	if(missing(probability))
		probability <- TRUE
    if (missing(xlab))
		xlab <- "data"
    if (missing(ylab)) 
		ylab <- ifelse(probability, "Density", "Frequency")
    if (missing(main)) 
		main <- ifelse(probability, "Histogram and theoretical densities", 
					   "Histogram and theoretical frequencies")
	
    mydata <- ft[[1]]$data
	
    if (missing(xlim))
    {
        xmin <- min(mydata)
        xmax <- max(mydata)
        xlim <- range(mydata)
    }
    else
    {
        xmin <- xlim[1]
        xmax <- xlim[2]
    }
	if (missing(ylim))
	{
		options(warn=-1)
		reshist <- hist(mydata, plot=FALSE, probability=probability)
		options(warn=getOption("warn"))
		if(!probability)
			ylim <- c(0, max(reshist$counts))
		else
			ylim <- c(0, max(reshist$density))
	}
	
    verif.ftidata <- function(fti)
    {
        if (any(fti$data != mydata))
			stop("All compared fits must have been obtained with the same dataset")
		invisible()
    }
	lapply(ft, verif.ftidata)

	n <- length(mydata)
	sfin <- seq(xmin, xmax, by=(xmax-xmin)/100)
	scalefactor <- ifelse(probability, 1, n * diff(reshist$breaks))

	
# computation of each fitted distribution
    comput.fti <- function(i, ...)
    {
        fti <- ft[[i]]
        para <- c(as.list(fti$estimate),as.list(fti$fix.arg))
        distname <- fti$distname
        ddistname <- paste("d",distname,sep="")
		
		do.call(ddistname, c(list(x=sfin), as.list(para))) * scalefactor
    }
    fitteddens <- sapply(1:nft, comput.fti, ...)
	
	
	reshist <- hist(mydata, main=main, xlab=xlab, ylab=ylab, xlim=xlim, 
					ylim=range(ylim, fitteddens), col=datacol, probability=probability, ...)
	for(i in 1:nft)
		lines(sfin, fitteddens[,i], lty=fitlty[i], col=fitcol[i], ...)
    
	if(demp)
		lines(density(mydata)$x, density(mydata)$y * scalefactor, col=dempcol)
	
    if (addlegend)
    {
        if (missing(legendtext) && !demp) 
		{
			legendtext <- paste("fit", 1:nft)
		}else if (missing(legendtext) && demp) 
		{
			legendtext <- c(paste("fit", 1:nft), "emp.")
			fitlty <- c(fitlty, 1)
			fitcol <- c(fitcol, dempcol)
		}else if(demp)
		{
			legendtext <- c(legendtext, "emp.")
			fitlty <- c(fitlty, 1)
			fitcol <- c(fitcol, dempcol)
		}	

        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext,
			   lty=fitlty, col=fitcol, ...)
    }
}
