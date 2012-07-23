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


qqcomp <- function(ft, addlegend=TRUE, legendtext, datapch, datacol, xlogscale=FALSE, ylogscale=FALSE,
	fitcol, fitlty, xlab, ylab, xlim, ylim, main, xlegend = "bottomright", ylegend = NULL,
	..., use.ppoints = TRUE, a.ppoints = 0.5)
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
	if (missing(datapch)) datapch <- 21
    if (missing(datacol)) datacol <- "black"    
	if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (length(fitcol) != nft)
		stop("if specified, fitcol must be a vector of length
		 the number of fitted distributions to represent")
    if (missing(fitlty)) fitlty <- 1:nft
    if (length(fitlty) != nft)
		stop("if specified, fitlty must be a vector of length
		 the number of fitted distributions to represent")
    if (missing(xlab))
		xlab <- "data"
    if (missing(ylab)) 
		ylab <- "Probability"
    if (missing(main)) 
		main <- "Empirical and theoretical quantiles"
	
    mydata <- ft[[1]]$data
	
    if (missing(xlim))
    {
        xmin <- min(mydata)
        xmax <- max(mydata)
        xlim <- range(mydata)
    }else
    {
        xmin <- xlim[1]
        xmax <- xlim[2]
    }

    verif.ftidata <- function(fti)
    {
        if (any(fti$data != mydata))
			stop("All compared fits must have been obtained with the same dataset")
		invisible()
    }
	lapply(ft, verif.ftidata)

	n <- length(mydata)
	if (use.ppoints)
		obsp <- ppoints(n, a = a.ppoints)
    else
		obsp <- (1:n) / n
	if (missing(ylim))
		ylim <- range(obsp)

	if(ylim[1] == 0)
		ylim[1] <- tail(sort(obsp), -1)[1] * 0.1
	probseq <- seq(ylim[1], ylim[2],by=(ylim[2]-ylim[1])/101)
	
	
	logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")

# computation of each fitted distribution
    comput.fti <- function(i, ...)
    {
        fti <- ft[[i]]
        para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
        distname <- fti$distname

        qdistname <- paste("q", distname, sep="")
		do.call(qdistname, c(list(p=probseq), as.list(para)))
	}
    fittedquant <- sapply(1:nft, comput.fti, ...)
			
	
	resquant <- plot(sort(mydata), obsp, main=main, xlab=xlab, ylab=ylab, log=logxy,
			pch=datapch, xlim=range(xlim, fittedquant), ylim=ylim, col=datacol, ...)
	for(i in 1:nft)
		lines(fittedquant[,i], probseq, lty=fitlty[i], col=fitcol[i], ...)
    
    if (addlegend)
    {
        if (missing(legendtext)) 
			legendtext <- paste("fit",1:nft)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
    }
}
