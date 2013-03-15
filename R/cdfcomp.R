#############################################################################
#   Copyright (c) 2011 Marie Laure Delignette-Muller, Christophe Dutang
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
### plot cumulative distribution functions for various fits
### of continuous distribution(s) (fitdist results)
### on a same dataset
###
###         R functions
###


cdfcomp <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    datapch, datacol, fitlty, fitcol, addlegend = TRUE, legendtext, 
    xlegend = "bottomright", ylegend = NULL, horizontals = TRUE, verticals = FALSE, 
    use.ppoints = TRUE, a.ppoints = 0.5, lines01 = FALSE, discrete, ...)
{
    if(inherits(ft, "fitdist"))
    {
        ft <- list(ft)
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
    if (missing(datacol)) datacol <- "black"
    if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (missing(fitlty)) fitlty <- 1:nft
    fitcol <- rep(fitcol, length.out=nft)
    fitlty <- rep(fitlty, length.out=nft)

    if (missing(xlab))
        xlab <- ifelse(xlogscale, "data in log scale", "data")
    if (missing(ylab)) ylab <- "CDF"
    if (missing(main)) main <- paste("Empirical and theoretical CDFs")

    mydata <- ft[[1]]$data
	distname <- ft[[1]]$distname
    n <- length(mydata)
    s <- sort(mydata)
	
    if ((xlogscale == TRUE) & min(mydata) <= 0)
        stop("log transformation of data requires only positive
    values")

    if(missing(xlim))
    {
        xmin <- min(mydata)
        xmax <- max(mydata)
        xlim <- c(xmin, xmax)
    }
    else
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
    lapply( ft,verif.ftidata)
	
	# initiate discrete if not given 
	if(missing(discrete))
	{
		if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
			discrete <- TRUE
		else 
			discrete <- FALSE
	}
	if(!is.logical(discrete))
		stop("wrong argument 'discrete'.")
	
    
    # plot of data (ecdf)
    if(xlogscale && !discrete)
        sfin <- seq(log10(xmin), log10(xmax), by=(log10(xmax)-log10(xmin))/100)
    else if(!xlogscale && !discrete)
        sfin <- seq(xmin, xmax, by=(xmax-xmin)/100)
	else
	{
		sfin <- seq(xmin, xmax, by=1)
	}

    if (use.ppoints && !discrete)
        obsp <- ppoints(n,a = a.ppoints)
    else
        # previous version with no vizualisation of ex-aequos
        # obsp <- ecdf(s)(s) 
        obsp <- (1:n) / n
    
    
    # computation of each fitted distribution
    comput.fti <- function(i)
    {
        fti <- ft[[i]]
        para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
        distname <- fti$distname
        pdistname <- paste("p",distname,sep="")
        if(xlogscale && !discrete)
        {
            do.call(pdistname, c(list(q=10^sfin), as.list(para)))
        }else
        {
            do.call(pdistname, c(list(q=sfin), as.list(para)))
        }
    }
    fittedprob <- sapply(1:nft, comput.fti)  	
	if(missing(ylim))
        ylim <- range(obsp, fittedprob) 
	else
		ylim <- range(ylim) #in case of users enter a bad ylim
    
    logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
    #main plotting
    plot(s, obsp, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
         log=logxy, pch=datapch, col=datacol, type="p", ...)

    # optional add of horizontal and vbertical lines for step function
    if (horizontals)
    {
        xhleft <- s[-length(s)]
        xhright <- s[-1L]
        yh <- obsp[-length(s)]
        segments(xhleft, yh, xhright, yh, col=datacol,...)
        segments(s[length(s)], 1, xmax, 1, col=datacol, lty = 2, ...)
        
        segments(xmin, 0, s[1], 0, col=datacol, lty = 2, ...)
        if (verticals)
        {
           xv <-xhright
           yvdown <- yh
           yvup <- obsp[-1L] 
           segments(xv, yvdown, xv, yvup, col=datacol,...)
           segments(s[1], 0, s[1], obsp[1], col=datacol, ...)
        }
    }
    #plot fitted cdfs
    if(!xlogscale)
        for(i in 1:nft)
            lines(sfin, fittedprob[,i], lty=fitlty[i], col=fitcol[i], 
				  type=ifelse(discrete, "s", "l"), ...)
    if(xlogscale)
        for(i in 1:nft)
            lines(10^sfin, fittedprob[,i], lty=fitlty[i], col=fitcol[i], 
				  type=ifelse(discrete, "s", "l"), ...)

    if(lines01)
        abline(h=c(0, 1), lty="dashed", col="grey")
    
    if(addlegend)
    {
        if(missing(legendtext)) 
            legendtext <- paste("fit", 1:nft)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, 
               lty=fitlty, col=fitcol,...)
    }
}
