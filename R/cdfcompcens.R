#############################################################################
#   Copyright (c) 2011 Marie Laure Delignette-Muller
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

cdfcompcens <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    datacol, fitlty, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
    ylegend = NULL, lines01 = FALSE,Turnbull.confint = FALSE,...)
{

  if(inherits(ft, "fitdistcens"))
  {
    ft <- list(ft)
  }else if(!is.list(ft))
  {
    stop("argument ft must be a list of 'fitdistcens' objects")
  }else
  {
    if(any(sapply(ft, function(x) !inherits(x, "fitdistcens"))))        
      stop("argument ft must be a list of 'fitdistcens' objects")
  }
  
    nft <- length(ft)
    if (missing(datacol)) datacol <- "black"
    if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (missing(fitlty)) fitlty <- 1:nft
    fitcol <- rep(fitcol, length.out=nft)
    fitlty <- rep(fitlty, length.out=nft)

    if (missing(xlab))
        xlab <- ifelse(xlogscale, "censored data in log scale", "censored data")
    if (missing(ylab)) ylab <- "CDF"
    if (missing(main)) main <- paste("Empirical and theoretical CDFs")


    censdata <- ft[[1]]$censdata
    if(missing(xlim))
    {
        xmin <- min(c(censdata$left, censdata$right), na.rm=TRUE)
        xmax <- max(c(censdata$left, censdata$right), na.rm=TRUE)
        xlim <- c(xmin, xmax)
    }
    else
    {
        xmin <- xlim[1]
        xmax <- xlim[2]
    }
    if ((xlogscale == TRUE) & xmin <= 0) 
        stop("log transformation of data requires only positive values")


    verif.ftidata <- function(fti)
    {
        if (any(fti$censdata$left != censdata$left, na.rm=TRUE) | 
            any(fti$censdata$right != censdata$right, na.rm=TRUE))
            stop("All compared fits must have been obtained with the same dataset")
    }
    l <- lapply( ft, verif.ftidata)
    rm(l)

    # plot of data (ecdf) using Turnbull algorithm
    survdata <- Surv(time = censdata$left, time2 = censdata$right, type="interval2")
    survfitted <- survfit(survdata ~ 1)
    
    logxy <- paste(ifelse(xlogscale, "x", ""), ifelse(ylogscale, "y", ""), sep="")
    #main plotting
    if(missing(ylim))
    {
        if (Turnbull.confint)
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                log=logxy, col=datacol, xlim = xlim, ...)
        else
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                log=logxy, col=datacol, conf.int = FALSE, xlim = xlim, ...)
    }
    else
    {
        if (Turnbull.confint)
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                log=logxy, col=datacol, xlim = xlim, ylim=ylim, ...)
        else
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                log=logxy, col=datacol, conf.int = FALSE, xlim = xlim, ylim = ylim, ...)
    }
    # plot of each fitted distribution
    plot.fti <- function(i, ...)
    {
        fti <- ft[[i]]
        para=c(as.list(fti$estimate), as.list(fti$fix.arg))
        distname <- fti$distname
        pdistname <- paste("p", distname, sep="")
        if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois")))
            warning(" Be careful, variables are considered continuous in this function!")
        if (xlogscale == TRUE)
        {
            sfin <- seq(log10(xmin), log10(xmax), by=(log10(xmax)-log10(xmin))/100)
            theopfin <- do.call(pdistname, c(list(q=10^sfin), as.list(para)))
            lines(10^sfin, theopfin, lty=fitlty[i], col=fitcol[i], ...)
        }
        else
        {
            sfin <- seq(xmin, xmax, by=(xmax-xmin)/100)
            theopfin <- do.call(pdistname, c(list(q=sfin), as.list(para)))
            lines(sfin, theopfin, lty=fitlty[i], col=fitcol[i], ...)
        }
    }
    s <- sapply(1:nft, plot.fti, ...)
    rm(s)

    if(lines01)
        abline(h=c(0, 1), lty="dashed", col="grey")

    if (addlegend)
    {
        if (missing(legendtext)) 
            legendtext <- paste("fit", 1:nft)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
    }
}
