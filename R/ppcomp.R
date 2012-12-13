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
### P-P plot for various fits
### of continuous distribution(s) (fitdist results)
### on a same dataset
###
###         R functions
###



ppcomp <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
    fitpch, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", ylegend = NULL, 
    use.ppoints = TRUE, a.ppoints = 0.5, line01 = TRUE, line01col = "black", line01lty = 1,
    ynoise = TRUE, ...)
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
    if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (missing(fitpch)) fitpch <- 21
    fitcol <- rep(fitcol, length.out=nft)
    fitpch <- rep(fitpch, length.out=nft)
    
    if (missing(xlab))
        xlab <- "Theoretical probabilities"
    if (missing(ylab)) 
        ylab <- "Empirical probabilities"
    if (missing(main)) 
        main <- "P-P plot"
    
    mydata <- ft[[1]]$data

    verif.ftidata <- function(fti)
    {
        if (any(fti$data != mydata))
            stop("All compared fits must have been obtained with the same dataset")
        invisible()
    }
    lapply(ft, verif.ftidata)

    n <- length(mydata)
    sdata <- sort(mydata)
    if (use.ppoints)
        obsp <- ppoints(n, a = a.ppoints)
    else
        obsp <- (1:n) / n
    
    # computation of each fitted distribution
    comput.fti <- function(i, ...)
    {
        fti <- ft[[i]]
        para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
        distname <- fti$distname
        pdistname <- paste("p", distname, sep="")
        do.call(pdistname, c(list(q=sdata), as.list(para)))
    }
    fittedprob <- sapply(1:nft, comput.fti, ...)
    
    if (missing(xlim))
        xlim <- range(fittedprob)
    if (missing(ylim))
        ylim <- range(obsp)

    logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
    #main plotting
    resquant <- plot(fittedprob[,1], obsp, main=main, xlab=xlab, ylab=ylab, log=logxy,
            pch=fitpch[1], xlim=xlim, ylim=ylim, col=fitcol[1], ...)
    
    #plot fitted quantiles
    if(nft > 1 && !ynoise)
        for(i in 2:nft)
            points(fittedprob[,i], obsp, pch=fitpch[i], col=fitcol[i], ...)
    if(nft > 1 && ynoise)
        for(i in 2:nft)
            points(fittedprob[,i], obsp*(1 + rnorm(n, 0, 0.01)), pch=fitpch[i], col=fitcol[i], ...)
    
    if(line01)
        abline(0, 1, lty=line01lty, col=line01col)
    
    if(addlegend)
    {
        if (missing(legendtext)) 
            legendtext <- paste("fit",1:nft)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, pch=fitpch, col=fitcol, ...)
    }
}
