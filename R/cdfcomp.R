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

cdfcomp <-
function(ft,xlogscale=FALSE,addlegend=TRUE,legendtext,datapch,datacol,fitcol,fitlty,xlab,ylab,xlim,
main,xlegend = "bottomright",ylegend = NULL,horizontals = TRUE,verticals = FALSE, ...){


    if(inherits(ft, "fitdist"))
            ft <- list(ft)
    nft <- length(ft)
    if (missing(datapch)) datapch <- 16
    if (missing(datacol)) datacol <- "black"
    if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (length(fitcol)!=nft)
        stop("if specified, fitcol must be a vector of length
        the number of fitted distributions to represent")
    if (missing(fitlty)) fitlty <- 1:nft
    if (length(fitlty)!=nft)
        stop("if specified, fitlty must be a vector of length
        the number of fitted distributions to represent")
    if (missing(xlab))
    {
        if (xlogscale == TRUE) xlab <- "data in log scale"
        else xlab <- "data"
    }
    if (missing(ylab)) ylab <- "CDF"
    if (missing(main)) main <- paste("Empirical and theoretical CDFs")

# verification of the content of the list ft
    verif.fti <- function(fti)
    {
        if (!inherits(fti, "fitdist"))
        stop("argument ft must be a list of 'fitdist' objects")
    }


    if (is.list(ft))
    {
        l<-lapply( ft,verif.fti)
        rm(l)
    }
    else
    {
        stop("argument ft must be a list of 'fitdist' objects")
    }

    data <- ft[[1]]$data
    if ((xlogscale == TRUE) & min(data)<=0)
        stop("log transformation of data requires only positive
    values")

    if (missing(xlim))
    {
        xmin <- min(data)
        xmax <- max(data)
        xlim <- c(xmin,xmax)
    }
    else
    {
        xmin <- xlim[1]
        xmax <- xlim[2]
    }

    verif.ftidata <- function(fti)
    {
        if (any(fti$data != data))
            stop("All compared fits must have been obtained with the same dataset")
    }
    l<-lapply( ft,verif.ftidata)
    rm(l)

    # plot of data (ecdf)
    n <- length(data)
    s <- sort(data)
    obsp <- ecdf(s)(s)
    if (xlogscale == TRUE)
        plot(s,obsp,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=c(0,1),log="x",pch=datapch,col=datacol,...)
    else
        plot(s,obsp,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=c(0,1),pch=datapch,col=datacol,...)
    # optional add of horizontal and vbertical lines for step function
    if (horizontals)
    {
        xhleft <- s[-length(s)]
        xhright <- s[-1L]
        yh <- obsp[-length(s)]
        segments(xhleft,yh,xhright,yh,col=datacol,...)
        segments(s[length(s)],1,xmax,1,col=datacol,lty = 2, ...)
        
        segments(xmin,0,s[1],0,col=datacol,lty = 2, ...)
        if (verticals)
        {
           xv <-xhright
           yvdown <- yh
           yvup <- obsp[-1L] 
           segments(xv,yvdown,xv,yvup,col=datacol,...)
           segments(s[1],0,s[1],obsp[1],col=datacol, ...)
        }
    }

    # plot of each fitted distribution
    plot.fti <- function(i,...)
    {
        fti <- ft[[i]]
        para=c(as.list(fti$estimate),as.list(fti$fix.arg))
        distname <- fti$distname
        pdistname<-paste("p",distname,sep="")
        if (is.element(distname,c("binom","nbinom","geom","hyper","pois")))
            warning(" Be careful, variables are considered continuous in this function!")
        if (xlogscale == TRUE)
        {
            sfin<-seq(log10(xmin),log10(xmax),by=(log10(xmax)-log10(xmin))/100)
            theopfin<-do.call(pdistname,c(list(q=10^sfin),as.list(para)))
            lines(10^sfin,theopfin,lty=fitlty[i],col=fitcol[i],...)
        }
        else
        {
            sfin<-seq(xmin,xmax,by=(xmax-xmin)/100)
            theopfin<-do.call(pdistname,c(list(q=sfin),as.list(para)))
            lines(sfin,theopfin,lty=fitlty[i],col=fitcol[i],...)
        }
    }
    s<-sapply(1:nft,plot.fti,...)
    rm(s)

    if (addlegend==TRUE)
    {
        if (missing(legendtext)) legendtext <- paste("fit",1:nft)
#       next lines replaced by default argument xlegend fixed to "bottomright"
#        if (missing(xlegend))
#        {
#            if ((xlogscale == TRUE))
#            {
#                xlegendlog10 <- log10(xmin) + (log10(xmax) - log10(xmin))*2/3
#                xlegend <- 10^xlegendlog10
#            }
#            else
#                xlegend <- xmin + (xmax - xmin)*2/3
#        }
#        if (missing(ylegend)) ylegend <- 0.5

        legend(x=xlegend,y=ylegend,bty="n",legend=legendtext,lty=fitlty,col=fitcol,...)
    }
}
