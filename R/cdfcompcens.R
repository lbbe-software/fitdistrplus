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

cdfcompcens <-
function(ft,xlogscale=FALSE,addlegend=TRUE,legendtext,datacol,fitcol,fitlty,xlab,ylab,
main,xlegend = "bottomright",ylegend = NULL, ...){


    if(inherits(ft, "fitdistcens"))
            ft <- list(ft)
    nft <- length(ft)
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
        if (xlogscale == TRUE) xlab <- "censored data in log scale"
        else xlab <- "censored data"
    }
    if (missing(ylab)) ylab <- "CDF"
    if (missing(main)) main <- paste("Empirical and theoretical CDFs")

# verification of the content of the list ft
    verif.fti <- function(fti)
    {
        if (!inherits(fti, "fitdistcens"))
        stop("argument ft must be a list of 'fitdistcens' objects")
    }


    if (is.list(ft))
    {
        l<-lapply( ft,verif.fti)
        rm(l)
    }
    else
    {
        stop("argument ft must be a list of 'fitdistcens' objects")
    }

    censdata <- ft[[1]]$censdata
    xmin <- min(censdata$left,na.rm=TRUE) 
    xmax <- max(censdata$right,na.rm=TRUE)
    if ((xlogscale == TRUE) & xmin <= 0) 
        stop("log transformation of data requires only positive values")


    verif.ftidata <- function(fti)
    {
        if (any(fti$censdata$left != censdata$left,na.rm=TRUE) | 
            any(fti$censdata$right != censdata$right,na.rm=TRUE))
            stop("All compared fits must have been obtained with the same dataset")
    }
    l<-lapply( ft,verif.ftidata)
    rm(l)

    # plot of data (ecdf) using Turnbull algorithm
    survdata <- Surv(time = censdata$left, time2 = censdata$right, type="interval2")
    survfitted <- survfit(survdata ~ 1)
    if (xlogscale == TRUE)
        plot(survfitted,fun="event",xlab=xlab,ylab=ylab,main=main,log="x",col=datacol,...)
    else
        plot(survfitted,fun="event",xlab=xlab,ylab=ylab,main=main,col=datacol,...)


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
