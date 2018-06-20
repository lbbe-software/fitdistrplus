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
    datacol, filldatacol, fitlty, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
    ylegend = NULL, lines01 = FALSE,Turnbull.confint = FALSE, NPMLE.method = "Wang", add = FALSE,...)
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
  
  if ((Turnbull.confint == TRUE) & (NPMLE.method == "Wang"))
  {
    warning("When Turnbull.confint is TRUE NPMLE.method is forced to Turnbull." )
    NPMLE.method <- "Turnbull"
  } 
  
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("cdfcompcens is not yet available when using weights")
  
    nft <- length(ft)
    if (missing(datacol)) datacol <- "black"
    if (missing(filldatacol)) filldatacol <- "lightgrey"
    if (missing(fitcol)) fitcol <- 2:(nft+1)
    if (missing(fitlty)) fitlty <- 1:nft
    fitcol <- rep(fitcol, length.out=nft)
    fitlty <- rep(fitlty, length.out=nft)

    if (missing(xlab))
        xlab <- ifelse(xlogscale, "censored data in log scale", "censored data")
    if (missing(ylab)) ylab <- "CDF"
    if (missing(main)) main <- paste("Empirical and theoretical CDFs")

    censdata <- ft[[1]]$censdata
    logxy <- paste(ifelse(xlogscale, "x", ""), ifelse(ylogscale, "y", ""), sep="")
    
    verif.ftidata <- function(fti)
    {
        if (any(fti$censdata$left != censdata$left, na.rm=TRUE) | 
            any(fti$censdata$right != censdata$right, na.rm=TRUE))
            stop("All compared fits must have been obtained with the same dataset")
    }
    l <- lapply( ft, verif.ftidata)
    rm(l)
  ####### Plot of the data ############################
  if (NPMLE.method == "Turnbull")
  # Turnbull plot
  {
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
    
    # plot of data (ecdf) using Turnbull algorithm
    survdata <- Surv(time = censdata$left, time2 = censdata$right, type="interval2")
    survfitted <- survfit(survdata ~ 1)
    
    #main plotting
    if(missing(ylim))
    {
        if (Turnbull.confint)
        {
          if (!add)
          plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
               log=logxy, col=datacol, xlim = xlim, ...)
          else
          lines(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                 log=logxy, col=datacol, xlim = xlim, ...)
          
        }else
        {
          if (!add)
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
               log=logxy, col=datacol, conf.int = FALSE, xlim = xlim, ...)
          else
            lines(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                 log=logxy, col=datacol, conf.int = FALSE, xlim = xlim, ...)
          
        }
    }
    else
    {
        if (Turnbull.confint)
        {
          if (!add)
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
               log=logxy, col=datacol, xlim = xlim, ylim=ylim, ...)
          else
            lines(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                 log=logxy, col=datacol, xlim = xlim, ylim=ylim, ...)
          
        } else
        {
          if (!add)
            plot(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
               log=logxy, col=datacol, conf.int = FALSE, xlim = xlim, ylim = ylim, ...)
          else
            lines(survfitted, fun="event", xlab=xlab, ylab=ylab, main=main, 
                 log=logxy, col=datacol, conf.int = FALSE, xlim = xlim, ylim = ylim, ...)
          
        }
    }
  } else # if NPMLE.method == "Wang"
  # Wang plot
  {
    db <- censdata
    db$left[is.na(db$left)] <- -Inf
    db$right[is.na(db$right)] <- Inf
    f <- npsurv(db)$f
    bounds <- c(f$right, f$left)
    finitebounds <- bounds[is.finite(bounds)]
    upper <- max(finitebounds)
    lower <- min(finitebounds)
    width <- upper - lower
    
    if(missing(xlim))
    {
      if (xlogscale == TRUE)
      {
        xmin <- lower * (upper / lower)^(-0.1)
        xmax <- upper * (upper / lower)^0.1
      } else
      {
        xmin <- lower - width * 0.1
        xmax <- upper + width * 0.1
      }
      xlim <- c(xmin, xmax)
    }
    else
    {
      xmin <- xlim[1]
      xmax <- xlim[2]
    }
    if(missing(ylim))
    {
      ylim <- c(0,1)
    }
    
    if ((xlogscale == TRUE) & xmin <= 0) 
      stop("log transformation of data requires only positive values")
    
    if (xlogscale == TRUE)
    {
      xmininf <- lower * (upper / lower)^(-10) # 10 to be very large
      xmaxinf <- upper * (upper / lower)^10
    } else
    {
      xmininf <- lower - width * 10
      xmaxinf <- upper + width * 10
    }
    k <- length(f$left)
    Fnpsurv <- cumsum(f$p) 
    
    ## calculation of points for Q and P in graphs
    Fbefore <- c(0, Fnpsurv[-k])
    df <- data.frame(left = f$left, right = f$right)
    dfb <- df[(df$left != -Inf) & (df$right != Inf), ]
    
    # Definition of vertices of each rectangle
    Qi.left <- df$left # dim k
    Qi.left4plot <- Qi.left
    if (Qi.left4plot[1] == - Inf) Qi.left4plot[1] <- xmininf
    Qi.right <- df$right
    Qi.right4plot <- Qi.right
    if (Qi.right4plot[k] == Inf) Qi.right4plot[k] <- xmaxinf
    Pi.low <- Fbefore
    Pi.up <- Fnpsurv
    
    # Plot of the ECDF
    if (!add)
    plot(1, 1, type = "n", xlab=xlab, ylab=ylab, main=main, 
         log = logxy, xlim = xlim, ylim = ylim, ...)
    
    # the line at right of the rectangles
    dright <- c(f$left[1], rep(f$right, rep(2,k)), f$right[k]) 
    Fright <- rep(c(0,Fnpsurv), rep(2,k+1))
    lines(dright, Fright, col = datacol)
    ### the line at left of the rectangles
    dleft = rep(c(f$left,f$right[k]), rep(2,k+1))
    Fleft = c(0,rep(Fnpsurv, rep(2,k)),1)
    lines(dleft, Fleft, col = datacol)
    
    # Add of the filled rectangles
    # ca donne un rendu bizarre - plutot ajouter un argument fill.datacol
    # rgbdatacol <- col2rgb(datacol)
    # lightdatacol <- rgb(rgbdatacol[1], rgbdatacol[2], rgbdatacol[3], maxColorValue = 255, 
    #                     alpha = 10)
    for(i in 1:k) {
      rect(xleft = Qi.left4plot, ybottom = Pi.low, xright = Qi.right4plot, ytop = Pi.up,
           border = datacol, col = filldatacol)
    }
  }
  
  ################## plot of each fitted distribution
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
      # check legend parameters if added
      if(missing(legendtext)) 
      {
        legendtext <- sapply(ft, function(x) x$distname)
        if(length(legendtext) != length(unique(legendtext)))
          legendtext <- paste(legendtext, sapply(ft, function(x) toupper(x$method)), sep="-")
        if(length(legendtext) != length(unique(legendtext)))
          legendtext <- paste(legendtext, 1:nft, sep="-")
      }
      legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
    }
    invisible()
}
