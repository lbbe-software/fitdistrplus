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
                        datacol, fillrect, fitlty, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
                        ylegend = NULL, lines01 = FALSE,Turnbull.confint = FALSE, NPMLE.method = "Wang", add = FALSE,
                        plotstyle = "graphics", ...)
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
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  if ((plotstyle == "ggplot") & (NPMLE.method == "Turnbull"))
  {
    warning("When NPMLE.method is Turnbull, plotstyle is forced to graphics." )
    plotstyle <- "graphics"
  } 
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("cdfcompcens is not yet available when using weights")
  
  nft <- length(ft)
  if (missing(datacol)) datacol <- "black"
  if (missing(fillrect)) fillrect <- "lightgrey"
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
  
  
  # calculations for Wang method, for both graphics and ggplot displays
  if (NPMLE.method == "Wang") {
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
    
    # Definition of vertices of each rectangle
    Qi.left <- df$left # dim k
    Qi.left4plot <- Qi.left
    if (is.infinite(Qi.left4plot[1]) | is.nan(Qi.left4plot[1])) Qi.left4plot[1] <- xmininf
    Qi.right <- df$right
    Qi.right4plot <- Qi.right
    if (is.infinite(Qi.right4plot[k]) | is.nan(Qi.right4plot[k])) Qi.right4plot[k] <- xmaxinf
    # keep only 16 significants digits for R configured with noLD (--disable-long-double)
    Pi.low <- signif(Fbefore, 16)
    Pi.up <- signif(Fnpsurv, 16)
    
    # the line at right of the rectangles
    dright <- c(f$left[1], rep(f$right, rep(2,k)), f$right[k]) 
    Fright <- rep(c(0,Fnpsurv), rep(2,k+1))
    # the line at left of the rectangles
    dleft <- rep(c(f$left,f$right[k]), rep(2,k+1))
    Fleft <- c(0,rep(Fnpsurv, rep(2,k)),1)
  }
  
  
  # check legend parameters if added
  if(missing(legendtext)) 
  {
    legendtext <- sapply(ft, function(x) x$distname)
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, sapply(ft, function(x) toupper(x$method)), sep="-")
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, 1:nft, sep="-")
  }
  
  
  if(plotstyle == "graphics") {
    ######## plot if plotstyle=='graphics' ########
    
    
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
      # Plot of the ECDF
      if (!add)
        plot(1, 1, type = "n", xlab=xlab, ylab=ylab, main=main, 
             log = logxy, xlim = xlim, ylim = ylim, ...)
      
      # the line at right of the rectangles
      lines(dright, Fright, col = datacol)
      # the line at left of the rectangles
      lines(dleft, Fleft, col = datacol)
      
      # Add of the filled rectangles
      # ca donne un rendu bizarre - plutot ajouter un argument fill.datacol
      # rgbdatacol <- col2rgb(datacol)
      # lightdatacol <- rgb(rgbdatacol[1], rgbdatacol[2], rgbdatacol[3], maxColorValue = 255, 
      #                     alpha = 10)
      for(i in 1:k) {
        rect(xleft = Qi.left4plot, ybottom = Pi.low, xright = Qi.right4plot, ytop = Pi.up,
             border = datacol, col = fillrect)
      }
    }
    
    ################## plot of each fitted distribution
    plot.fti <- function(i, ...)
    {
      fti <- ft[[i]]
      para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
      distname <- fti$distname
      pdistname <- paste("p", distname, sep="")
      if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois")))
        warning(" Be careful, variables are considered continuous in this function!")
      if (xlogscale == TRUE)
      {
        sfin <- 10^seq(log10(xmin), log10(xmax), by=(log10(xmax)-log10(xmin))/100)
      }
      else
      {
        sfin <- seq(xmin, xmax, by=(xmax-xmin)/100)
      }
      theopfin <- do.call(pdistname, c(list(sfin), as.list(para)))
      lines(sfin, theopfin, lty=fitlty[i], col=fitcol[i], ...)
    }
    
    s <- sapply(1:nft, plot.fti, ...)
    rm(s)
    
    if(lines01)
      abline(h=c(0, 1), lty="dashed", col="grey")
    
    if (addlegend)
    {
      legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
    }
    
    invisible()
    
    
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work with plotstyle = 'ggplot'. Please install it", call. = FALSE)
    
  } else {
    ######## plot if plotstyle=='ggplot' ########
    if (NPMLE.method == "Wang") {
      
      # recode the legend position according to available positions in ggplot2
      if(xlegend %in% c("topleft", "bottomleft"))
        xlegend <- "left"
      if(xlegend %in% c("topright", "bottomright"))
        xlegend <- "right"
      if(xlegend == "center")
        xlegend <- "right"
      
      # the line at right of the rectangles
      dsegmright <- cbind(dright, Fright)[2:9,]
      dsegmright <- cbind(dsegmright[-8, ], dsegmright[-1,])
      dsegmright <- as.data.frame(dsegmright)
      colnames(dsegmright) <- c("x1", "y1", "x2", "y2")
      
      # the line at left of the rectangles
      dsegmleft <- cbind(dleft, Fleft)[2:9,]
      dsegmleft <- cbind(dsegmleft[-8, ], dsegmleft[-1,])
      dsegmleft <- as.data.frame(dsegmleft)
      colnames(dsegmleft) <- c("x1", "y1", "x2", "y2")
      
      drect <- data.frame(x1=Qi.left4plot, x2=Qi.right4plot, y1=Pi.low, y2=Pi.up)
      if (xlogscale == TRUE) {
        sfin <- rep(10^seq(log10(xmin), log10(xmax), by=(log10(xmax)-log10(xmin))/100), times = nft)
      } else {
        sfin <- rep(seq(xmin, xmax, by=(xmax-xmin)/100), times = nft)
      }
      theopfin <- vector(mode = "numeric", length = length(sfin))
      ind <- vector(mode = "character", length = length(sfin))
      len <- length(sfin) / nft
      for(i in 1:nft) {
        fti <- ft[[i]]
        para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
        distname <- fti$distname
        if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois")))
          warning(" Be careful, variables are considered continuous in this function!")
        pdistname <- paste("p", distname, sep="")
        theopfin[((i - 1) * len + 1):(i * len)] <- do.call(pdistname, c(list(sfin[((i - 1) * len + 1):(i * len)]), as.list(para)))
        ind[((i - 1) * len + 1):(i * len)] <- distname
      }
      dline <- data.frame(x = sfin, y = theopfin, ind = ind)
      dline$ind <- factor(dline$ind, levels = unique(dline$ind))   # reorder levels in the appearance order of the input
      
      ggcdfcompcens <- ggplot2::ggplot() + 
        ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) + 
        ggplot2::ggtitle(main) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
        {if(lines01) ggplot2::geom_hline(ggplot2::aes(yintercept=0), color="grey", linetype="dashed")} +
        {if(lines01) ggplot2::geom_hline(ggplot2::aes(yintercept=1), color="grey", linetype="dashed")} +
        ggplot2::geom_rect(data=drect, mapping=ggplot2::aes_(xmin=quote(x1), xmax=quote(x2), ymin=quote(y1), ymax=quote(y2)), colour = datacol, fill = fillrect, alpha=0.5) +
        ggplot2::geom_segment(data=dsegmright, mapping=ggplot2::aes_(x=quote(x1), y=quote(y1), xend=quote(x2), yend=quote(y2)), colour = datacol) +
        ggplot2::geom_segment(data=dsegmleft, mapping=ggplot2::aes_(x=quote(x1), y=quote(y1), xend=quote(x2), yend=quote(y2)), colour = datacol) +
        ggplot2::geom_line(data=dline, ggplot2::aes_(quote(x), quote(y), group = quote(ind), colour = quote(ind), linetype = quote(ind))) +
        ggplot2::theme_bw() +
        {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend), plot.title = ggplot2::element_text(hjust = 0.5)) else ggplot2::theme(legend.position = "none", plot.title = ggplot2::element_text(hjust = 0.5))} +
        ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
        ggplot2::scale_linetype_manual(values = fitlty, labels = legendtext) +
        ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
        ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
        {if(xlogscale) ggplot2::scale_x_continuous(trans='log10')} +
        {if(ylogscale) ggplot2::scale_y_continuous(trans='log10')}
      
      return(ggcdfcompcens)
    }
  }
}
