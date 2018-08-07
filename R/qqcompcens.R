#############################################################################
#   Copyright (c) 2018 Marie Laure Delignette-Muller, Christophe Dutang, 
#                      Aurelie Siberchicot
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
### QQ plot for various fits
### of continuous distribution(s) (fitdistcens results)
### on a same dataset
###
###         R functions
###


qqcompcens <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, fillrect,
                       fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", ylegend = NULL, 
                       line01 = TRUE, line01col = "black", line01lty = 1, ynoise = TRUE, plotstyle = "graphics", ...)
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
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("qqcompcens is not yet available when using weights")
  
  censdata <- ft[[1]]$censdata
  
  # check data
  verif.ftidata <- function(fti)
  {
    if (any(fti$censdata$left != censdata$left, na.rm=TRUE) | 
        any(fti$censdata$right != censdata$right, na.rm=TRUE))
      stop("All compared fits must have been obtained with the same dataset")
  }
  l <- lapply( ft, verif.ftidata)
  rm(l)
  
  if (xlogscale != ylogscale)
  {
    xlogscale <- ylogscale <- TRUE
    warning("As a Q-Q plot should use the same scale on x and y axes, 
            both axes were put in a logarithmic scale.")
  }
  logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
  
  # manage default parameters
  nft <- length(ft)
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  fitcol <- rep(fitcol, length.out=nft)
  if (missing(fillrect)) 
    if ((nft == 1) | plotstyle == "ggplot") fillrect <- "lightgrey" else fillrect <- NA
  
  
  # check legend parameters if added
  if(missing(legendtext)) 
  {
    legendtext <- sapply(ft, function(x) x$distname)
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, sapply(ft, function(x) toupper(x$method)), sep="-")
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, 1:nft, sep="-")
  }
  
  if (missing(xlab))
    xlab <- "Theoretical quantiles"
  if (missing(ylab)) 
    ylab <- "Empirical quantiles"
  if (missing(main)) 
    main <- "Q-Q plot"
  
  # computation from censdata
  db <- censdata
  db$left[is.na(db$left)] <- -Inf
  db$right[is.na(db$right)] <- Inf
  f <- npsurv(db)$f
  bounds <- c(f$right, f$left)
  finitebounds <- bounds[is.finite(bounds)]
  
  if(missing(xlim) & missing(ylim))
  {
    upper <- max(finitebounds)
    lower <- min(finitebounds)
    width <- upper - lower
    if (xlogscale == TRUE)
    {
      xmin <- lower * (upper / lower)^(-0.1)
      xmax <- upper * (upper / lower)^0.1
      xmininf <- lower * (upper / lower)^(-10) # 10 to be very large
      xmaxinf <- upper * (upper / lower)^10
    } else
    {
      xmin <- lower - width * 0.1
      xmax <- upper + width * 0.1
      xmininf <- lower - width * 10
      xmaxinf <- upper + width * 10
    }
    xlim <- c(xmin, xmax)
    ylim <- c(xmin, xmax)
  } else # at least xlim or ylim are specified
  {
    if (missing(xlim) | missing(ylim))
    {
      warning("By default the same limits are applied to x and y axes.
            You should specify both if you want different x and y limits")
      if (missing(xlim)) xlim <- ylim else ylim <- xlim
    }
    lower <- min(c(xlim, ylim))
    upper <- max(c(xlim, ylim))
    width <- upper - lower
    if (xlogscale == TRUE)
    {
      xmininf <- lower * (upper / lower)^(-10) # 10 to be very large
      xmaxinf <- upper * (upper / lower)^10
    } else
    {
      xmininf <- lower - width * 10
      xmaxinf <- upper + width * 10
    }
  }    
  
  k <- length(f$left)
  Fnpsurv <- cumsum(f$p) 
  Fbefore <- c(0, Fnpsurv[-k])
  df <- data.frame(left = f$left, right = f$right)
  
  # Definition of vertices of each rectangle
  Qi.left <- df$left # dim k
  Qi.left4plot <- Qi.left
  if (Qi.left4plot[1] == - Inf) Qi.left4plot[1] <- xmininf
  Qi.right <- df$right
  Qi.right4plot <- Qi.right
  if (Qi.right4plot[k] == Inf) Qi.right4plot[k] <- xmaxinf
  Pi.low <- Fbefore
  Pi.up <- Fnpsurv
  nPi <- length(Pi.low)
  
  lrect <- vector(mode = "list", length = nft)
  for(i in 1:nft)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    qdistname <- paste("q", distname, sep="")
    
    if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois")))
      warning(" Be careful, variables are considered continuous in this function!")
    
    Qitheo.left <- do.call(qdistname, c(list(Pi.low), as.list(para)))
    Qitheo.right <- do.call(qdistname, c(list(Pi.up), as.list(para)))
    Qitheo.left4plot <- Qitheo.left
    if (Qitheo.left4plot[1] == - Inf) Qitheo.left4plot[1] <- xmininf
    Qitheo.right4plot <- Qitheo.right
    if (Qitheo.right4plot[k] == Inf) Qitheo.right4plot[k] <- xmaxinf
    lrect[[i]] <- data.frame(Qitheo.left4plot = Qitheo.left4plot, 
                             Qi.left4plot = Qi.left4plot, 
                             Qitheo.right4plot = Qitheo.right4plot, 
                             Qi.right4plot = Qi.right4plot, ind = legendtext[i])
  }
  
  if(plotstyle == "graphics") {
    ######## plot if plotstyle=='graphics' ########
    
    # main plot
    plot(1, 1, type = "n", main = main, xlim = xlim, ylim = ylim,
         xlab = xlab, ylab = ylab, log = logxy)
    
    # plot of rectangles
    plot.fti <- function(i, ...)
    {
      Qitheo.left4plot <- lrect[[i]]$Qitheo.left4plot 
      Qi.left4plot <- lrect[[i]]$Qi.left4plot 
      Qitheo.right4plot <- lrect[[i]]$Qitheo.right4plot 
      Qi.right4plot <- lrect[[i]]$Qi.right4plot
      
      if (ynoise & nft > 1)
      {
        if (xlogscale == TRUE)
        {
          noise2mult <- runif(nPi, 0.99, 1.01)
          rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot * noise2mult, 
               xright = Qitheo.right4plot, 
               ytop = Qi.right4plot * noise2mult, 
               border = fitcol[i], col = fillrect[i])
        }
        else
        {
          noise2add <- runif(nPi, -width*0.01, width*0.01)
          rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot + noise2add, 
               xright = Qitheo.right4plot, 
               ytop = Qi.right4plot + noise2add, 
               border = fitcol[i], col = fillrect[i])
        }
      } else # ! ynoise
      {
        rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot, xright = Qitheo.right4plot, 
             ytop = Qi.right4plot, 
             border = fitcol[i], col = fillrect[i])
      }
    }
    s <- sapply(1:nft, plot.fti, ...)
    rm(s)
    
    if(line01)
      abline(0, 1, lty = line01lty, col = line01col)
    
    if (addlegend)
    {
      legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, col=fitcol, lty = 1, ...)
    }
    invisible()
    
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work with plotstyle = 'ggplot'. Please install it", call. = FALSE)
    
  } else {
    ######## plot if plotstyle=='ggplot' ########
    
    drect <-  do.call("rbind", lrect)
    ind <- as.factor(drect$ind)
    fitcol <- rep(fitcol, table(ind))
    fillrect <- if(length(fillrect) > 1) {rep(fillrect, table(ind))} else {fillrect}
    
    ggqqcompcens <- ggplot2::ggplot(drect) + 
      ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)  +
      ggplot2::ggtitle(main) + ggplot2::xlab(xlab) + ggplot2::ylab(ylab) +
      ggplot2::geom_rect(data=drect, mapping=ggplot2::aes_(xmin=quote(Qitheo.left4plot), xmax=quote(Qitheo.right4plot), ymin=quote(Qi.left4plot), ymax=quote(Qi.right4plot)), colour = fitcol, fill = fillrect, alpha=0.5) +
      ggplot2::theme_bw() +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
      {if(line01) ggplot2::geom_abline(ggplot2::aes(slope = 1, intercept = 0), color = line01col, linetype = line01lty)} +
      {if(xlogscale) ggplot2::scale_x_continuous(trans='log10')} +
      {if(ylogscale) ggplot2::scale_y_continuous(trans='log10')} + 
      ggplot2::facet_wrap(~ind)
    
    return(ggqqcompcens)
  }
}
