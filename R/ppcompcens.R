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
### PP plot for various fits
### of continuous distribution(s) (fitdistcens results)
### on a same dataset
###
###         R functions
###


ppcompcens <- function(ft, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, fillrect,
                   fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", ylegend = NULL, 
                   line01 = TRUE, line01col = "black", line01lty = 1, ynoise = TRUE, ...)
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
    if (nft == 1) fillrect <- "lightpink" else fillrect <- NA

  
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
    xlab <- "Theoretical probabilities"
  if (missing(ylab)) 
    ylab <- "Empirical probabilities"
  if (missing(main)) 
    main <- "P-P plot"
  
  # computation from censdata
  db <- censdata
  db$left[is.na(db$left)] <- -Inf
  db$right[is.na(db$right)] <- Inf
  f <- npsurv(db)$f

  if(missing(xlim) & missing(ylim))
  {
    if (xlogscale == TRUE)
    {
      stop("You must define the limits of axes using xlim or ylim in order 
           to use a logarithmic scale for a P-P plot with censored data.")
    } else
    {
      xlim <- ylim <- c(0, 1)
    }
  } else # at least xlim or ylim are specified
  {
    if (missing(xlim) | missing(ylim))
    {
      warning("By default the same limits are applied to x and y axes.
            You should specifiy both if you want different x and y limits")
      if (missing(xlim)) xlim <- ylim else ylim <- xlim
    }
  }    
 
  k <- length(f$left)
  Fnpsurv <- cumsum(f$p) 
  Fbefore <- c(0, Fnpsurv[-k])
  df <- data.frame(left = f$left, right = f$right)

  # Definition of vertices of each rectangle
  Qi.left <- df$left # dim k
  Qi.right <- df$right
  nQi <- length(Qi.left)
  Pi.low <- Fbefore
  Pi.up <- Fnpsurv
  
  
  ######## plot with graphics ########
  # main plot
  plot(1, 1, type = "n", main = main, xlim = xlim, ylim = ylim,
       xlab = xlab, ylab = ylab, log = logxy)
  
  # plot of rectangles
  plot.fti <- function(i, ...)
  {
    fti <- ft[[i]]
    para=c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    pdistname <- paste("p", distname, sep="")
    
    if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
      warning(" Be careful, variables are considered continuous in this function!")
    Pitheo.low <- do.call(pdistname, c(list(Qi.left), as.list(para)))
    Pitheo.up <- do.call(pdistname, c(list(Qi.right), as.list(para)))
    if (ynoise & nft > 1)
    {
      if (xlogscale == TRUE)
      {
        noise2mult <- runif(nQi, 0.99, 1.01)
        rect(xleft = Pitheo.low, ybottom = Pi.low * noise2mult, 
             xright = Pitheo.up, ytop = Pi.up * noise2mult, 
             border = fitcol[i], col = fillrect)
      }
      else
      {
        noise2add <- runif(nQi, -0.01, 0.01)
        rect(xleft = Pitheo.low, ybottom = Pi.low + noise2add, 
             xright = Pitheo.up, ytop = Pi.up + noise2add, 
             border = fitcol[i], col = fillrect)
      }
    } else # ! ynoise
    {
      rect(xleft = Pitheo.low, ybottom = Pi.low, xright = Pitheo.up, ytop = Pi.up, 
           border = fitcol[i], col = fillrect)
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
}
