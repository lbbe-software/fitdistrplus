#############################################################################
#   Copyright (c) 2011 Marie Laure Delignette-Muller, Christophe Dutang, Aurélie Siberchicot
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
                    datapch, datacol, fitlty, fitcol, addlegend = TRUE, legendtext, xlegend = "bottomright", 
                    ylegend = NULL, horizontals = TRUE, verticals = FALSE, do.points = TRUE, 
                    use.ppoints = TRUE, a.ppoints = 0.5, lines01 = FALSE, discrete, add = FALSE, plotstyle = "graphics", ...)
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
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("cdfcomp is not yet available when using weights")
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # manage default parameters
  nft <- length(ft)
  if (missing(datapch)) datapch <- 16
  if (missing(datacol)) datacol <- "black"
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitlty)) fitlty <- 1:nft
  fitcol <- rep(fitcol, length.out=nft)
  fitlty <- rep(fitlty, length.out=nft)
  
  if (missing(xlab))
    xlab <- ifelse(xlogscale, "data in log scale", "data")
  if (missing(ylab)) 
    ylab <- "CDF"
  if (missing(main)) 
    main <- paste("Empirical and theoretical CDFs")
  
  # check legend parameters if added
  if(missing(legendtext)) 
    legendtext <- paste("fit", 1:nft)
  
  # initiate discrete if not given 
  if(missing(discrete))
  {
    discrete <- ft[[1]]$discrete
  }
  if(!is.logical(discrete))
    stop("wrong argument 'discrete'.")
  
  
  # check data
  mydata <- ft[[1]]$data
  verif.ftidata <- function(fti)
  {
    if (any(fti$data != mydata))
      stop("All compared fits must have been obtained with the same dataset")
    invisible()
  }
  lapply(ft, verif.ftidata)
  
  # check xlim
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
  
  # some variable definitions
  distname <- ft[[1]]$distname
  n <- length(mydata)
  sdata <- sort(mydata)
  largedata <- (n > 1e4)
  logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
  
  if ((xlogscale == TRUE) & min(mydata) <= 0)
    stop("log transformation of data requires only positive values")
  
  # plot of data (ecdf)
  if(xlogscale && !discrete)
    sfin <- seq(log10(xmin), log10(xmax), by=(log10(xmax)-log10(xmin))/100)
  else # (!xlogscale && !discrete) and discrete
    sfin <- seq(xmin, xmax, length.out=101)
  
  
  # previous version with no vizualisation of ex-aequos
  # obsp <- ecdf(sdata)
  if (use.ppoints && !discrete)
    obsp <- ppoints(n, a = a.ppoints)
  else
    obsp <- (1:n) / n
  
  # computation of each fitted distribution
  comput.fti <- function(i)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    pdistname <- paste("p", distname, sep = "")
    if(xlogscale && !discrete)
    {
      do.call(pdistname, c(list(q=10^sfin), as.list(para)))
    }else
    {
      do.call(pdistname, c(list(q=sfin), as.list(para)))
    }
  }
  fittedprob <- sapply(1:nft, comput.fti)  	
  if(NCOL(fittedprob) != nft || NROW(fittedprob) != length(sfin))
    stop("problem when computing fitted CDFs.")
  
  
  # check ylim
  if(missing(ylim))
    ylim <- range(obsp, fittedprob) 
  else
    ylim <- range(ylim) #in case of users enter a bad ylim
  
  # optional add of horizontal and vertical lines for step function
  xhleft <- sdata[-length(sdata)]
  xhright <- sdata[-1L]
  yh <- obsp[-length(sdata)]
  xv <- xhright
  yvdown <- yh
  yvup <- obsp[-1L] 
  
  if(xlogscale)
    sfin <- 10^sfin
  
  
  if(plotstyle == "graphics") {
    ######## plot if plotstyle=='graphics' ########
    
    #main plot
    if(!add) #create a new graphic
    {
      if(!largedata && do.points)
        plot(sdata, obsp, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
             log=logxy, pch=datapch, col=datacol, type="p", ...)
      else if(largedata)
        plot(sdata, obsp, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
             log=logxy, col=datacol, type="s", ...)
      else if(!do.points)
        plot(sdata, obsp, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
             log=logxy, col=datacol, type="n", ...)     
      else
        stop("internal error in cdfcomp().")
    }else #add to the current graphic
    {
      #do not need parameters: main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, log=logxy,  
      if(!largedata && do.points)
        points(sdata, obsp, pch=datapch, col=datacol, type="p", ...)
      else if(largedata)
        points(sdata, obsp, col=datacol, type="s", ...)
      #else if(!do.points) nothing to plot
      
    }
    
    # optional add of horizontal and vertical lines for step function
    if (!largedata && horizontals) {
      segments(xhleft, yh, xhright, yh, col=datacol,...)
      segments(sdata[length(sdata)], 1, xmax, 1, col=datacol, lty = 2, ...)
      segments(xmin, 0, sdata[1], 0, col=datacol, lty = 2, ...)
      
      if (verticals) {
        segments(xv, yvdown, xv, yvup, col=datacol,...)
        segments(sdata[1], 0, sdata[1], obsp[1], col=datacol, ...)
      }
    }
    
    # plot fitted cdfs
    for(i in 1:nft)
      lines(sfin, fittedprob[,i], lty=fitlty[i], col=fitcol[i], type=ifelse(discrete, "s", "l"), ...)
    
    if(lines01)
      abline(h=c(0, 1), lty="dashed", col="grey")
    
    if(addlegend)
      legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol,...)
    
    invisible()
    
    
  } else if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 needed for this function to work with plotstyle = 'ggplot'. Please install it", call. = FALSE)
    
  } else {
    ######## plot if plotstyle=='ggplot' ########
    
    # recode the legend position according to available positions in ggplot2
    if(xlegend %in% c("topleft", "bottomleft"))
      xlegend <- "left"
    if(xlegend %in% c("topright", "bottomright"))
      xlegend <- "right"
    
    # structure the fittedprob in a relevant data.frame
    fittedprob <- as.data.frame(fittedprob)
    colnames(fittedprob) <- unlist(lapply(ft, function(X) X["distname"]))
    fittedprob <- stack(fittedprob)
    fittedprob$sfin <- sfin   # sfin is recycled in the standard fashion
    fittedprob$ind <- factor(fittedprob$ind, levels = unique(fittedprob$ind))   # reorder levels in the appearance order of the input
    step <- data.frame(values = obsp, ind = "step", sfin = sdata)
    horiz <- data.frame(x = xhleft, y = yh, xend = xhright, yend = yh, ind = "horiz")
    horiz0 <- data.frame(x = xmin, y = 0, xend = sdata[1], yend = 0, ind = "horiz0")
    horiz1 <- data.frame(x = sdata[length(sdata)], y = 1, xend = xmax, yend = 1, ind = "horiz1")
    verti <- data.frame(x = sdata[1], y = 0, xend = sdata[1], yend = obsp[1], ind = "verti")
    
    ggcdfcomp <-
      ggplot2::ggplot(data = fittedprob, ggplot2::aes_(quote(sfin), quote(values), group = quote(ind), colour = quote(ind))) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::ggtitle(main) +
      ggplot2::coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(ylim[1], ylim[2])) +
      
      {if(!largedata && do.points) ggplot2::geom_point(data = step, ggplot2::aes_(quote(sfin), quote(values)), show.legend = FALSE, colour = datacol, shape = datapch)} +
      {if(largedata) ggplot2::geom_step(data = step, ggplot2::aes_(quote(sfin), quote(values)), show.legend = FALSE, colour = datacol, shape = datapch)} +
      
      {if(!largedata && horizontals && !verticals) ggplot2::geom_segment(data = horiz, ggplot2::aes_(x=quote(x), y=quote(y), xend=quote(xend), yend=quote(yend)), show.legend = FALSE, colour = datacol)} +
      {if(!largedata && horizontals && verticals) ggplot2::geom_step(data = step, ggplot2::aes_(quote(sfin), quote(values)), show.legend = FALSE, colour = datacol)} +
      
      {if(!largedata && horizontals) ggplot2::geom_segment(data = horiz1, ggplot2::aes_(x=quote(x), y=quote(y), xend=quote(xend), yend=quote(yend)), show.legend = FALSE, colour = datacol, linetype = 2)} +
      {if(!largedata && horizontals) ggplot2::geom_segment(data = horiz0, ggplot2::aes_(x=quote(x), y=quote(y), xend=quote(xend), yend=quote(yend)), show.legend = FALSE, colour = datacol, linetype = 2)} +
      {if(!largedata && horizontals && verticals) ggplot2::geom_segment(data = verti, ggplot2::aes_(x=quote(x), y=quote(y), xend=quote(xend), yend=quote(yend)), show.legend = FALSE, colour = datacol)} +
      
      {if(discrete) ggplot2::geom_step(data = fittedprob, ggplot2::aes_(linetype = quote(ind), colour = quote(ind)), size = 0.4)} +
      {if(!discrete) ggplot2::geom_line(data = fittedprob, ggplot2::aes_(linetype = quote(ind), colour = quote(ind)), size = 0.4)} +
      
      {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend)) else ggplot2::theme(legend.position = "none")} +
      ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
      ggplot2::scale_linetype_manual(values = fitlty, labels = legendtext) +
      ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
      
      {if(lines01) ggplot2::geom_hline(ggplot2::aes(yintercept=0), color="grey", linetype="dashed")} +
      {if(lines01) ggplot2::geom_hline(ggplot2::aes(yintercept=1), color="grey", linetype="dashed")} +
      {if(xlogscale) ggplot2::scale_x_continuous(trans='log10')} +
      {if(ylogscale) ggplot2::scale_y_continuous(trans='log10')}
    
    return(ggcdfcomp)
  }
}
