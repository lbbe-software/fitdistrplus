#############################################################################
#   Copyright (c) 2012 Christophe Dutang, Aurelie Siberchicot
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
### plot density functions for various fits
### of continuous distribution(s) (fitdist results)
### on a same dataset
###
###         R functions
###


denscomp <- function(ft, xlim, ylim, probability = TRUE, main, xlab, ylab, datacol, fitlty, fitcol, 
                     addlegend = TRUE, legendtext, xlegend = "topright", ylegend = NULL, 
                     demp = FALSE, dempcol = "grey", plotstyle = "graphics", ...)
{
  if(inherits(ft, "fitdist"))
  {
    ft <- list(ft)
  } else if(!is.list(ft))
  {
    stop("argument ft must be a list of 'fitdist' objects")
  } else
  {
    if(any(sapply(ft, function(x) !inherits(x, "fitdist"))))        
      stop("argument ft must be a list of 'fitdist' objects")
  }
  
  # In the future developments, it will be necessary to check that all the fits share the same weights
  if(!is.null(ft[[1]]$weights))
    stop("denscomp is not yet available when using weights")
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # parameters used in 'hist' function 
  argshistPlotFalse <- c("breaks", "nclass", "include.lowest", "right")
  argshistPlotTrue <- c(argshistPlotFalse, "density", "angle", "border", "axes", "labels")
  
  # manage default parameters
  nft <- length(ft)
  if (missing(datacol)) datacol <- NULL    
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitlty)) fitlty <- 1:nft
  fitcol <- rep(fitcol, length.out=nft)
  fitlty <- rep(fitlty, length.out=nft)
  
  if (missing(xlab))
    xlab <- "data"
  if (missing(ylab)) 
    ylab <- ifelse(probability, "Density", "Frequency")
  if (missing(main)) 
    main <- ifelse(probability, "Histogram and theoretical densities", 
                   "Histogram and theoretical frequencies")
  
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
    xlim <- range(mydata)
  }else
  {
    xmin <- xlim[1]
    xmax <- xlim[2]
  }
  
  # some variable definitions
  n <- length(mydata)
  sfin <- seq(xmin, xmax, length.out=101)
  reshist <- hist(mydata, plot = FALSE, ...)
  scalefactor <- ifelse(probability, 1, n * diff(reshist$breaks))
  binwidth <- min(diff(reshist$breaks))
  
  # computation of each fitted distribution
  comput.fti <- function(i)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    ddistname <- paste("d", distname, sep="")
    
    do.call(ddistname, c(list(x=sfin), as.list(para))) * scalefactor
  }
  fitteddens <- sapply(1:nft, comput.fti)
  if(NCOL(fitteddens) != nft || NROW(fitteddens) != length(sfin))
    stop("problem when computing fitted densities.")
  
  # check ylim
  if (missing(ylim))
  {
    if(!probability)
      ylim <- c(0, max(reshist$counts))
    else
      ylim <- c(0, max(reshist$density))
    ylim <- range(ylim, fitteddens)	
  }else
    ylim <- range(ylim) # in case of users enter a bad ylim
  
  # check legend parameters if added
  if (missing(legendtext) && !demp) 
  {
    legendtext <- paste("fit", 1:nft)
  }else if (missing(legendtext) && demp) 
  {
    legendtext <- c(paste("fit", 1:nft), "emp.")
    fitlty <- c(fitlty, 1)
    fitcol <- c(fitcol, dempcol)
  }else if(demp)
  {
    legendtext <- c(legendtext, "emp.")
    fitlty <- c(fitlty, 1)
    fitcol <- c(fitcol, dempcol)
  }   
  
  if(plotstyle == "graphics") {
    ######## plot if plotstyle=='graphics' ########
    
    #main plotting
    reshist <- hist(mydata, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, col = datacol, probability = probability, ...)
    
    #plot fitted densities
    for(i in 1:nft)
      lines(sfin, fitteddens[,i], lty=fitlty[i], col=fitcol[i], ...)
    
    #plot empirical density
    if(demp)
      lines(density(mydata)$x, density(mydata)$y * scalefactor, col=dempcol)
    
    if (addlegend)
      legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
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
    
    # the default colors of the bars is the same as panel.background.fill in theme_grey()
    if(is.null(datacol))
      datacol <- "grey92"
    
    # structure the fitteddens in a relevant data.frame
    fitteddens <- as.data.frame(fitteddens)
    colnames(fitteddens) <- unlist(lapply(ft, function(X) X["distname"]))
    fitteddens <- stack(fitteddens)
    fitteddens$sfin <- sfin   # sfin is recycled in the standard fashion
    fitteddens$ind <- factor(fitteddens$ind, levels = unique(fitteddens$ind))   # reorder levels in the appearance order of the input
    if(demp) # bind empirical data if demp is TRUE
      fitteddens <- rbind(fitteddens, data.frame(values = density(mydata)$y * scalefactor, ind = "demp", sfin = density(mydata)$x))
    
    histdata <- data.frame(values = mydata, ind = "hist", sfin = mydata) # the added data must have the same column names as the main data to be compatible with ggplot
    binwidth <- min(diff(reshist$breaks))
    
    ggdenscomp <-
      ggplot2::ggplot(fitteddens, ggplot2::aes_(quote(sfin), quote(values), group = quote(ind), colour = quote(ind))) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::ggtitle(main) +
      ggplot2::coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(ylim[1], ylim[2])) +
      {if(probability) ggplot2::geom_histogram(data = histdata, ggplot2::aes_(quote(values), quote(..density..)), binwidth = binwidth, boundary = 0, show.legend = FALSE, col = "black", alpha = 1, fill = datacol)
        else ggplot2::geom_histogram(data = histdata, ggplot2::aes_(quote(values), quote(..count..)), binwidth = binwidth, boundary = 0, show.legend = FALSE, col = "black", alpha = 1, fill = datacol)} +
      ggplot2::geom_line(data = fitteddens, ggplot2::aes_(linetype = quote(ind), colour = quote(ind)), size = 0.4) +
      ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
      {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend)) else ggplot2::theme(legend.position = "none")} +
      ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
      ggplot2::scale_linetype_manual(values = fitlty, labels = legendtext)
    return(ggdenscomp)
  }
}
