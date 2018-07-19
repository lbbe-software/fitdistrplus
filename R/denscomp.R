#############################################################################
#   Copyright (c) 2012 Christophe Dutang, Aurelie Siberchicot, 
#                      Marie Laure Delignette-Muller
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


denscomp <- function(ft, xlim, ylim, probability = TRUE, main, xlab, ylab, 
                     datacol, fitlty, fitcol, addlegend = TRUE, legendtext, 
                     xlegend = "topright", ylegend = NULL, demp = FALSE, 
                     dempcol = "black", plotstyle = "graphics", 
                     discrete, fitnbpts = 101, fittype="l", ...)
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
    stop("denscomp is not yet available when using weights")
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # parameters used in 'hist' function 
  ###### Where are they used - to remove ? !!!!!!!!!!!!!!!!!!!!!!!!
  argshistPlotFalse <- c("breaks", "nclass", "include.lowest", "right")
  argshistPlotTrue <- c(argshistPlotFalse, "density", "angle", "border", "axes", "labels")
  
  # manage default parameters
  nft <- length(ft)
  if (missing(datacol)) datacol <- NULL    
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitlty)) fitlty <- 1:nft
  fitcol <- rep(fitcol, length.out=nft)
  fitlty <- rep(fitlty, length.out=nft)
  fittype <- match.arg(fittype[1], c("p", "l", "o"))
  
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
  # initiate discrete if not given 
  if(missing(discrete))
  {
    discrete <- any(sapply(ft, function(x) x$discrete))
  }
  if(!is.logical(discrete))
    stop("wrong argument 'discrete'.")
  if(!is.logical(demp))
    stop("wrong argument 'discrete'.")
  
  # some variable definitions
  n <- length(mydata)
  if(!discrete)
    sfin <- seq(xmin, xmax, length.out = fitnbpts[1])
  else
    sfin <- unique(round(seq(xmin, xmax, length.out = fitnbpts[1]), digits = 0))
  reshist <- hist(mydata, plot = FALSE, ...)
  if (!discrete)
  {
    if (probability)
    {
      scalefactor <- 1
    } else
    {
      if (length(unique(diff(reshist$breaks))) > 1) # wrong histogram and not possibleto compute a scalefactor
        stop("You should not use probability = FALSE with non-equidistant breaks for the histogram !") else
      scalefactor <- n * diff(reshist$breaks)[1]
    }
#    previous writing that gave incorrect output in case of probability = 1 and non-equidistant breaks
#    scalefactor <- ifelse(probability, 1, n * diff(reshist$breaks))
  } else
  {
    scalefactor <- ifelse(probability, 1, n)
  }
#  binwidth <- min(diff(reshist$breaks))
  
  # computation of each fitted distribution
  comput.fti <- function(i)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    ddistname <- paste("d", distname, sep="")
    
    do.call(ddistname, c(list(sfin), as.list(para))) * scalefactor
  }
  fitteddens <- sapply(1:nft, comput.fti)
  if(NCOL(fitteddens) != nft || NROW(fitteddens) != length(sfin))
    stop("problem when computing fitted densities.")
  
  # check ylim
  if (missing(ylim))
  {
    if(!probability)
      if (discrete) 
      {
        ylim <- c(0, max(as.numeric(table(mydata))))
      } else
      {
        ylim <- c(0, max(reshist$counts))
      }
    else # so if probability
    {
      if (discrete) 
      {
        ylim <- c(0, max(as.numeric(table(mydata))/length(mydata)))
      } else
      {
        ylim <- c(0, max(reshist$density))
      }
    }
    ylim <- range(ylim, fitteddens)	
  }else
    ylim <- range(ylim) # in case of users enter a bad ylim
  
  # check legend parameters if added
  if(missing(legendtext)) 
  {
    legendtext <- sapply(ft, function(x) x$distname)
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, sapply(ft, function(x) toupper(x$method)), sep="-")
    if(length(legendtext) != length(unique(legendtext)))
      legendtext <- paste(legendtext, 1:nft, sep="-")
  }
  
  # forces demp to TRUE if discrete is TRUE
  if(discrete)
    demp <- TRUE
  
  #add empirical density/fmp to legend vectors
  if(demp)
  {
    legendtext <- c(legendtext, "emp.")
    fitlty <- c(fitlty, 1)
    fitcol <- c(fitcol, dempcol)
  }
  
  if(plotstyle == "graphics") {
    ######## plot if plotstyle=='graphics' ########
    
    if(!discrete)
    {
      #main plotting
      reshist <- hist(mydata, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, col = datacol, probability = probability, ...)
      
      #plot fitted densities (line)
      for(i in 1:nft)
          lines(sfin, fitteddens[,i], lty=fitlty[i], col=fitcol[i], ...)
 
      #plot empirical density
      if(demp)
        lines(density(mydata)$x, density(mydata)$y * scalefactor, col=dempcol)
      
      if (addlegend)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
    }else # so if discrete
    {
      #main plotting
      # plotting of an empty histogramm
      reshist <- hist(mydata, main = main, xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, border = "white",
                      probability = probability, ...)
       
      eps <- diff(range(sfin))/200
      if(fittype %in% c("l", "o"))
      {
        #plot fitted mass probability functions (line)
        for(i in 1:nft)
          lines(sfin+(i)*eps, fitteddens[,i], lty=fitlty[i], col=fitcol[i], type="h", ...)
        #plot empirical mass probabilty function
        if(demp)
        {
          empval <- sort(unique(mydata))
          empprob <- as.numeric(table(mydata))/length(mydata) * scalefactor
          lines(empval, empprob, col=dempcol, type="h")
        }
      }
      if(fittype %in% c("p", "o"))  
      {
        #plot fitted mass probability functions (point)
        for(i in 1:nft)
          points(sfin+(i)*eps, fitteddens[,i], col=fitcol[i], pch=1)
        #plot empirical density
        if(demp)
        {
          empval <- sort(unique(mydata))
          empprob <- as.numeric(table(mydata))/length(mydata) * scalefactor
          points(empval, empprob, col=dempcol, pch=1)
        }  
      }
      
      if (addlegend && fittype %in% c("l", "o"))
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, lty=fitlty, col=fitcol, ...)
      if (addlegend && fittype == "p")
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, pch=1, col=fitcol, ...)
    }
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
    
    if (!discrete)
    {
      # structure the fitteddens in a relevant data.frame
      fitteddens <- as.data.frame(fitteddens)
      colnames(fitteddens) <- unlist(lapply(ft, function(X) X["distname"]))
      fitteddens <- stack(fitteddens)
      fitteddens$sfin <- sfin   # sfin is recycled in the standard fashion
      fitteddens$ind <- factor(fitteddens$ind, levels = unique(fitteddens$ind))   # reorder levels in the appearance order of the input
      if(demp) # bind empirical data if demp is TRUE
        fitteddens <- rbind(fitteddens, data.frame(values = density(mydata)$y * scalefactor, ind = "demp", sfin = density(mydata)$x))
      
      histdata <- data.frame(values = mydata, ind = "hist", sfin = mydata) # the added data must have the same column names as the main data to be compatible with ggplot
      
        ggdenscomp <-
        ggplot2::ggplot(fitteddens, ggplot2::aes_(quote(sfin), quote(values), group = quote(ind), colour = quote(ind))) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::ggtitle(main) +
        ggplot2::coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(ylim[1], ylim[2])) +
          {if(probability) ggplot2::geom_histogram(data = histdata, ggplot2::aes_(quote(values), quote(..density..)), breaks = reshist$breaks, boundary = 0, show.legend = FALSE, col = "black", alpha = 1, fill = datacol)
            else ggplot2::geom_histogram(data = histdata, ggplot2::aes_(quote(values), quote(..count..)), breaks = reshist$breaks, boundary = 0, show.legend = FALSE, col = "black", alpha = 1, fill = datacol)} +
        ggplot2::geom_line(data = fitteddens, ggplot2::aes_(linetype = quote(ind), colour = quote(ind)), size = 0.4) +
        ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
        ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend)) else ggplot2::theme(legend.position = "none")} +
        ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
        ggplot2::scale_linetype_manual(values = fitlty, labels = legendtext)
      return(ggdenscomp)
      
    } else
    {
      eps <- diff(range(sfin))/200
      
      # structure the fitteddens in a relevant data.frame
      fitteddens <- as.data.frame(fitteddens)
      colnames(fitteddens) <- unlist(lapply(ft, function(X) X["distname"]))
      fitteddens <- stack(fitteddens)
      fitteddens$ind <- factor(fitteddens$ind, levels = unique(fitteddens$ind))   # reorder levels in the appearance order of the input
      fitteddens$sfin <- sfin + sapply(fitteddens$ind, function(X) which(X == levels(fitteddens$ind))) *eps   # sfin is recycled in the standard fashion
      
      if(demp) # bind empirical data if demp is TRUE
        fitteddens <- rbind(fitteddens, data.frame(values = as.numeric(table(mydata))/length(mydata) * scalefactor, 
                                                   ind = "demp", 
                                                   sfin = as.numeric(names(table(mydata)))))
      
      ggdenscomp <-
        ggplot2::ggplot(fitteddens, ggplot2::aes_(quote(sfin), quote(values), group = quote(ind), colour = quote(ind))) +
        ggplot2::xlab(xlab) +
        ggplot2::ylab(ylab) +
        ggplot2::ggtitle(main) +
        ggplot2::coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(ylim[1], ylim[2])) +
        {if(fittype %in% c("l", "o")) ggplot2::geom_segment(data = fitteddens, ggplot2::aes_(x = quote(sfin), xend = quote(sfin), y = 0, yend = quote(values), linetype = quote(ind)))} +
        {if(fittype %in% c("p", "o")) ggplot2::geom_point(data = fitteddens, ggplot2::aes_(x = quote(sfin), y = quote(values), colour = quote(ind)), shape = 1)} +
        ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
        ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend)) else ggplot2::theme(legend.position = "none")} +
        ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
        ggplot2::scale_linetype_manual(values = fitlty, labels = legendtext)
      return(ggdenscomp)
      
    }
  }
}
