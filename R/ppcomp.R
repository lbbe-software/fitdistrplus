#############################################################################
#   Copyright (c) 2012 Christophe Dutang, Aurélie Siberchicot
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
                   ynoise = TRUE, plotstyle = "graphics", ...)
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
    stop("ppcomp is not yet available when using weights")
  
  # check the 'plotstyle' argument
  plotstyle <- match.arg(plotstyle[1], choices = c("graphics", "ggplot"), several.ok = FALSE)
  
  # check data
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
  largedata <- (n > 1e4)
  logxy <- paste(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""), sep="")
  
  # manage default parameters
  nft <- length(ft)
  if (missing(fitcol)) fitcol <- 2:(nft+1)
  if (missing(fitpch)) fitpch <- ifelse(largedata, 1, 21)
  fitcol <- rep(fitcol, length.out=nft)
  fitpch <- rep(fitpch, length.out=nft)
  if (missing(legendtext)) 
    legendtext <- paste("fit",1:nft)
  
  if (missing(xlab))
    xlab <- "Theoretical probabilities"
  if (missing(ylab)) 
    ylab <- "Empirical probabilities"
  if (missing(main)) 
    main <- "P-P plot"
  
  if (use.ppoints)
    obsp <- ppoints(n, a = a.ppoints)
  else
    obsp <- (1:n) / n
  
  # computation of each fitted distribution
  comput.fti <- function(i)
  {
    fti <- ft[[i]]
    para <- c(as.list(fti$estimate), as.list(fti$fix.arg))
    distname <- fti$distname
    pdistname <- paste("p", distname, sep="")
    do.call(pdistname, c(list(q=sdata), as.list(para)))
  }
  fittedprob <- sapply(1:nft, comput.fti)
  if(NCOL(fittedprob) != nft || NROW(fittedprob) != length(sdata))
    stop("problem when computing fitted probabilities.")
  
  # check limits
  if (missing(xlim))
    xlim <- range(fittedprob)
  if (missing(ylim))
    ylim <- range(obsp)
  
  if(plotstyle == "graphics") {
    ######## plot if plotstyle=='graphics' ########
    
    #main plotting
    if(!largedata)
      resquant <- plot(fittedprob[,1], obsp, main=main, xlab=xlab, ylab=ylab, log=logxy,
                       pch = fitpch[1], xlim=xlim, ylim=ylim, col=fitcol[1], type="p", ...)
    else
      resquant <- plot(fittedprob[,1], obsp, main=main, xlab=xlab, ylab=ylab, log=logxy,
                       lty = fitpch[1], xlim=xlim, ylim=ylim, col=fitcol[1], type="l", ...)
    
    #plot other fitted probabilities
    if(nft > 1 && !ynoise && !largedata)
      for(i in 2:nft)
        points(fittedprob[,i], obsp, pch=fitpch[i], col=fitcol[i], ...)
    if(nft > 1 && ynoise && !largedata)
      for(i in 2:nft)
        points(fittedprob[,i], obsp*(1 + rnorm(n, 0, 0.01)), pch=fitpch[i], col=fitcol[i], ...)
    if(largedata)
      for(i in 2:nft)
        lines(fittedprob[,i], obsp, col=fitcol[i], lty = fitpch[i], ...)
    
    if(line01)
      abline(0, 1, lty=line01lty, col=line01col)
    
    if(addlegend)
    {
      if(!largedata)
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, col=fitcol, pch = fitpch, ...)
      else
        legend(x=xlegend, y=ylegend, bty="n", legend=legendtext, col=fitcol, lty = fitpch, ...)
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
    
    # structure the fittedprob in a relevant data.frame
    fittedprob <- as.data.frame(fittedprob)
    colnames(fittedprob) <- unlist(lapply(ft, function(X) X["distname"]))
    fittedprob <- stack(fittedprob)
    fittedprob$obsp <- obsp   # obsp is recycled in the standard fashion
    fittedprob$ind <- factor(fittedprob$ind, levels = unique(fittedprob$ind))   # reorder levels in the appearance order of the input
    if(nft > 1 && ynoise && !largedata) {
      fittedprob$obsp <- fittedprob$obsp*(1 + rnorm(n*nft, 0, 0.01))
    }
    
    ggppcomp <-
      ggplot2::ggplot(data = fittedprob, ggplot2::aes_(quote(values), quote(obsp), group = quote(ind), colour = quote(ind), shape = quote(ind))) +
      ggplot2::xlab(xlab) +
      ggplot2::ylab(ylab) +
      ggplot2::ggtitle(main) +
      ggplot2::coord_cartesian(xlim = c(xlim[1], xlim[2]), ylim = c(ylim[1], ylim[2])) +
      {if(!largedata) ggplot2::geom_point() else ggplot2::geom_line(ggplot2::aes_(linetype = quote(ind)))} +
      
      {if(addlegend) ggplot2::theme(legend.position = c(xlegend, ylegend)) else ggplot2::theme(legend.position = "none")} +
      ggplot2::scale_color_manual(values = fitcol, labels = legendtext) +
      ggplot2::scale_shape_manual(values = fitpch, labels = legendtext) +
      ggplot2::scale_linetype_manual(values = fitpch, labels = legendtext) +
      
      ggplot2::guides(colour = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(shape = ggplot2::guide_legend(title = NULL)) +
      ggplot2::guides(linetype = ggplot2::guide_legend(title = NULL)) +
      
      {if(line01) ggplot2::geom_abline(intercept = 0, slope = 1)} +
      {if(xlogscale) ggplot2::scale_x_continuous(trans='log10')} +
      {if(ylogscale) ggplot2::scale_y_continuous(trans='log10')}
    
    return(ggppcomp)
  }
}
