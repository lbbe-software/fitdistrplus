#############################################################################
#   Copyright (c) 2016 Marie Laure Delignette-Muller, Christophe Dutang
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
### plot cumulative distribution functions with confidence interval band
###
###         R functions
###


CIcdfplot <- function(b, CI.output, CI.type = "two.sided", CI.level = 0.95, CI.col = "red", CI.lty = 2, 
                    CI.fill = NULL, CI.only = FALSE, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, 
                    xlab, ylab, datapch, datacol, fitlty, fitcol, horizontals = TRUE, verticals = FALSE, 
                    do.points = TRUE, use.ppoints = TRUE, a.ppoints = 0.5, lines01 = FALSE, ...)
{
  if(inherits(b, "bootdist"))
  {
    cens <- FALSE
  } else
  if(inherits(b, "bootdistcens"))
  {
    cens <- TRUE
  } else
  {
    stop("argument b must be a 'bootdist' or a `bootdistcens` object")
  }
  if(missing(CI.output))
    stop("argument CI.output must be specified: either 'probability' or 'quantile'.")
  CI.output <- match.arg(CI.output, c("probability", "quantile"))
  CI.type <- match.arg(CI.type, c("two.sided", "less", "greater"))
  CI.level <- CI.level[1]
  #compute lower and upper value for the area
  if (!cens)
  {
    mydat <- b$fitpart$data
    n <- length(mydat)
    xmin <- min(mydat)
    xmax <- max(mydat)
  } else
  {
    censdata <- b$fitpart$censdata
    n <- nrow(censdata)
    xmin <- min(c(censdata$left, censdata$right), na.rm=TRUE)
    xmax <- max(c(censdata$left, censdata$right), na.rm=TRUE)
  }
  if (missing(xlim)) xlim <- c(xmin, xmax)
  lowx <- min(xlim[1], ifelse(xmin < 0, xmin*1.5, xmin*.5))
  uppx <- max(xlim[2], ifelse(xmax < 0, xmax*.5, xmax*1.5))

  if(missing(ylim)) ylim <- c(0, 1)
  
  if(!is.logical(CI.only))
    stop("argument CI.only must be a logical")
  
  #default values (same as cdfcomp())
  if (missing(datapch)) datapch <- 16
  if (missing(datacol)) datacol <- "black"
  if (missing(fitcol)) fitcol <- 2
  if (missing(fitlty)) fitlty <- 1
   if (missing(xlab)) 
  {
     if (!cens)
       xlab <- ifelse(xlogscale, "data in log scale", "data")
     else
       xlab <- ifelse(xlogscale, "censored data in log scale", "censored data")
  }
  if (missing(ylab)) ylab <- "CDF"
  if (missing(main)) main <- ifelse(CI.only, "Theoretical CDF with CI", "Empirical and theoretical CDF with CI")

  #get name and cdf name
  distname <- b$fitpart$distname
  pdistname <- paste("p",distname,sep="")
  qdistname <- paste("q",distname,sep="")
  if (!exists(pdistname, mode="function") && CI.output == "probability")
    stop(paste("The ", pdistname, " function must be defined"))
  if (!exists(qdistname, mode="function") && CI.output == "quantile")
    stop(paste("The ", qdistname, " function must be defined"))
  
  #compute c.d.f. values on bootstraped parameters
  if(CI.output == "probability")
  {
    cdfval <- function(x)
    {  
      calcp <- function(i)
      {
        parai <- c(as.list(b$estim[i, ]), as.list(b$fitpart$fix.arg))
        do.call(pdistname, c(list(x), as.list(parai)))
      }
      
      res <- t(sapply(1:b$nbboot, calcp))
      rownames(res) <- 1:b$nbboot
      colnames(res) <- paste0("x=", x)
      res
    }
    x <- seq(lowx, uppx, length=501)
    
    #compute quantiles on c.d.f. 
    if (CI.type == "two.sided")
    {
      alpha <- (1-CI.level)/2
      CIband <- t(apply(cdfval(x), MARGIN=2, quantile, probs=c(alpha, 1-alpha), na.rm=TRUE))
      colnames(CIband) <- format.perc(c(alpha, 1-alpha), 3)
    }else if (CI.type == "less")
    {
      CIband <- as.matrix(apply(cdfval(x), MARGIN=2, quantile, probs=CI.level, na.rm=TRUE))
      colnames(CIband) <- format.perc(CI.level, 3)
    }else
    {
      CIband <- as.matrix(apply(cdfval(x), MARGIN=2, quantile, probs=1-CI.level, na.rm=TRUE))
      colnames(CIband) <- format.perc(1-CI.level, 3)
    }
  }else #CI.output == "quantile"
  {
    qval <- function(p)
    {  
      calcp <- function(i)
      {
        parai <- c(as.list(b$estim[i, ]), as.list(b$fitpart$fix.arg))
        do.call(qdistname, c(list(p), as.list(parai)))
      }
      
      res <- t(sapply(1:b$nbboot, calcp))
      rownames(res) <- 1:b$nbboot
      colnames(res) <- paste0("p=", p)
      res
    }
    #compute lower and upper value for the area
#    p <- seq(sqrt(.Machine$double.eps), 1- sqrt(.Machine$double.eps), length=101)
    p <- seq(0.0001, 1- 0.0001, length=501)
    
    #compute quantiles on c.d.f. 
    if (CI.type == "two.sided")
    {
      alpha <- (1-CI.level)/2
      CIband <- t(apply(qval(p), MARGIN=2, quantile, probs=c(alpha, 1-alpha), na.rm=TRUE))
      colnames(CIband) <- format.perc(c(alpha, 1-alpha), 3)
    }else if (CI.type == "less")
    {
      CIband <- as.matrix(apply(qval(p), MARGIN=2, quantile, probs=1-CI.level, na.rm=TRUE))
      colnames(CIband) <- format.perc(CI.level, 3)
    }else
    {
      CIband <- as.matrix(apply(qval(p), MARGIN=2, quantile, probs=CI.level, na.rm=TRUE))
      colnames(CIband) <- format.perc(1-CI.level, 3)
    }
  }
  
  #temp var to open a graphic (if needed)
  logxy <- paste0(ifelse(xlogscale,"x",""), ifelse(ylogscale,"y",""))
  
  ##### plot ####
  #open graphic window
  plot(0, 0, main=main, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim,
         log=logxy, type="n")
  
  if (!is.null(CI.fill)) # first fill the band
  {
    if(CI.output == "probability")
    {
      if(CI.type == "two.sided")
        polygon(c(x, rev(x)), c(CIband[,2], rev(CIband[,1])), col=CI.fill, border=CI.fill, ...)
      else if(CI.type == "less")
        polygon(c(x, uppx, uppx), c(CIband, 1, 0), col=CI.fill, border=CI.fill, ...)
      else #if(CI.type == "greater")
        polygon(c(x, lowx, lowx), c(CIband, 1, 0), col=CI.fill, border=CI.fill, ...)
      
    }else #CI.output == "quantile"
    {
      if(CI.type == "two.sided")
        polygon(c(CIband[,2], rev(CIband[,1])), c(p, rev(p)), col=CI.fill, border=CI.fill, ...)
      else if(CI.type == "less")
        polygon(c(CIband, uppx, uppx), c(p, 1, 0), col=CI.fill, border=CI.fill, ...)
      else #if(CI.type == "greater")
        polygon(c(CIband, lowx, lowx), c(p, 1, 0), col=CI.fill, border=CI.fill, ...)
    }
  }

  # add lines for the bounds of the CI
  if(CI.output == "probability")
  {
    matlines(x, CIband, col=CI.col, lty=CI.lty, ...) 
  }else #CI.output == "quantile"
  {
    matlines(CIband, p, col=CI.col, lty=CI.lty, ...) 
  }
        
  if(!CI.only) # add the empirical and fitted distributions
  {
    if (!cens)
    {
      cdfcomp(b$fitpart, xlim=xlim, ylim=ylim, xlogscale = xlogscale, ylogscale = ylogscale, 
              main=main, xlab=xlab, ylab=ylab, datapch=datapch, datacol=datacol, fitlty=fitlty, 
              fitcol=fitcol, horizontals = horizontals, verticals = verticals, do.points = do.points, 
              use.ppoints = use.ppoints, a.ppoints = a.ppoints, lines01 = lines01, addlegend = FALSE,
              add=TRUE)
      
    } else
    {
      cdfcompcens(b$fitpart, xlim=xlim, ylim=ylim, xlogscale = xlogscale, ylogscale = ylogscale, 
              main=main, xlab=xlab, ylab=ylab, datacol=datacol, fitlty=fitlty,  
              fitcol=fitcol, lines01 = lines01, Turnbull.confint = FALSE, addlegend = FALSE, add=TRUE)
      
    }
  }
}