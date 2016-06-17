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


cdfband <- function(b, xlim, ylim, xlogscale = FALSE, ylogscale = FALSE, main, xlab, ylab, 
                    datapch, datacol, fitlty, fitcol, addlegend = FALSE, legendtext, xlegend = "bottomright", 
                    ylegend = NULL, horizontals = TRUE, verticals = FALSE, do.points = TRUE, 
                    use.ppoints = TRUE, a.ppoints = 0.5, lines01 = FALSE, CI.type = "two.sided", 
                    CI.level = 0.95, CI.col="red", CI.lty=2, CI.fill=FALSE, ...)
{
  if(!inherits(b, "bootdist"))
  {
    stop("argument b must be a 'bootdist' objects")
  }
  
  CI.type <- match.arg(CI.type, c("two.sided", "less", "greater"))
  CI.level <- CI.level[1]
  
  if (missing(datapch)) datapch <- 16
  if (missing(datacol)) datacol <- "black"
  if (missing(fitcol)) fitcol <- 2
  if (missing(fitlty)) fitlty <- 1
  
  cdfval <- function(x)
  {  
    calcp <- function(i)
    {
      parai <- c(as.list(b$estim[i, ]), as.list(b$fitpart$fix.arg))
      do.call(pdistname, c(list(q=x), as.list(parai)))
    }
    res <- t(sapply(1:b$nbboot, calcp))
    rownames(res) <- 1:b$nbboot
    colnames(res) <- paste0("x=", x)
    res
  }
  lowx <- ifelse(min(b$fitpart$data) < 0, min(b$fitpart$data)*1.5, min(b$fitpart$data)*.5)
  uppx <- ifelse(max(b$fitpart$data) < 0, max(b$fitpart$data)*.5, max(b$fitpart$data)*1.5)
  x <- seq(lowx, uppx, length=101)
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
  
  cdfcomp(b1$fitpart, xlim=xlim, ylim=ylim, xlogscale = xlogscale, ylogscale = ylogscale, 
          main=main, xlab=xlab, ylab=ylab, datapch=datapch, datacol=datacol, fitlty=fitlty, 
          fitcol=fitcol, addlegend = addlegend, legendtext=legendtext, xlegend = xlegend, 
          ylegend = ylegend, horizontals = horizontals, verticals = verticals, do.points = do.points, 
          use.ppoints = use.ppoints, a.ppoints = a.ppoints, lines01 = lines01)
  
  if(!CI.fill)
  {
    matlines(x, CIband, col=CI.col, lty=CI.lty, ...) 
  }
  else
  {
    if(CI.type == "two.sided")
    {
      polygon(c(x, rev(x)), c(CIband[,2], rev(CIband[,1])), col=CI.col, border=CI.col)
    }
    if(CI.type == "less")
    {
      polygon(c(x, uppx, uppx), c(CIband, 1, 0), col=CI.col, border=CI.col)
    }
    if(CI.type == "greater")
    {
      polygon(c(x, lowx, lowx), c(CIband, 1, 0), col=CI.col, border=CI.col)
    }
    
    plot(ecdf(b$fitpart$data), col=datacol, add=TRUE, do.points=do.points, verticals = verticals, pch=datapch)
  }
}