#############################################################################
#   Copyright (c) 2015 Christophe Dutang and Marie Laure Delignette-Muller
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
### (log)-likelihood surface/line
###
###         R functions
### 

llsurface <- function(plot.min, plot.max, plot.arg, plot.np=50, fix.arg=NULL, obs, distr, 
    loglik=TRUE, plot.type=c("contour", "filled.contour", "persp", "image"), 
    enhance=TRUE, ...)
{
  stopifnot(is.vector(plot.arg) || length(plot.arg) == 2)
  stopifnot(is.list(fix.arg) || is.null(fix.arg))
  plot.type <- match.arg(plot.type, c("contour", "filled.contour", "persp", "image"))
  
  
  #get distrib name
  ddistname <- paste("d",distr,sep="") 
  
  #sanity check for argument names
  argdistname <- names(formals(ddistname))
  
  m <- match(plot.arg, argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'plot.arg' must specify names which are arguments to 'distr'.")
  m <- match(names(fix.arg), argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'fix.arg' must specify names which are arguments to 'distr'.")
  
  #function to plot
  if(loglik)
  {
    f2plot <- function(x, y)
    {
      par <- list(x,y)
      names(par) <- plot.arg
      loglikelihood(par, fix.arg=fix.arg, obs=obs, ddistnam=ddistname)
    }
  }else
  {
    f2plot <- function(x, y)
    {
      par <- list(x,y)
      names(par) <- plot.arg
      likelihood(par, fix.arg=fix.arg, obs=obs, ddistnam=ddistname)
    }
  }
  
  #create x, y and z matrix.
  p1 <- seq(plot.min[1], plot.max[1], length=plot.np)
  p2 <- seq(plot.min[2], plot.max[2], length=plot.np)
  z <- t(sapply(p1, function(x) sapply(p2, function(y) f2plot(x, y))))
  
  if(enhance)
  {
    nrz <- nrow(z)
    ncz <- ncol(z)
    # Create a function interpolating colors in the range of specified colors
    jet.colors <- colorRampPalette( c("blue", "green", "yellow", "orange", "red") )
    # Generate the desired number of colors from this palette
    nbcol <- 100
    color <- jet.colors(nbcol)
    # Compute the z-value at the facet centres
    z[is.infinite(z)] <- sort(z, decreasing=TRUE)[2]
    zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
    # Recode facet z-values into color indices
    facetcol <- cut(zfacet, nbcol)
    col <- color[facetcol]
    
    switch(plot.type,
           contour = contour(p1, p2, z, col=color, ...),
           filled.contour = filled.contour(p1, p2, z, col=col, ...),
           image = image(p1, p2, z, col=color, ...),
           persp = persp(p1, p2, z, zlab="Log-likelihood", col=col, ...))
  }else
  {
    switch(plot.type,
           contour = contour(p1, p2, z, ...),
           filled.contour = filled.contour(p1, p2, z, ...),
           image = image(p1, p2, z, ...),
           persp = persp(p1, p2, z, zlab="Log-likelihood",...))
  }
  invisible()
}


llcurve <- function(plot.min, plot.max, plot.arg, plot.np=50, fix.arg=NULL, obs, distr, 
                      loglik=TRUE, plot.type=c("line", "points"), enhance=TRUE, ...)
{
  stopifnot(is.vector(plot.arg) || length(plot.arg) == 1)
  stopifnot(is.list(fix.arg) || is.null(fix.arg))
  plot.type <- match.arg(plot.type, c("line", "points"))
  
  
  #get distrib name
  ddistname <- paste("d",distr,sep="") 
  
  #sanity check for argument names
  argdistname <- names(formals(ddistname))
  
  m <- match(plot.arg, argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'plot.arg' must specify names which are arguments to 'distr'.")
  m <- match(names(fix.arg), argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'fix.arg' must specify names which are arguments to 'distr'.")
  
  #function to plot
  if(loglik)
  {
    f2plot <- function(x)
    {
      par <- list(x)
      names(par) <- plot.arg
      loglikelihood(par, fix.arg=fix.arg, obs=obs, ddistnam=ddistname)
    }
  }else
  {
    f2plot <- function(x)
    {
      par <- list(x)
      names(par) <- plot.arg
      likelihood(par, fix.arg=fix.arg, obs=obs, ddistnam=ddistname)
    }
  }
  
  #create x, y and z matrix.
  p1 <- seq(plot.min[1], plot.max[1], length=plot.np)
  y <- sapply(p1, function(x) f2plot(x))
  
  if(enhance)
  {
    switch(plot.type,
           line = plot(p1, y, col="red", type="l", ...),
           points = plot(p1, y, col="red", type="p", ...))
  }else
  {
    switch(plot.type,
           line = plot(p1, y, type="l", ...),
           points = plot(p1, y, type="p", ...))
  }
  invisible()
}


#local definition of loglikelihood
loglikelihood <- function(par, fix.arg, obs, ddistnam) 
  sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )

#local definition of likelihood
likelihood <- function(par, fix.arg, obs, ddistnam) 
  prod(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) 


