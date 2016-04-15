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

llplot <- function(mlefit, loglik = TRUE, expansion = 1, lseq = 50,
                   back.col = TRUE, nlev = 10, pal.col = terrain.colors(100), ...)
{
  if (!inherits(mlefit, "fitdist"))
    stop("Use only with 'fitdist' objects")
  if(!is.null(mlefit$weights))
    stop("The plot of the fit is not yet available when using weights")
  if(mlefit$method !="mle")
    stop("This plot is only available for distribution fits using maximum likelihood")
  
  data <- mlefit$data
  distr <- mlefit$distname
  np <- length(mlefit$estimate)
  if (np == 1)
  {
    estim.value <- mlefit$estimate
    estim.sd <- mlefit$sd
    plot.arg <- names(mlefit$estimate)
    fix.arg <- mlefit$fix.arg
    llcurve(data, distr, plot.arg = plot.arg, 
            min.arg = estim.value - estim.sd * 2 *expansion, 
            max.arg = estim.value + estim.sd * 2 *expansion, 
            lseq = lseq, fix.arg = fix.arg, loglik = loglik, ...)
  } else # so if np > 1
    if (np == 2)
    {
      estim.value <- mlefit$estimate
      estim.sd <- mlefit$sd
      plot.arg <- names(mlefit$estimate)
      fix.arg <- mlefit$fix.arg
      llsurface(data, distr, plot.arg = plot.arg, 
              min.arg = estim.value - estim.sd * 2 *expansion, 
              max.arg = estim.value + estim.sd * 2 *expansion, 
              lseq = lseq, fix.arg = fix.arg, loglik = loglik,
              col.back = back.col, nlev = nlev, col.pal = pal.col, ...)
      
    } else # so if np > 2
    {
      def.par <- par(no.readonly = TRUE)
      ncombi <- choose(np, 2)
      lay <- lower.tri(matrix(0, (np - 1), (np - 1)), TRUE)
      lay[which(lay, TRUE)] <- 1:ncombi
      layout(lay)
      par(mar = c(5, 4, 0.2, 0.2))
      for (i in 1:(np - 1))
        for (j in (i+1):np)
        {
          plot.arg <- names(mlefit$estimate)[c(i, j)]
          estim.value <- mlefit$estimate[c(i, j)]
          estim.sd <- mlefit$sd[c(i, j)]
          fix.arg <- c(mlefit$fix.arg, as.list(mlefit$estimate[-c(i,j)]))
          llsurface(data, distr, plot.arg = plot.arg, 
                    min.arg = estim.value - estim.sd * 2 *expansion, 
                    max.arg = estim.value + estim.sd * 2 *expansion, 
                    lseq = lseq, fix.arg = fix.arg, loglik = loglik,
                    col.back = back.col, nlev = nlev, col.pal = pal.col, ...)
        }
      par(def.par)
      
    }
}

llsurface <- function(data, distr, plot.arg, min.arg, max.arg,   lseq = 50, fix.arg = NULL,  
                      loglik = TRUE, col.back = TRUE, nlev = 10, col.pal = terrain.colors(100), ...)
{
  stopifnot(is.vector(plot.arg) || length(plot.arg) == 2)
  stopifnot(is.list(fix.arg) || is.null(fix.arg))
  
  #get distrib name
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  
  #sanity check for argument names
  argdistname <- names(formals(ddistname))
  
  m <- match(plot.arg, argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'plot.arg' must specify names which are arguments to 'distr'.")
  m <- match(names(fix.arg), argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'fix.arg' must specify names which are arguments to 'distr'.")
  
  #function to plot
  if (is.vector(data))
  {
    if(loglik)
    {
      f2plot <- function(x, y)
      {
        par <- list(x,y)
        names(par) <- plot.arg
        loglikelihood(par, fix.arg = fix.arg, obs = data, ddistnam = ddistname)
        #sum(log(do.call(ddistname, c(list(data), par, as.list(fix.arg)) ) ) )
      }
    }else
    {
      f2plot <- function(x, y)
      {
        par <- list(x,y)
        names(par) <- plot.arg
        likelihood(par, fix.arg = fix.arg, obs= data, ddistnam=ddistname)
        # prod(do.call(ddistname, c(list(data), as.list(par), as.list(fix.arg)) ) ) 
      }
    }
  } else stop("This function is not yet available for censored data")
  # TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #create x, y and z matrix.
  p1 <- seq(min.arg[1], max.arg[1], length=lseq)
  p2 <- seq(min.arg[2], max.arg[2], length=lseq)
  z <- outer(p1, p2, Vectorize(f2plot, c("x","y")))
  # vectorize is necessary to vectorize the function f2plot
  
  if (col.back)
  {
    image(p1, p2, z, col = col.pal, xlab = plot.arg[1], ylab = plot.arg[2])
    if (nlev > 0)
      contour(p1, p2, z, nlevels = nlev, add = TRUE)
  } else
  {
    contour(p1, p2, z, nlevels = nlev, xlab = plot.arg[1], ylab = plot.arg[2])
  }
  
  invisible()
}


llcurve <- function(data, distr, plot.arg, min.arg, max.arg,   lseq = 50, fix.arg = NULL,  
                    loglik = TRUE, ...)
{
  stopifnot(is.vector(plot.arg) || length(plot.arg) == 1)
  stopifnot(is.list(fix.arg) || is.null(fix.arg))

  
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  
  #sanity check for argument names
  argdistname <- names(formals(ddistname))
  
  m <- match(plot.arg, argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'plot.arg' must specify names which are arguments to 'distr'.")
  m <- match(names(fix.arg), argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'fix.arg' must specify names which are arguments to 'distr'.")
  
  if (is.vector(data))
  {
    #function to plot
    if(loglik)
    {
      f2plot <- function(x)
      {
        par <- list(x)
        names(par) <- plot.arg
        loglikelihood(par, fix.arg=fix.arg, obs = data, ddistnam = ddistname)
      }
    }else
    {
      f2plot <- function(x)
      {
        par <- list(x)
        names(par) <- plot.arg
        likelihood(par, fix.arg=fix.arg, obs = data, ddistnam = ddistname)
      }
    }
  } else stop("This function is not yet available for censored data")
  # TO DO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #create x, y and z matrix.
  p1 <- seq(min.arg[1], max.arg[1], length = lseq)
  y <- sapply(p1, function(x) f2plot(x))
  plot(p1, y, type="l", xlab = plot.arg, ylab = ifelse(loglik, "loglikelihood", "likelihood"), ...)
  invisible()
}


#local definition of loglikelihood
loglikelihood <- function(par, fix.arg, obs, ddistnam) 
  sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )

#local definition of likelihood
likelihood <- function(par, fix.arg, obs, ddistnam) 
  prod(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) 


