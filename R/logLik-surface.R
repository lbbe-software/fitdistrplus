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
                   back.col = TRUE, nlev = 10, pal.col = terrain.colors(100),
                   fit.show = FALSE, fit.pch = 4, ...)
{
  if (!inherits(mlefit, c("fitdist", "fitdistcens")))
    stop("Use only with 'fitdist' or 'fitdistcens' objects")
  if(inherits(mlefit, "fitdist"))
  {
    if(mlefit$method !="mle")
      stop("This plot is only available for distribution fits using maximum likelihood")
    data <- mlefit$data
    
  } else # censored data
  {
    data <- mlefit$censdata
  }
  
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
            lseq = lseq, fix.arg = fix.arg, loglik = loglik, weights = mlefit$weights, ...)
    if (fit.show) points(estim.value, ifelse(loglik, mlefit$loglik, exp(mlefit$loglik)), 
                         pch = fit.pch)
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
              back.col = back.col, nlev = nlev, pal.col = pal.col,
              weights = mlefit$weights, ...)
    if (fit.show) points(estim.value[1], estim.value[2], pch = fit.pch)
      
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
                    back.col = back.col, nlev = nlev, pal.col = pal.col, 
                    weights = mlefit$weights, ...)
          if (fit.show) points(estim.value[1], estim.value[2], pch = fit.pch)
        }
      par(def.par)
    }
  invisible()
}

llsurface <- function(data, distr, plot.arg, min.arg, max.arg,   lseq = 50, fix.arg = NULL,  
            loglik = TRUE, back.col = TRUE, nlev = 10, pal.col = terrain.colors(100), weights = NULL, ...)
{
  stopifnot(is.vector(plot.arg) || length(plot.arg) == 2)
  stopifnot(is.list(fix.arg) || is.null(fix.arg))

  if(!is.null(weights))
  {
    if(any(weights < 0))
      stop("weights should be a vector of numerics greater than 0")
    if(length(weights) != NROW(data))
      stop("weights should be a vector with a length equal to the observation number")
  } else
  {
    weights <- rep(1, NROW(data))
  }
  
    
  if (is.vector(data))
  {
    cens <- FALSE
  } else
  {
    cens <- TRUE
    # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector, 
    #   icens (interval censored) = dataframe with left and right 
    # and ncens (not censored) = vector
    censdata <- data
    irow.lcens <- is.na(censdata$left) # rows corresponding to lcens data
    lcens <- censdata[irow.lcens, ]$right
    if (any(is.na(lcens)) )
      stop("An observation cannot be both right and left censored, coded with two NA values")
    irow.rcens <- is.na(censdata$right)  # rows corresponding to rcens data
    rcens <- censdata[irow.rcens, ]$left
    irow.ncens <- censdata$left==censdata$right & !is.na(censdata$left) & 
      !is.na(censdata$right)  # rows corresponding to ncens data
    ncens<-censdata[irow.ncens, ]$left
    irow.icens <- censdata$left!=censdata$right & !is.na(censdata$left) & 
      !is.na(censdata$right)  # rows corresponding to icens data
    icens<-censdata[irow.icens, ]
  }
  
  #get distrib name
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  if (!exists(ddistname, mode="function"))
  stop(paste("The ", ddistname, " function must be defined"))
  if (cens)
  {
    pdistname <- paste("p", distname, sep="")
    if (!exists(pdistname, mode="function"))
      stop(paste("The ", pdistname, " function must be defined"))
    
  }
  
  
  #sanity check for argument names
  argdistname <- names(formals(ddistname))
  
  m <- match(plot.arg, argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'plot.arg' must specify names which are arguments to 'distr'.")
  m <- match(names(fix.arg), argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'fix.arg' must specify names which are arguments to 'distr'.")
  
  #function to plot
  if (!cens)
  {
    if(loglik)
    {
      f2plot <- function(x, y)
      {
        par <- list(x,y)
        names(par) <- plot.arg
        loglikelihood(par, fix.arg = fix.arg, obs = data, ddistnam = ddistname, weights = weights)
        #sum(log(do.call(ddistname, c(list(data), par, as.list(fix.arg)) ) ) )
      }
    }else
    {
      f2plot <- function(x, y)
      {
        par <- list(x,y)
        names(par) <- plot.arg
        likelihood(par, fix.arg = fix.arg, obs= data, ddistnam = ddistname, weights = weights)
        # prod(do.call(ddistname, c(list(data), as.list(par), as.list(fix.arg)) ) ) 
      }
    }
  } else # for censored data
  {
    if(loglik)
    {
      f2plot <- function(x, y)
      {
        par <- list(x,y)
        names(par) <- plot.arg
        loglikelihoodcens(par, fix.arg = fix.arg, rcens = rcens, lcens = lcens, icens = icens,
                          ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, weights = weights,
                          irow.ncens = irow.ncens, irow.lcens = irow.lcens, 
                          irow.rcens = irow.rcens, irow.icens = irow.icens)
      }
    }else
    {
      f2plot <- function(x, y)
      {
        par <- list(x,y)
        names(par) <- plot.arg
        likelihoodcens(par, fix.arg = fix.arg, rcens = rcens, lcens = lcens, icens = icens,
                       ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, weights = weights,
                       irow.ncens = irow.ncens, irow.lcens = irow.lcens, 
                       irow.rcens = irow.rcens, irow.icens = irow.icens)
      }
    }
  }  
  #create x, y and z matrix.
  p1 <- seq(min.arg[1], max.arg[1], length=lseq)
  p2 <- seq(min.arg[2], max.arg[2], length=lseq)
  z <- outer(p1, p2, Vectorize(f2plot, c("x","y")))
  # vectorize is necessary to vectorize the function f2plot
  
  if (back.col)
  {
    image(p1, p2, z, col = pal.col, xlab = plot.arg[1], ylab = plot.arg[2])
    if (nlev > 0)
      contour(p1, p2, z, nlevels = nlev, add = TRUE)
  } else
  {
    contour(p1, p2, z, nlevels = nlev, xlab = plot.arg[1], ylab = plot.arg[2])
  }
  
  invisible()
}


llcurve <- function(data, distr, plot.arg, min.arg, max.arg,   lseq = 50, fix.arg = NULL,  
                    loglik = TRUE, weights = NULL, ...)
{
  stopifnot(is.vector(plot.arg) || length(plot.arg) == 1)
  stopifnot(is.list(fix.arg) || is.null(fix.arg))

  if(!is.null(weights))
  {
    if(any(weights < 0))
      stop("weights should be a vector of numerics greater than 0")
    if(length(weights) != NROW(data))
      stop("weights should be a vector with a length equal to the observation number")
  } else
  {
    weights <- rep(1, NROW(data))
  }
  
  
  if (is.vector(data))
  {
    cens <- FALSE
  } else
  {
    cens <- TRUE
    # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector, 
    #   icens (interval censored) = dataframe with left and right 
    # and ncens (not censored) = vector
    censdata <- data
    irow.lcens <- is.na(censdata$left) # rows corresponding to lcens data
    lcens <- censdata[irow.lcens, ]$right
    if (any(is.na(lcens)) )
      stop("An observation cannot be both right and left censored, coded with two NA values")
    irow.rcens <- is.na(censdata$right)  # rows corresponding to rcens data
    rcens <- censdata[irow.rcens, ]$left
    irow.ncens <- censdata$left==censdata$right & !is.na(censdata$left) & 
      !is.na(censdata$right)  # rows corresponding to ncens data
    ncens<-censdata[irow.ncens, ]$left
    irow.icens <- censdata$left!=censdata$right & !is.na(censdata$left) & 
      !is.na(censdata$right)  # rows corresponding to icens data
    icens<-censdata[irow.icens, ]
  }
  
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
 
  if (cens)
  {
    pdistname <- paste("p", distname, sep="")
    if (!exists(pdistname, mode="function"))
      stop(paste("The ", pdistname, " function must be defined"))
    
  }
  
  #sanity check for argument names
  argdistname <- names(formals(ddistname))
  
  m <- match(plot.arg, argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'plot.arg' must specify names which are arguments to 'distr'.")
  m <- match(names(fix.arg), argdistname)
  if (any(is.na(m))) #check unexpected names
    stop("'fix.arg' must specify names which are arguments to 'distr'.")
  
  if (!cens)
  {
    #function to plot
    if(loglik)
    {
      f2plot <- function(x)
      {
        par <- list(x)
        names(par) <- plot.arg
        loglikelihood(par, fix.arg=fix.arg, obs = data, ddistnam = ddistname, weights = weights)
      }
    }else
    {
      f2plot <- function(x)
      {
        par <- list(x)
        names(par) <- plot.arg
        likelihood(par, fix.arg=fix.arg, obs = data, ddistnam = ddistname, weights = weights)
      }
    }
  } else # for censored data
  {
    if(loglik)
    {
      f2plot <- function(x)
      {
        par <- list(x)
        names(par) <- plot.arg
        loglikelihoodcens(par, fix.arg = fix.arg, rcens = rcens, lcens = lcens, icens = icens,
                          ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, weights = weights,
                          irow.ncens = irow.ncens, irow.lcens = irow.lcens, 
                          irow.rcens = irow.rcens, irow.icens = irow.icens)
      }
    }else
    {
      f2plot <- function(x)
      {
        par <- list(x)
        names(par) <- plot.arg
        likelihoodcens(par, fix.arg = fix.arg, rcens = rcens, lcens = lcens, icens = icens,
                          ncens = ncens, ddistnam = ddistname, pdistnam = pdistname, weights = weights,
                          irow.ncens = irow.ncens, irow.lcens = irow.lcens, 
                          irow.rcens = irow.rcens, irow.icens = irow.icens)
      }
    }
  }
  
  #create x, y matrix.
  p1 <- seq(min.arg[1], max.arg[1], length = lseq)
  y <- sapply(p1, function(x) f2plot(x))
  plot(p1, y, type="l", xlab = plot.arg, ylab = ifelse(loglik, "loglikelihood", "likelihood"), ...)
  invisible()
}


#local definition of loglikelihood
loglikelihood <- function(par, fix.arg, obs, ddistnam, weights) 
  sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )

#local definition of likelihood
likelihood <- function(par, fix.arg, obs, ddistnam, weights) 
  prod(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) )^weights ) 

#local definition of loglikelihood for censored data
loglikelihoodcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam, 
                              weights, irow.ncens, irow.lcens, irow.rcens, irow.icens)
{
  p1 <- log(do.call(ddistnam, c(list(ncens), as.list(par), as.list(fix.arg))))
  p2 <- log(do.call(pdistnam, c(list(lcens), as.list(par), as.list(fix.arg)))) 
  p3 <- log(1-do.call(pdistnam, c(list(rcens), as.list(par), as.list(fix.arg))))
  p4 <- log(do.call(pdistnam, c(list(icens$right), as.list(par), as.list(fix.arg))) - 
              do.call(pdistnam, c(list(icens$left), as.list(par), as.list(fix.arg))) )
  sum(weights[irow.ncens] * p1) + 
    sum(weights[irow.lcens] * p2) + 
    sum(weights[irow.rcens] * p3) + 
    sum(weights[irow.icens] * p4) 
}

#local definition of likelihood for censored data
likelihoodcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam, 
                           weights, irow.ncens, irow.lcens, irow.rcens, irow.icens)
{
  p1 <- do.call(ddistnam, c(list(ncens), as.list(par), as.list(fix.arg)))
  p2 <- do.call(pdistnam, c(list(lcens), as.list(par), as.list(fix.arg)))
  p3 <- 1-do.call(pdistnam, c(list(rcens), as.list(par), as.list(fix.arg)))
  p4 <- do.call(pdistnam, c(list(icens$right), as.list(par), as.list(fix.arg))) - 
              do.call(pdistnam, c(list(icens$left), as.list(par), as.list(fix.arg))) 
  prod(p1^weights[irow.ncens]) * 
    prod(p2^weights[irow.lcens]) * 
    prod(p3^weights[irow.rcens]) * 
    prod(p4^weights[irow.icens]) 
}

