#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### bootstrap in fitdistrplus
###
###         R functions
### 


bootdist <- function (f, bootmethod="param", niter=1001, silent=TRUE, 
                      parallel = c("no", "snow", "multicore"), ncpus)
{ 
  if (niter<10) 
    stop("niter must be an integer above 10")
  bootmethod <- match.arg(bootmethod, c("param", "nonparam"))
  
  parallel <- match.arg(parallel, c("no", "snow", "multicore"))
  if (parallel == "multicore" & .Platform$OS.type == "windows")
  {
    parallel <- "snow"
    warning("As the multicore option is not supported on Windows it was replaced by snow")
  }
  if ((parallel == "snow" | parallel == "multicore") & missing(ncpus)) 
    stop("You have to specify the number of available processors to parallelize 
             the bootstrap")
  
  if (!inherits(f, "fitdist"))
    stop("Use only with 'fitdist' objects")
  
  if(!is.null(f$weights))
    stop("Bootstrap is not yet available when using weights")
  
  #simulate bootstrap data
  if (bootmethod == "param") 
  { # parametric bootstrap
    rdistname <- paste("r", f$distname, sep="")
    if (!exists(rdistname, mode="function"))
      stop(paste("The ", rdistname, " function must be defined"))
    rdata <- do.call(rdistname, c(list(niter*f$n), as.list(f$estimate), f$fix.arg))
    dim(rdata) <- c(f$n, niter)
  }else 
  { # non parametric bootstrap
    rdata <- sample(f$data, size=niter*f$n, replace=TRUE)
    dim(rdata) <- c(f$n, niter)
  }
  
  #compute bootstrap estimates
  foncestim <- switch(f$method, "mle"=mledist, "qme"=qmedist, "mme"=mmedist, 
                      "mge"=mgedist, "mse"=msedist)
  start <- as.list(f$estimate) #a named vector is no longer is accepted as starting values.
  if(is.function(f$fix.arg.fun))
    fix.arg <- f$fix.arg.fun
  else 
    fix.arg <- f$fix.arg
  if (is.null(f$dots) && !is.function(fix.arg))
  {    
    func <- function(iter) 
    {
      res <- try(do.call(foncestim, list(data=rdata[, iter], distr=f$distname, start=start, 
                                         fix.arg=fix.arg, checkstartfix=TRUE)), silent=silent)
      if(inherits(res, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(res$estimate, res$convergence))
    }
  }else if (is.null(f$dots) && is.function(fix.arg))
  {    
    func <- function(iter) 
    {
      fix.arg.iter <- fix.arg(rdata[, iter])
      res <- try(do.call(foncestim, list(data=rdata[, iter], distr=f$distname, start=start, 
                                         fix.arg=fix.arg.iter, checkstartfix=TRUE)), silent=silent)
      if(inherits(res, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(res$estimate, res$convergence))
    }
  }else if(!is.null(f$dots) && !is.function(fix.arg))
  {    
    func <- function(iter) 
    {
      res <- try(do.call(foncestim, c(list(data=rdata[, iter], distr=f$distname, start=start, 
                                           fix.arg=fix.arg, checkstartfix=TRUE), f$dots)), silent=silent)
      if(inherits(res, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(res$estimate, res$convergence))
      
    }
  }else if(!is.null(f$dots) && is.function(fix.arg))
  {    
    func <- function(iter) 
    {
      fix.arg.iter <- fix.arg(rdata[, iter])
      res <- try(do.call(foncestim, c(list(data=rdata[, iter], distr=f$distname, start=start, 
                                           fix.arg=fix.arg.iter, checkstartfix=TRUE), f$dots)), silent=silent)
      if(inherits(res, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(res$estimate, res$convergence))
      
    }
  }else
    stop("wrong implementation in bootdist")
  
  owarn <- getOption("warn")
  oerr <- getOption("show.error.messages")
  options(warn=ifelse(silent, -1, 0), show.error.messages=!silent)
  
  # parallel or sequential computation
  if (parallel != "no") 
  {
    if (parallel == "snow") type <- "PSOCK"
    else if (parallel == "multicore") type <- "FORK"
    clus <- parallel::makeCluster(ncpus, type = type)
    resboot <- parallel::parSapply(clus, 1:niter, func)
    parallel::stopCluster(clus)
  }
  else
  {
    resboot <- sapply(1:niter, func)
  }
  
  options(warn=owarn, show.error.messages=oerr)
  
  rownames(resboot) <- c(names(start), "convergence")
  if (length(resboot[, 1]) > 2)
  {
    estim <- data.frame(t(resboot)[, -length(resboot[, 1])])
    bootCI <- cbind(apply(resboot[-length(resboot[, 1]), ], 1, median, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.025, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.975, na.rm=TRUE))
    colnames(bootCI) <- c("Median", "2.5%", "97.5%")
  }else 
  {
    estim <- as.data.frame(t(resboot)[, -length(resboot[, 1])])
    names(estim) <- names(f$estimate)
    bootCI <- c(median(resboot[-length(resboot[, 1]), ], na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.025, na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.975, na.rm=TRUE)) 
    names(bootCI) <- c("Median", "2.5%", "97.5%") 
  } 
  
  # code of convergence of the optimization function for each iteration
  converg <- t(resboot)[, length(resboot[, 1])]
  
  res <- structure(list(estim=estim, converg=converg, 
                        method=bootmethod, nbboot=niter, CI=bootCI, fitpart=f), 
                   class="bootdist")
  res    
}

print.bootdist <- function(x, ...)
{
  if (!inherits(x, "bootdist"))
    stop("Use only with 'bootdist' objects")
  if (x$method=="param") 
    cat("Parameter values obtained with parametric bootstrap \n")
  else
    cat("Parameter values obtained with nonparametric bootstrap \n")
  print(head(x$estim), ...)    
  nconverg <- length(x$converg[x$converg==0])
  if (nconverg < length(x$converg))
  {
    cat("\n")
    cat("The estimation method converged only for", nconverg, "among", 
        length(x$converg), "iterations \n")
  }
  
}

plot.bootdist <- function(x, main="Bootstrapped values of parameters", enhance=FALSE, 
                          trueval=NULL, rampcol=NULL, nbgrid = 100, nbcol = 100, ...)
{
  if (!inherits(x, "bootdist"))
    stop("Use only with 'bootdist' objects")
  if (dim(x$estim)[2]==1) 
  {
    stripchart(x$estim, method="jitter", main=main, 
               xlab="Bootstrapped values of the parameter", ...)
  }else 
  {
    if(!is.null(trueval)) #true value supplied
      stopifnot(length(trueval) == NCOL(x$estim))
    if(!is.logical(enhance))
      stop("wrong argument enhance for plot.bootdist.")
    if (!enhance)
    {
      if(is.null(trueval)) #no true value supplied
        pairs(data.matrix(x$estim), main=main, ...)
      else #some true value supplied
        pairs4boot(x$estim, main=main, trueval=trueval, enhance=FALSE, ...)
    }else 
    {
      if(is.null(rampcol))
        rampcol <- c("green", "yellow", "orange", "red")
      pairs4boot(x$estim, main=main, trueval=trueval, col4ramp = rampcol, 
                 nbgrid = nbgrid, nbcol = nbcol, ...)
    }
  }
}

summary.bootdist <- function(object, ...)
{
  if (!inherits(object, "bootdist"))
    stop("Use only with 'bootdist' objects")
  
  class(object) <- c("summary.bootdist", class(object))  
  object
}

print.summary.bootdist <- function(x, ...)
{
  
  if (!inherits(x, "summary.bootdist"))
    stop("Use only with 'summary.bootdist' objects")
  
  if (x$method=="param") 
    cat("Parametric bootstrap medians and 95% percentile CI \n")
  else
    cat("Nonparametric bootstrap medians and 95% percentile CI \n")
  print(x$CI)
  
  nconverg <- length(x$converg[x$converg==0])
  if (nconverg < length(x$converg))
  {
    cat("\n")
    cat("The estimation method converged only for", nconverg, "among", 
        length(x$converg), "iterations \n")
  }
}

density.bootdist <- function(..., bw = "nrd0", adjust = 1, kernel = "gaussian")
{
  x <- list(...)
  if(inherits(x, "bootdist"))
  {
    x <- list(x)
  }else if(!is.list(x))
  {
    stop("argument x must be a list of 'bootdist' objects")
  }else
  {
    if(any(sapply(x, function(y) !inherits(y, "bootdist"))))        
      stop("argument x must be a list of 'bootdist' objects")
  }
  nx <- length(x)
  nbpar <- NCOL(x[[1]]$estim)
  denslist <- lapply(1:nx, function(j)
  {
    dres <- lapply(1:nbpar, function(i)
    {
      #compute empirical density, mean, sd
      res <- density(x[[j]]$estim[,i], bw=bw, adjust=adjust, kernel=kernel)
      res$mean <- mean(x[[j]]$estim[,i], na.rm=TRUE)
      res$sd <- sd(x[[j]]$estim[,i], na.rm=TRUE)
      res$distname <- x[[j]]$fitpart$distname
      res
    }
    )
    #name list of densities
    names(dres) <- colnames(x[[j]]$estim)
    dres
  }
  )
  nbboot <- sapply(1:nx, function(j) x[[j]]$nbboot)
  
  structure(denslist, 
            distname=x[[1]]$fitpart$distname,
            nbobject=nx,
            nbboot=nbboot,
            n=x[[1]]$fitpart$n,
            class="density.bootdist")
}

plot.density.bootdist <- function(x, mar=c(4,4,2,1), lty=NULL, col=NULL, lwd=NULL, ...)
{
  if (!inherits(x, "density.bootdist"))
    stop("Use only with 'density.bootdist' objects")
  nbpar <- length(x[[1]])
  nft <-  length(x)
  if(is.null(lty))
    lty <- 1:nft
  else
    lty <- rep(lty, length.out=nft)
  if(is.null(col))
    col <- 1:nft
  else
    col <- rep(col, length.out=nft)
  if(is.null(lwd))
    lwd <- 1.5
  else
    lwd <- rep(lwd, length.out=nft)
  xname <- names(x[[1]])
  m <- ceiling(nbpar/2)
  par(mfrow=c(m, 2), mar=mar)
  for(i in 1:nbpar)
  {
    xlim <- range(sapply(x, function(d) d[[i]]$x))
    ylim <- range(sapply(x, function(d) d[[i]]$y))
    nbboot <- sapply(x, function(d) d[[i]]$n)
    mymean <- signif(sapply(x, function(d) d[[i]]$mean), 3)
    mysd <- signif(sapply(x, function(d) d[[i]]$sd), 3)
    ylim[2] <- ylim[2]*1.2
    
    main <- paste0(x[[1]][[1]]$distname, " distribution - ", xname[i])
    plot(x[[1]][[i]], xlim=xlim, ylim=ylim, lwd=lwd[1], main=main, xlab=xname[i],
         col=col[1], lty=lty[1], ylab="Density of bootstrapped values") 
    if(nft > 1)
    {
      for(j in 2:nft)
      {
        lines(x[[j]][[i]], lty=lty[j], lwd=lwd*(1+(j-1)/10), col=col[j]) 
      }
      myleg <- paste("n=", nbboot, ", mean=", mymean, ", sd=", mysd)
      legend("topleft", lty=lty, col=col, lwd=lwd, legend=myleg, bty="n") 
    }
  }
  invisible(NULL)
}

print.density.bootdist <- function(x, ...)
{
  if (!inherits(x, "density.bootdist"))
    stop("Use only with 'density.bootdist' objects")
  
  nbboot <- paste(attr(x, "nbboot"), collapse=", ")
  cat("\nBootstrap values for: ", attr(x, "distname"), " for ",
      attr(x, "nbobject"), " object(s) with ", 
      nbboot, " bootstrap values (original sample size ",
      attr(x, "n"), ").", sep = "")
  invisible(x)
}

