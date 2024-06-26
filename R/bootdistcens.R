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
### bootstrap in fitdistrplus with censored data
###
###         R functions
### 

bootdistcens <- function (f, niter=1001, silent=TRUE, 
                          parallel = c("no", "snow", "multicore"), ncpus)
{ 
  if (niter<10) 
    stop("niter must be an integer above 10")
  
  parallel <- match.arg(parallel, c("no", "snow", "multicore"))
  if (parallel == "multicore" & .Platform$OS.type == "windows")
  {
    parallel <- "snow"
    warning("As the multicore option is not supported on Windows it was replaced by snow")
  }
  if ((parallel == "snow" | parallel == "multicore") & missing(ncpus)) 
    stop("You have to specify the number of available processors to parallelize 
           the bootstrap")
  
  if (!inherits(f, "fitdistcens"))
    stop("Use only with 'fitdistcens' objects")
  
  if(!is.null(f$weights))
    stop("Bootstrap is not yet available when using weights")
  
  # non parametric bootstrap
  n <- length(f$censdata[, 1])
  numrow <- seq(1, n)
  rnumrow <- sample(numrow, size=niter*n, replace=TRUE)
  dim(rnumrow) <- c(n, niter)
  start <- as.list(f$estimate) #a named vector is no longer is accepted as starting values.
  if(is.function(f$fix.arg.fun))
    fix.arg <- f$fix.arg.fun
  else 
    fix.arg <- f$fix.arg
  if (is.null(f$dots) && !is.function(fix.arg))
  {
    funcmle <- function(iter) 
    {
      mle <- try(do.call(mledist, list(data=data.frame(left=f$censdata[rnumrow[, iter], ]$left, 
                                                       right=f$censdata[rnumrow[, iter], ]$right), distr=f$distname, start=start, 
                                       fix.arg=fix.arg, checkstartfix=TRUE)), silent=silent)
      if(inherits(mle, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(mle$estimate, mle$convergence))
      
    }
  }else if (is.null(f$dots) && is.function(fix.arg))
  {
    funcmle <- function(iter) 
    {
      bootdata <- data.frame(left=f$censdata[rnumrow[, iter], ]$left, 
                             right=f$censdata[rnumrow[, iter], ]$right)
      fix.arg.iter <- fix.arg(cens2pseudo(bootdata)$pseudo)
      mle <- try(do.call(mledist, list(data=bootdata, distr=f$distname, start=start, 
                                       fix.arg=fix.arg.iter, checkstartfix=TRUE)), silent=silent)
      if(inherits(mle, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(mle$estimate, mle$convergence))
    }
    
  }else if(!is.null(f$dots) && !is.function(fix.arg)) 
  {
    funcmle <- function(iter) 
    {
      mle <- try(do.call(mledist, c(list(data=data.frame(left=f$censdata[rnumrow[, iter], ]$left, 
                                                         right=f$censdata[rnumrow[, iter], ]$right), distr=f$distname, start=start), 
                                    fix.arg=fix.arg, f$dots, checkstartfix=TRUE)), silent=silent)
      if(inherits(mle, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(mle$estimate, mle$convergence))
    }
  }else if(!is.null(f$dots) && is.function(fix.arg)) 
  {
    funcmle <- function(iter) 
    {
      bootdata <- data.frame(left=f$censdata[rnumrow[, iter], ]$left, 
                             right=f$censdata[rnumrow[, iter], ]$right)
      fix.arg.iter <- fix.arg(cens2pseudo(bootdata)$pseudo)
      mle <- try(do.call(mledist, c(list(data=bootdata, distr=f$distname, start=start), 
                                    fix.arg=fix.arg.iter, f$dots, checkstartfix=TRUE)), silent=silent)
      if(inherits(mle, "try-error"))
        return(c(rep(NA, length(start)), 100))
      else
        return(c(mle$estimate, mle$convergence))
    }
  }else
    stop("wrong implementation in bootdistcens")
  
  owarn <- getOption("warn")
  oerr <- getOption("show.error.messages")
  options(warn=ifelse(silent, -1, 0), show.error.messages=!silent)
  
  # parallel or sequential computation
  if (parallel != "no") 
  {
    if (parallel == "snow") type <- "PSOCK"
    else if (parallel == "multicore") type <- "FORK"
    clus <- parallel::makeCluster(ncpus, type = type)
    resboot <- parallel::parSapply(clus, 1:niter, funcmle)
    parallel::stopCluster(clus)
  }
  else
  {
    resboot <- sapply(1:niter, funcmle)
  }
  
  options(warn=owarn, show.error.messages=oerr)
  
  rownames(resboot) <- c(names(start), "convergence")
  if (length(resboot[, 1])>2) {
    estim <- data.frame(t(resboot)[, -length(resboot[, 1])])
    bootCI <- cbind(apply(resboot[-length(resboot[, 1]), ], 1, median, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.025, na.rm=TRUE), 
                    apply(resboot[-length(resboot[, 1]), ], 1, quantile, 0.975, na.rm=TRUE))
    colnames(bootCI) <- c("Median", "2.5%", "97.5%")
  }
  else {
    estim <- as.data.frame(t(resboot)[, -length(resboot[, 1])])
    names(estim) <- names(f$estimate)
    bootCI <- c(median(resboot[-length(resboot[, 1]), ], na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.025, na.rm=TRUE), 
                quantile(resboot[-length(resboot[, 1]), ], 0.975, na.rm=TRUE)) 
    names(bootCI) <- c("Median", "2.5%", "97.5%")
  }       
  res <- structure(list(estim=estim, converg=t(resboot)[, length(resboot[, 1])], 
                        method="nonparam", nbboot=niter, CI=bootCI, fitpart=f), 
                   class="bootdistcens")
  res
}

print.bootdistcens <- function(x, ...){
  if (!inherits(x, "bootdistcens"))
    stop("Use only with 'bootdistcens' objects")
  cat("Parameter values obtained with nonparametric bootstrap \n")
  print(x$estim, ...)    
  nconverg <- length(x$converg[x$converg==0])
  if (nconverg < length(x$converg))
  {
    cat("\n")
    cat("The estimation method converged only for", nconverg, "among", 
        length(x$converg), "iterations \n")
  }
  
}

plot.bootdistcens <- function(x, ...){
  if (!inherits(x, "bootdistcens"))
    stop("Use only with 'bootdistcens' objects")
  if (dim(x$estim)[2]==1) {
    stripchart(x$estim, method="jitter", 
               xlab="Bootstrapped values of the parameter", ...)
  }
  else {
    if (dim(x$estim)[2]==2)
      plot(x$estim, 
           main="Bootstrapped values of the two parameters", ...)
    else 
      plot(x$estim, 
           main="Bootstrapped values of parameters", ...)
  }
}

summary.bootdistcens <- function(object, ...){
  if (!inherits(object, "bootdistcens"))
    stop("Use only with 'bootdistcens' objects")
  
  class(object) <- c("summary.bootdistcens", class(object))  
  object
  
}


print.summary.bootdistcens <- function(x, ...){
  if (!inherits(x, "summary.bootdistcens"))
    stop("Use only with 'summary.bootdistcens' objects")
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

density.bootdistcens <- function(..., bw = "nrd0", adjust = 1, kernel = "gaussian")
{
  x <- list(...)
  if(inherits(x, "bootdistcens"))
  {
    x <- list(x)
  }else if(!is.list(x))
  {
    stop("argument x must be a list of 'bootdistcens' objects")
  }else
  {
    if(any(sapply(x, function(y) !inherits(y, "bootdistcens"))))        
      stop("argument x must be a list of 'bootdistcens' objects")
  }
  print(str(x))
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
            class="density.bootdistcens")
}


plot.density.bootdistcens <- function(x, mar=c(4,4,2,1), lty=NULL, col=NULL, lwd=NULL, ...)
{
  if (!inherits(x, "density.bootdistcens"))
    stop("Use only with 'density.bootdistcens' objects")
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

print.density.bootdistcens <- function(x, ...)
{
  if (!inherits(x, "density.bootdistcens"))
    stop("Use only with 'density.bootdistcens' objects")
  
  nbboot <- paste(attr(x, "nbboot"), collapse=", ")
  cat("\nBootstrap values for: ", attr(x, "distname"), " for ",
      attr(x, "nbobject"), " object(s) with ", 
      nbboot, " bootstrap values (original sample size ",
      attr(x, "n"), ").", sep = "")
  invisible(x)
}
