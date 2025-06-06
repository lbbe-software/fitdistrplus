#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller                                                                                                  
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
### maximum goodness-of-fit estimation for censored or non-censored data
### and continuous distributions
### (at this time only available for non censored data)
###
###         R functions
### 

mgedist <- function (data, distr, gof = "CvM", start=NULL, fix.arg=NULL, optim.method="default",
                     lower=-Inf, upper=Inf, custom.optim=NULL, silent=TRUE, gradient=NULL, 
                     checkstartfix=FALSE, calcvcov=FALSE, ...)
  # data may correspond to a vector for non censored data or to
  # a dataframe of two columns named left and right for censored data 
{
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else 
    distname <- distr
  
  if (is.element(distname,c("binom","nbinom","geom","hyper","pois"))) 
    stop("Maximum goodness-of-fit estimation method is not intended to fit discrete distributions")
  
  pdistname <- paste("p",distname,sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", pdistname, " function must be defined"))
  
  ddistname <- paste("d",distname,sep="")    
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  argddistname <- names(formals(ddistname))
  
  if(is.null(custom.optim))
    optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
  
  gof <- match.arg(gof, c("CvM", "KS", "AD", "ADR", "ADL", "AD2R", "AD2L", "AD2"))
  
  start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
  if(is.vector(start.arg)) #backward compatibility
    start.arg <- as.list(start.arg)
  
  my3dots <- list(...)
  if ("weights" %in% names(my3dots))
    stop("Weights is not allowed for maximum GOF estimation")
  
  if (is.vector(data)) {
    cens <- FALSE
    if (!(is.numeric(data) & length(data)>1)) 
      stop("data must be a numeric vector of length greater than 1 for non censored data
            or a dataframe with two columns named left and right and more than one line for censored data")
    checkUncensoredNAInfNan(data)
  }
  else {
    cens <- TRUE
    censdata <- data
    if (!(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[,1])>1))
      stop("data must be a numeric vector of length greater than 1 for non censored data
        or a dataframe with two columns named left and right and more than one line for censored data")
    pdistname<-paste("p",distname,sep="")
    if (!exists(pdistname,mode="function"))
      stop(paste("The ",pdistname," function must be defined to apply maximum likelihood to censored data"))
    
  }
  
  if (cens) {
    # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector,
    #   icens (interval censored) = dataframe with left and right 
    # and ncens (not censored) = vector
    lcens<-censdata[is.na(censdata$left),]$right
    if (any(is.na(lcens)) )
      stop("An observation cannot be both right and left censored, coded with two NA values")
    rcens<-censdata[is.na(censdata$right),]$left
    ncens<-censdata[censdata$left==censdata$right & !is.na(censdata$left) & 
                      !is.na(censdata$right),]$left
    icens<-censdata[censdata$left!=censdata$right & !is.na(censdata$left) & 
                      !is.na(censdata$right),]
    # Definition of a data set for calculation of starting values
    data<-c(rcens,lcens,ncens,(icens$left+icens$right)/2)
  }
  
  if(!checkstartfix) #pre-check has not been done by fitdist() or bootdist()
  {
    # manage starting/fixed values: may raise errors or return two named list
    arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=data, 
                                distname=distname)
    
    #check inconsistent parameters
    hasnodefaultval <- sapply(formals(ddistname), is.name)
    arg_startfix <- checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg, 
                                   argddistname, hasnodefaultval)
    #arg_startfix contains two names list (no longer NULL nor function)  
    
    #set fix.arg.fun
    if(is.function(fix.arg))
      fix.arg.fun <- fix.arg
    else
      fix.arg.fun <- NULL
  }else #pre-check has been done by fitdist() or bootdist()
  {
    arg_startfix <- list(start.arg=start, fix.arg=fix.arg)
    fix.arg.fun <- NULL
  }
  
  #unlist starting values as needed in optim()
  vstart <- unlist(arg_startfix$start.arg)
  #sanity check
  if(is.null(vstart))
    stop("Starting values could not be NULL with checkstartfix=TRUE")
  
  #erase user value
  #(cannot coerce to vector as there might be different modes: numeric, character...)
  fix.arg <- arg_startfix$fix.arg
  
  ############# MGE fit using optim or custom.optim ##########
  
  # definition of the function to minimize depending on the argument gof
  # for non censored data
  if (!cens) 
  {
    # the argument names are:
    # - par for parameters (like in optim function)
    # - fix.arg for optional fixed parameters
    # - obs for observations (previously dat but conflicts with genoud data.type.int argument)
    # - pdistnam for distribution name
    if (gof == "CvM")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        1/(12*n^2) + mean( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
      }
    }else if (gof == "KS")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam) 
      {
        n <- length(obs)
        s <- sort(obs)
        obspu <- seq(1,n)/n
        obspl <- seq(0,n-1)/n
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        max(pmax(abs(theop-obspu), abs(theop-obspl)))
      }
    }else if (gof == "AD")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        ilogpi <- log(theop * (1 - rev(theop))) * (2 * 1:n - 1)
        idx <- is.finite(ilogpi)
        ifelse(sum(idx) > 0, - sum(idx)/n - mean( ilogpi[idx] )/n, .Machine$integer.max)
      }
    }else if (gof == "ADR")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        ilogpi <- log(1 - rev(theop)) * (2 * 1:n - 1)
        idx <- is.finite(ilogpi)
        ifelse(sum(idx) > 0, sum(idx)/2/n - 2 * sum(theop[idx])/n - mean ( ilogpi[idx] )/n, .Machine$integer.max)
      }
    }else if (gof == "ADL")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        ilogpi <- (2 * 1:n - 1) * log(theop)
        idx <- is.finite(ilogpi)
        ifelse(sum(idx) > 0, -3*sum(idx)/2/n + 2 * sum(theop[idx])/n - mean ( ilogpi[idx] )/n, .Machine$integer.max)
      }
    }else if (gof == "AD2R")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        logpi <- log(1 - theop)
        i1pi2 <- (2 * 1:n - 1) / (1 - rev(theop))
        idx <- is.finite(logpi) & is.finite(i1pi2)
        ifelse(sum(idx) > 0, 2 * sum(logpi[idx])/n + mean( i1pi2[idx] )/n, .Machine$integer.max)
      }
    }else if (gof == "AD2L")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        logpi <- log(theop)
        i1pi <- (2 * 1:n - 1) / theop 
        idx <- is.finite(logpi) & is.finite(i1pi)
        ifelse(sum(idx) > 0, 2 * sum( logpi[idx] )/n + mean ( i1pi[idx] )/n, .Machine$integer.max)
      }
    }else if (gof == "AD2")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      { 
        n <- length(obs)
        s <- sort(obs)
        theop <- do.call(pdistnam, c(list(s), as.list(par), as.list(fix.arg)))
        logpi <- log(theop * (1 - theop)) 
        i1pi <- (2 * 1:n - 1) / theop
        i1pi2 <- (2 * 1:n - 1) / (1 - rev(theop))
        idx <- is.finite(logpi) & is.finite(i1pi) & is.finite(i1pi2)
        ifelse(sum(idx) > 0, 2 * sum( logpi[idx] )/n + mean ( i1pi[idx]  + i1pi2[idx] )/n, .Machine$integer.max)
      }
    }else
      stop("wrong gof metric.")
  }else # if (!cens) 
    stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
  
  # Function to calculate the loglikelihood to return
  loglik <- function(par, fix.arg, obs, ddistnam) 
  {
    sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
  }
  
  owarn <- getOption("warn")        
  
  # Try to minimize the gof distance using the base R optim function
  if(is.null(custom.optim))
  {
    hasbound <- any(is.finite(lower) | is.finite(upper))
    
    # Choice of the optimization method  
    if (optim.method == "default")
    {
      meth <- ifelse(length(vstart) > 1, "Nelder-Mead", "BFGS") 
    }else
      meth <- optim.method
    
    if(meth == "BFGS" && hasbound && is.null(gradient))
    {
      meth <- "L-BFGS-B"
      txt1 <- "The BFGS method cannot be used with bounds without provided the gradient."
      txt2 <- "The method is changed to L-BFGS-B."
      warning(paste(txt1, txt2))
    }
    
    options(warn=ifelse(silent, -1, 0))
    #select optim or constrOptim
    if(hasbound) #finite bounds are provided
    {
      if(!is.null(gradient))
      {
        opt.fun <- "constrOptim"
      }else #gradient == NULL
      {
        if(meth == "Nelder-Mead")
          opt.fun <- "constrOptim"
        else if(meth %in% c("L-BFGS-B", "Brent"))
          opt.fun <- "optim"
        else
        {
          txt1 <- paste("The method", meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
          txt2 <- "Only optimization methods L-BFGS-B, Brent and Nelder-Mead can be used in such case."
          stop(paste(txt1, txt2))
        }
      }
      if(opt.fun == "constrOptim")
      {
        #recycle parameters
        npar <- length(vstart) #as in optim() line 34
        lower <- as.double(rep_len(lower, npar)) #as in optim() line 64
        upper <- as.double(rep_len(upper, npar))
        
        # constraints are : Mat %*% theta >= Bnd, i.e. 
        # +1 * theta[i] >= lower[i]; 
        # -1 * theta[i] >= -upper[i]
        
        #select rows from the identity matrix
        haslow <- is.finite(lower)
        Mat <- diag(npar)[haslow, ]
        #select rows from the opposite of the identity matrix
        hasupp <- is.finite(upper)
        Mat <- rbind(Mat, -diag(npar)[hasupp, ])
        colnames(Mat) <- names(vstart)
        rownames(Mat) <- paste0("constr", 1:NROW(Mat))
        
        #select the bounds
        Bnd <- c(lower[is.finite(lower)], -upper[is.finite(upper)])
        names(Bnd) <- paste0("constr", 1:length(Bnd))
        
        initconstr <- Mat %*% vstart - Bnd
        if(any(initconstr < 0))
          stop("Starting values must be in the feasible region.")
        
        opttryerror <- try(opt <- constrOptim(theta=vstart, f=fnobj, ui=Mat, ci=Bnd, grad=gradient,
                                              fix.arg=fix.arg, obs=data, pdistnam=pdistname, hessian=!is.null(gradient), method=meth, 
                                              ...), silent=TRUE)
        if(!inherits(opttryerror, "try-error"))
          if(length(opt$counts) == 1) #appears when the initial point is a solution
            opt$counts <- c(opt$counts, NA)
        
        
      }else #opt.fun == "optim"
      {
        opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                        pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                                        ...), silent=TRUE)       
      }
      
    }else #hasbound == FALSE
    {
      opt.fun <- "optim"
      opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                      pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                                      ...), silent=TRUE)       
    }
    options(warn=owarn)
    
    if (inherits(opttryerror,"try-error"))
    {
      warnings("The function optim encountered an error and stopped.")
      if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))          
      return(list(estimate = rep(NA,length(vstart)), convergence = 100, loglik = NA, 
                  hessian = NA))
    }
    
    if (opt$convergence>0) {
      warnings("The function optim failed to converge, with the error code ",
               opt$convergence)
    }
    if(is.null(names(opt$par)))
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                hessian = opt$hessian, optim.function=opt.fun, optim.method=meth, 
                fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights=NULL,
                counts=opt$counts, optim.message=opt$message,
                loglik=loglik(opt$par, fix.arg, data, ddistname), gof=gof)
  }else # Try to minimize the gof distance using a user-supplied optim function 
  {
    options(warn=ifelse(silent, -1, 0))
    if (!cens)
      opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, pdistnam=pdistname, par=vstart, ...),
                         silent=TRUE)
    else
      stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
    options(warn=owarn)
    
    if (inherits(opttryerror,"try-error"))
    {
      warnings("The customized optimization function encountered an error and stopped.")
      if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))          
      return(list(estimate = rep(NA,length(vstart)), convergence = 100, value = NA, 
                  hessian = NA))
    }
    
    if (opt$convergence>0) {
      warnings("The customized optimization function failed to converge, with the error code ",
               opt$convergence)
    }
    if(is.null(names(opt$par)))
      names(opt$par) <- names(vstart)
    argdot <- list(...)
    method.cust <- argdot$method
    res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                hessian = opt$hessian, optim.function=custom.optim, optim.method=method.cust,
                fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights=NULL,
                counts=opt$counts, optim.message=opt$message,
                loglik=loglik(opt$par, fix.arg, data, ddistname), gof=gof)
  }   
  return(res)                
}
