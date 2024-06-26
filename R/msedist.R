#############################################################################
#   Copyright (c) 2019 Christophe Dutang                                                                                                  
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

msedist <- function (data, distr, phidiv="KL", power.phidiv=NULL, start=NULL, 
                     fix.arg=NULL, optim.method="default", lower=-Inf, upper=Inf, custom.optim=NULL, 
                     weights=NULL, silent=TRUE, gradient=NULL, checkstartfix=FALSE, calcvcov=FALSE, ...)
  # data may correspond to a vector for non censored data or to
  # a dataframe of two columns named left and right for censored data 
{
  if (!is.character(distr)) 
    stop("distr must be a character string naming a distribution")
  else 
    distname <- distr
  
  pdistname <- paste("p",distname,sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", pdistname, " function must be defined"))
  
  ddistname <- paste("d",distname,sep="")    
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  argddistname <- names(formals(ddistname))
  
  if(is.null(custom.optim))
    optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
  
  phidiv <- match.arg(phidiv, c("KL", "J", "R", "H", "V"))
  if(phidiv == "KL") #Kullback-Leibler information
    stopifnot(is.null(power.phidiv))
  if(phidiv == "J") #Jeffreys divergence
    stopifnot(is.null(power.phidiv))
  if(phidiv == "R") #Renyi divergence
  {
    stopifnot(length(power.phidiv) == 1)
    stopifnot(power.phidiv > 0 && power.phidiv != 1)
  }
  if(phidiv %in% c("H", "V")) #Hellinger distance or Vajda information
  {
    stopifnot(length(power.phidiv) == 1)
    stopifnot(power.phidiv >= 1)
  } 
  
  start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
  if(is.vector(start.arg)) #backward compatibility
    start.arg <- as.list(start.arg)
  
  if(!is.null(weights))
  {
    if(any(weights < 0))
      stop("weights should be a vector of integers greater than 0")
    if(!isallintw(weights))
      stop("weights should be a vector of (strictly) positive integers")
    if(length(weights) != NROW(data))
      stop("weights should be a vector with a length equal to the observation number")
    warning("weights are not taken into account in the default initial values")
  }
  
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
    stop("not yet implemented")
  }
  
  if (cens) {
    stop("not yet implemented")
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
    if(is.null(weights))
      weights <- rep(1, NROW(data))
    if(length(weights) != NROW(data))
      stop("weights should be a vector with a length equal to the observation number")
    # objective function (to minimize) is the opposite of the weighted average of the log spacings
    # the argument names are:
    # - par for parameters (like in optim function)
    # - fix.arg for optional fixed parameters
    # - obs for observations (previously dat but conflicts with genoud data.type.int argument)
    # - pdistnam for distribution name
    #NB1: by keeping only positive value, we take into cdf that are step functions
    #NB2: weight has a duplicate first value to have same length as spacings diffFx
    
    if(phidiv == "KL")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      {
        sx <- c(-Inf, sort(obs), Inf)
        n <- length(sx)
        Fxi <- do.call(pdistnam, c(list(sx[-1]), as.list(par), as.list(fix.arg)))
        Fxim1 <- do.call(pdistnam, c(list(sx[-n]), as.list(par), as.list(fix.arg)))
        diffFx <- Fxi - Fxim1
        idxPositive <- diffFx > 0
        if(sum(idxPositive, na.rm = TRUE) == 0)
          return(NaN)
        wx <- c(weights[1], weights) 
        - mean(wx[idxPositive] * log(diffFx[idxPositive]))
      }
    }else if(phidiv == "J")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      {
        sx <- c(-Inf, sort(obs), Inf)
        n <- length(sx)
        Fxi <- do.call(pdistnam, c(list(sx[-1]), as.list(par), as.list(fix.arg)))
        Fxim1 <- do.call(pdistnam, c(list(sx[-n]), as.list(par), as.list(fix.arg)))
        diffFx <- Fxi - Fxim1
        idxPositive <- diffFx > 0
        if(sum(idxPositive, na.rm = TRUE) == 0)
          return(NaN)
        wx <- c(weights[1], weights) 
        - mean(wx[idxPositive] * log(diffFx[idxPositive]) * (1-diffFx[idxPositive]))
      }
    }else if(phidiv == "R")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      {
        sx <- c(-Inf, sort(obs), Inf)
        n <- length(sx)
        Fxi <- do.call(pdistnam, c(list(sx[-1]), as.list(par), as.list(fix.arg)))
        Fxim1 <- do.call(pdistnam, c(list(sx[-n]), as.list(par), as.list(fix.arg)))
        diffFx <- Fxi - Fxim1
        idxPositive <- diffFx > 0
        if(sum(idxPositive, na.rm = TRUE) == 0)
          return(NaN)
        wx <- c(weights[1], weights) 
        - mean(wx[idxPositive] * diffFx[idxPositive]^power.phidiv * sign(1-power.phidiv))
      }
    }else if(phidiv == "H")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      {
        sx <- c(-Inf, sort(obs), Inf)
        n <- length(sx)
        Fxi <- do.call(pdistnam, c(list(sx[-1]), as.list(par), as.list(fix.arg)))
        Fxim1 <- do.call(pdistnam, c(list(sx[-n]), as.list(par), as.list(fix.arg)))
        diffFx <- Fxi - Fxim1
        idxPositive <- diffFx > 0
        if(sum(idxPositive, na.rm = TRUE) == 0)
          return(NaN)
        wx <- c(weights[1], weights) 
        mean(wx[idxPositive] * abs( 1-diffFx[idxPositive]^(1/power.phidiv) )^power.phidiv )
      }
    }else if(phidiv == "V")
    {
      fnobj <- function(par, fix.arg, obs, pdistnam)
      {
        sx <- c(-Inf, sort(obs), Inf)
        n <- length(sx)
        Fxi <- do.call(pdistnam, c(list(sx[-1]), as.list(par), as.list(fix.arg)))
        Fxim1 <- do.call(pdistnam, c(list(sx[-n]), as.list(par), as.list(fix.arg)))
        diffFx <- Fxi - Fxim1
        idxPositive <- diffFx > 0
        if(sum(idxPositive, na.rm = TRUE) == 0)
          return(NaN)
        wx <- c(weights[1], weights) 
        mean(wx[idxPositive] * abs( 1-diffFx[idxPositive] )^power.phidiv )
      }
    }else
      stop("error: wrong phidiv")
    
  }
  else # if (!cens) 
    stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
  
  # Function to calculate the loglikelihood to return
  loglik <- function(par, fix.arg, obs, ddistnam) {
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
                                              fix.arg=fix.arg, obs=data, pdistnam=pdistname, 
                                              hessian=!is.null(gradient), method=meth, 
                                              ...), silent=TRUE)
        if(!inherits(opttryerror, "try-error"))
          if(length(opt$counts) == 1) #appears when the initial point is a solution
            opt$counts <- c(opt$counts, NA)
        
        
      }else #opt.fun == "optim"
      {
        opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                        pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, 
                                        upper=upper, ...), silent=TRUE)       
      }
      
    }else #hasbound == FALSE
    {
      opt.fun <- "optim"
      opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                      pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, 
                                      upper=upper, ...), silent=TRUE)       
    }
    options(warn=owarn)
    
    if (inherits(opttryerror,"try-error"))
    {
      warnings("The function optim encountered an error and stopped.")
      if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))          
      return(list(estimate = rep(NA,length(vstart)), convergence = 100, loglik = NA, 
                  hessian = NA))
    }
    
    if (opt$convergence>0)
    {
      warnings("The function optim failed to converge, with the error code ",
               opt$convergence)
    }
    if(is.null(names(opt$par)))
      names(opt$par) <- names(vstart)
    res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                hessian = opt$hessian, optim.function=opt.fun, optim.method=meth, 
                fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights=NULL,
                counts=opt$counts, optim.message=opt$message,
                loglik=loglik(opt$par, fix.arg, data, ddistname),
                phidiv=phidiv, power.phidiv=power.phidiv)
    
  }
  else # Try to minimize the gof distance using a user-supplied optim function 
  {
    options(warn=ifelse(silent, -1, 0))
    if (!cens)
      opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, 
                                             pdistnam=pdistname, par=vstart, ...),
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
    
    if (opt$convergence>0) 
    {
      warnings("The customized optimization function failed to converge, with the error code ",
               opt$convergence)
    }
    if(is.null(names(opt$par)))
      names(opt$par) <- names(vstart)
    argdot <- list(...)
    method.cust <- argdot$method
    res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                hessian = opt$hessian, optim.function=custom.optim, optim.method=method.cust,
                fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights=weights,
                counts=opt$counts, optim.message=opt$message,
                loglik=loglik(opt$par, fix.arg, data, ddistname),
                phidiv=phidiv, power.phidiv=power.phidiv)
  }   
  return(res)                
}
