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
### maximum likelihood estimation for censored or non-censored data
###
###         R functions
### 
### many ideas are taken from the fitdistr function of the MASS package and 
### the mle function of the stat package.

mledist <- function (data, distr, start=NULL, fix.arg=NULL, optim.method="default", 
    lower=-Inf, upper=Inf, custom.optim=NULL, weights=NULL, silent=TRUE, gradient=NULL, 
    checkstartfix=FALSE, ...)
    # data may correspond to a vector for non censored data or to
    # a dataframe of two columns named left and right for censored data 
{
    if (!is.character(distr)) 
        stop("distr must be a character string naming a distribution")
    else 
        distname <- distr
    ddistname <- paste("d", distname, sep="")
    argddistname <- names(formals(ddistname))
    
    if (!exists(ddistname, mode="function"))
        stop(paste("The ", ddistname, " function must be defined"))
    if(is.null(custom.optim))
      optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))

    start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
    if(is.vector(start.arg)) #backward compatibility
      start.arg <- as.list(start.arg)
    
    txt1 <- "data must be a numeric vector of length greater than 1 for non censored data"
    txt2 <- "or a dataframe with two columns named left and right and more than one line for censored data"
    if(!is.null(weights))
    {
      if(any(weights < 0))
        stop("weights should be a vector of integers greater than 0")
      if(!is.allint.w(weights))
        stop("weights should be a vector of (strictly) positive integers")
      if(length(weights) != NROW(data))
        stop("weights should be a vector with a length equal to the observation number")
      warning("weights are not taken into account in the default initial values")
    }
    
    if (is.vector(data)) {
        cens <- FALSE
        if (!(is.numeric(data) & length(data)>1)) 
            stop(paste(txt1, txt2))
    }
    else {
        cens <- TRUE
        censdata <- data
        if (!(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[, 1])>1))
            stop(paste(txt1, txt2))
        pdistname<-paste("p", distname, sep="")
        if (!exists(pdistname, mode="function"))
            stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))

    }
    
    if (cens) {
        #format data for calculation of starting values and fitting process
        dataformat <- cens2pseudo(censdata)
        data <- dataformat$pseudo
        rcens <- dataformat$rcens; lcens <- dataformat$lcens 
        icens <- dataformat$icens; ncens <- dataformat$ncens
        
        irow <- cens2idxrow(censdata)
        irow.rcens <- irow$rcens; irow.lcens <- irow$lcens
        irow.icens <- irow$icens; irow.ncens <- irow$ncens
    }
    
    if(!checkstartfix) #pre-check has not been done by fitdist() or bootdist()
    {
      # manage starting/fixed values: may raise errors or return two named list
      arg_startfix <- manageparam(start.arg=start, fix.arg=fix.arg, obs=data, 
                                  distname=distname)
      
      #check inconsistent parameters
      arg_startfix <- checkparamlist(arg_startfix$start.arg, arg_startfix$fix.arg, argddistname)
      #arg_startfix contains two names list (no longer NULL nor function)  
      
      #set fix.arg.fun
      if(is.function(fix.arg))
        fix.arg.fun <- fix.arg
      else
        fix.arg.fun <- NULL
    }else #pre-check has been done by fitdist<cens>() or bootdist<cens>()
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
    fix.arg <- unlist(arg_startfix$fix.arg)
    
    
    ############# closed-form formula for uniform distribution ##########
    if(distname == "unif")
    {
        par <- c(min=min(data), max=max(data))
        res <- list(estimate = par[!names(par) %in% names(fix.arg)], convergence = 0, loglik = NA, 
                    hessian = NA, optim.function= NA, fix.arg = fix.arg)
        return(res)
    }
    
   
    ############# MLE fit using optim or custom.optim ##########

    # definition of the function to minimize : - log likelihood
    # for non censored data
    if (!cens && is.null(weights)) {
        # the argument names are:
        # - par for parameters (like in optim function)
        # - fix.arg for optional fixed parameters
        # - obs for observations (previously dat but conflicts with genoud data.type.int argument)
        # - ddistnam for distribution name
        if ("log" %in% argddistname){
            fnobj <- function(par, fix.arg, obs, ddistnam){
                -sum(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg), log=TRUE) ) )
            }
        }
        else{
        fnobj <- function(par, fix.arg, obs, ddistnam) {
            -sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
            }
        }
    }
    else if(cens && is.null(weights)) #censored data
    {    
      argpdistname<-names(formals(pdistname))
        if (("log" %in% argddistname) & ("log.p" %in% argpdistname))
            fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam)
                -sum(do.call(ddistnam, c(list(x=ncens), as.list(par), as.list(fix.arg), list(log=TRUE)))) -
                sum(do.call(pdistnam, c(list(q=lcens), as.list(par), as.list(fix.arg), list(log=TRUE)))) -
                sum(do.call(pdistnam, c(list(q=rcens), as.list(par), as.list(fix.arg), list(lower.tail=FALSE), list(log=TRUE)))) -
                sum(log(do.call(pdistnam, c(list(q=icens$right), as.list(par), as.list(fix.arg))) - # without log=TRUE here
                do.call(pdistnam, c(list(q=icens$left), as.list(par), as.list(fix.arg))) )) # without log=TRUE here
        else
            fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam)
                -sum(log(do.call(ddistnam, c(list(x=ncens), as.list(par), as.list(fix.arg))))) -
                sum(log(do.call(pdistnam, c(list(q=lcens), as.list(par), as.list(fix.arg))))) -
                sum(log(1-do.call(pdistnam, c(list(q=rcens), as.list(par), as.list(fix.arg))))) -
                sum(log(do.call(pdistnam, c(list(q=icens$right), as.list(par), as.list(fix.arg))) - 
                do.call(pdistnam, c(list(q=icens$left), as.list(par), as.list(fix.arg))) ))
    }else if(!cens && !is.null(weights))
    {
        fnobj <- function(par, fix.arg, obs, ddistnam) {
          -sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
        }
    }else if(cens && !is.null(weights))
    {
      fnobjcens <- function(par, fix.arg, rcens, lcens, icens, ncens, ddistnam, pdistnam)
      {
        p1 <- log(do.call(ddistnam, c(list(x=ncens), as.list(par), as.list(fix.arg))))
        p2 <- log(do.call(pdistnam, c(list(q=lcens), as.list(par), as.list(fix.arg)))) 
        p3 <- log(1-do.call(pdistnam, c(list(q=rcens), as.list(par), as.list(fix.arg))))
        p4 <- log(do.call(pdistnam, c(list(q=icens$right), as.list(par), as.list(fix.arg))) - 
                    do.call(pdistnam, c(list(q=icens$left), as.list(par), as.list(fix.arg))) )
        -sum(weights[irow.ncens] * p1) - 
          sum(weights[irow.lcens] * p2) - 
          sum(weights[irow.rcens] * p3) - 
          sum(weights[irow.icens] * p4) 
      }
    }
   
    #get warning value    
    owarn <- getOption("warn")
    
    # Try to minimize the minus (log-)likelihood using the base R optim function
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
            
            if(!cens)
            {
              opttryerror <- try(opt <- constrOptim(theta=vstart, f=fnobj, ui=Mat, ci=Bnd, grad=gradient,
                    fix.arg=fix.arg, obs=data, ddistnam=ddistname, hessian=!is.null(gradient), method=meth, 
                    ...), silent=TRUE)
            }
            else #cens == TRUE
              opttryerror <- try(opt <- constrOptim(theta=vstart, f=fnobjcens, ui=Mat, ci=Bnd, grad=gradient,
                    ddistnam=ddistname, rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, pdistnam=pdistname,
                    fix.arg=fix.arg, hessian=!is.null(gradient), method=meth, ...), silent=TRUE)
            if(!inherits(opttryerror, "try-error"))
              if(length(opt$counts) == 1) #appears when the initial point is a solution
                opt$counts <- c(opt$counts, NA)
            
            
          }else #opt.fun == "optim"
          {
            if(!cens)
              opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                              ddistnam=ddistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                                              ...), silent=TRUE)       
            else #cens == TRUE
              opttryerror <- try(opt <- optim(par=vstart, fn=fnobjcens, fix.arg=fix.arg, gr=gradient,
                                              rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname, 
                                              pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                                              ...), silent=TRUE)   
          }
          
        }else #hasbound == FALSE
        {
          opt.fun <- "optim"
          if(!cens)
            opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                            ddistnam=ddistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                                            ...), silent=TRUE)       
          else #cens == TRUE
            opttryerror <- try(opt <- optim(par=vstart, fn=fnobjcens, fix.arg=fix.arg, gr=gradient,
                                            rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname, 
                                            pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                                            ...), silent=TRUE) 
        }
        options(warn=owarn)
        
        if (inherits(opttryerror, "try-error"))
        {
            warnings("The function optim encountered an error and stopped.")
            if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))          
            return(list(estimate = rep(NA, length(vstart)), convergence = 100, loglik = NA, 
                        hessian = NA, optim.function=opt.fun, fix.arg = fix.arg, 
                        optim.method=meth, fix.arg.fun = fix.arg.fun, counts=c(NA, NA)))
        }
        
        if (opt$convergence>0) {
            warnings("The function optim failed to converge, with the error code ", 
                     opt$convergence)
        }
        if(is.null(names(opt$par)))
          names(opt$par) <- names(vstart)
        res <- list(estimate = opt$par, convergence = opt$convergence, value=opt$value,  
                    hessian = opt$hessian, optim.function=opt.fun, optim.method=meth, 
                    fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights = weights, 
                    counts=opt$counts, optim.message=opt$message, loglik = -opt$value)
    }
    else # Try to minimize the minus (log-)likelihood using a user-supplied optim function 
    {
        options(warn=ifelse(silent, -1, 0))
        if (!cens)
            opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, 
                ddistnam=ddistname, par=vstart, ...), silent=TRUE)
        else
            opttryerror <-try(opt<-custom.optim(fn=fnobjcens, fix.arg=fix.arg, rcens=rcens, 
                lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname, pdistnam=pdistname, 
                par=vstart, ...), silent=TRUE)              
        options(warn=owarn)
        
        if (inherits(opttryerror, "try-error"))
        {
            warnings("The customized optimization function encountered an error and stopped.")
            if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))          
            return(list(estimate = rep(NA, length(vstart)), convergence = 100, loglik = NA, 
                        hessian = NA, optim.function=custom.optim, fix.arg = fix.arg, 
                        fix.arg.fun = fix.arg.fun, counts=c(NA, NA)))
        }
        
        if (opt$convergence>0) {
            warnings("The customized optimization function failed to converge, with the error code ", 
                     opt$convergence)
        }
        if(is.null(names(opt$par)))
          names(opt$par) <- names(vstart)
        argdot <- list(...)
        method.cust <- argdot$method
        res <- list(estimate = opt$par, convergence = opt$convergence, value=opt$value, 
                    hessian = opt$hessian, optim.function = custom.optim, optim.method = method.cust, 
                    fix.arg = fix.arg, fix.arg.fun = fix.arg.fun, weights = weights, 
                    counts=opt$counts, optim.message=opt$message, loglik = -opt$value)        
    }   
        
    return(res) 
}

## old function with previous name for censored data
mledistcens <- function(censdata, distr, start=NULL, optim.method="default", lower=-Inf, upper=Inf)
{
    stop("The function \"mledistcens\" is no more used. Now the same function \"mledist\" must be used both for censored and non censored data.")
}
