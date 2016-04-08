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
    lower=-Inf, upper=Inf, custom.optim=NULL, weights=NULL, silent=TRUE, ...)
    # data may correspond to a vector for non censored data or to
    # a dataframe of two columns named left and right for censored data 
{
    if (!is.character(distr)) 
        stop("distr must be a character string naming a distribution")
    else 
        distname <- distr
    ddistname <- paste("d", distname, sep="")
    
    if (!exists(ddistname, mode="function"))
        stop(paste("The ", ddistname, " function must be defined"))

    start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
    if(is.vector(start.arg)) #backward compatibility
      start.arg <- as.list(start.arg)
    
    txt1 <- "data must be a numeric vector of length greater than 1 for non censored data"
    txt2 <- "or a dataframe with two columns named left and right and more than one line for censored data"
    if(!is.null(weights))
    {
      if(any(weights < 0))
        stop("weights should be a vector of numerics greater than 0")
      if(length(weights) != NROW(data))
        stop("weights should be a vector with a length equal to the observation number")
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
        # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector, 
        #   icens (interval censored) = dataframe with left and right 
        # and ncens (not censored) = vector
        lcens<-censdata[is.na(censdata$left), ]$right
        if (any(is.na(lcens)) )
            stop("An observation cannot be both right and left censored, coded with two NA values")
        rcens<-censdata[is.na(censdata$right), ]$left
        ncens<-censdata[censdata$left==censdata$right & !is.na(censdata$left) & 
            !is.na(censdata$right), ]$left
        icens<-censdata[censdata$left!=censdata$right & !is.na(censdata$left) & 
            !is.na(censdata$right), ]
        # Definition of a data set for calculation of starting values
        data<-c(rcens, lcens, ncens, (icens$left+icens$right)/2)
    }
    
    # definition of starting/fixed values values
    argddistname <- names(formals(ddistname))
    chfixstt <- checkparam(start.arg=start.arg, fix.arg=fix.arg, argdistname=argddistname, 
                           errtxt=NULL, data10=head(data, 10), distname=distname)
    if(!chfixstt$ok)
      stop(chfixstt$txt)
    #unlist starting values as needed in optim()
    if(is.function(chfixstt$start.arg))
      vstart <- chfixstt$start.arg(data)
    else
      vstart <- unlist(chfixstt$start.arg)
    if(is.function(fix.arg)) #function
    { 
      fix.arg.fun <- fix.arg
      fix.arg <- fix.arg(data)
    }else
      fix.arg.fun <- NULL
    #otherwise fix.arg is a named list or NULL
    
    # end of the definition of starting/fixed values   
    
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
        -sum(weights * p1) - sum(weights * p2) - sum(weights * p3) - sum(weights * p4) 
      }
    }else
        stop("not yet implemented.")
   
    # Choice of the optimization method    
    if (optim.method == "default")
    {
        if(is.infinite(lower) && is.infinite(upper))
        { 
            if (length(vstart) > 1) 
                meth <- "Nelder-Mead"
            else 
                meth <- "BFGS"
        }else
            meth <- "L-BFGS-B"
    }else
        meth <- optim.method
        
    owarn <- getOption("warn")
    
    # Try to minimize the minus (log-)likelihood using the base R optim function
    if(is.null(custom.optim))
    {
        options(warn=ifelse(silent, -1, 0))
        if (!cens)
            opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, 
                ddistnam=ddistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                ...), silent=TRUE)        
        else 
            opttryerror <- try(opt <- optim(par=vstart, fn=fnobjcens, fix.arg=fix.arg, 
                rcens=rcens, lcens=lcens, icens=icens, ncens=ncens, ddistnam=ddistname, 
                pdistnam=pdistname, hessian=TRUE, method=meth, lower=lower, upper=upper, 
                ...), silent=TRUE)   
        options(warn=owarn)
        
        if (inherits(opttryerror, "try-error"))
        {
            warnings("The function optim encountered an error and stopped.")
            if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))          
            return(list(estimate = rep(NA, length(vstart)), convergence = 100, loglik = NA, 
                        hessian = NA, optim.function="optim", fix.arg = fix.arg, 
                        optim.method=meth, fix.arg.fun = fix.arg.fun, counts=c(NA, NA)))
        }
        
        if (opt$convergence>0) {
            warnings("The function optim failed to converge, with the error code ", 
                     opt$convergence)
        }
        if(is.null(names(opt$par)))
          names(opt$par) <- names(vstart)
        res <- list(estimate = opt$par, convergence = opt$convergence, loglik = -opt$value, 
                    hessian = opt$hessian, optim.function="optim", fix.arg = fix.arg, 
                    optim.method=meth, fix.arg.fun = fix.arg.fun, weights = weights, 
                    counts=opt$counts, optim.message=opt$message)
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
        argdot <- list(...)
        method.cust <- argdot[argdot == "method"]
        if(length(method.cust) == 0)
        {
          method.cust <- NULL
        }
        if(is.null(names(opt$par)))
          names(opt$par) <- names(vstart)
        res <- list(estimate = opt$par, convergence = opt$convergence, loglik = -opt$value, 
                      hessian = opt$hessian, optim.function = custom.optim, fix.arg = fix.arg,
                      method = method.cust, fix.arg.fun = fix.arg.fun, weights = weights, 
                      counts=opt$counts, optim.message=opt$message)        
    }   
        
    return(res) 
}

## old function with previous name for censored data
mledistcens <- function(censdata, distr, start=NULL, optim.method="default", lower=-Inf, upper=Inf)
{
    stop("The function \"mledistcens\" is no more used. Now the same function \"mledist\" must be used both for censored and non censored data.")
}
