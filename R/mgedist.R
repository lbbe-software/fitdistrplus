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
    lower=-Inf, upper=Inf, custom.optim=NULL, ...)
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
    
    # MGE fit 
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
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
            }
        else     
        if (gof == "KS")
            fnobj <- function(par, fix.arg, obs, pdistnam) 
            {
                n <- length(obs)
                s <- sort(obs)
                obspu <- seq(1,n)/n
                obspl <- seq(0,n-1)/n
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                max(pmax(abs(theop-obspu),abs(theop-obspl)))
            }
        else
        if (gof == "AD")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ) 
            }
        else
        if (gof == "ADR")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                n/2 - 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(1 - rev(theop)) )
            }
        else
        if (gof == "ADL")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                -3*n/2 + 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(theop) )
            }
        else  
        if (gof == "AD2R")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                2 * sum(log(1 - theop)) + mean ( (2 * 1:n - 1) / (1 - rev(theop)) )
            }
        else  
        if (gof == "AD2L")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                2 * sum(log(theop)) + mean ( (2 * 1:n - 1) / theop )
            }
         else  
        if (gof == "AD2")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                2 * sum(log(theop) + log(1 - theop) ) + 
                mean ( ((2 * 1:n - 1) / theop) + ((2 * 1:n - 1) / (1 - rev(theop))) )
            }
    }
    else # if (!cens) 
        stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
        
    # Function to calculate the loglikelihood to return
    loglik <- function(par, fix.arg, obs, ddistnam) {
        sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }
    
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

        
    # Try to minimize the gof distance using the base R optim function
    if(is.null(custom.optim))
    {
        if (!cens)
            opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, pdistnam=pdistname,
            hessian=TRUE, method=meth, lower=lower, upper=upper, ...), silent=TRUE)        
        else 
            stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
                
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
                    hessian = opt$hessian, gof=gof, optim.function="optim",
                    loglik=loglik(opt$par, fix.arg, data, ddistname), fix.arg = fix.arg, 
                    optim.method=meth, fix.arg.fun = fix.arg.fun)
    }
    else # Try to minimize the gof distance using a user-supplied optim function 
    {
        if (!cens)
            opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, pdistnam=pdistname, par=vstart, ...),
            silent=TRUE)
        else
            stop("Maximum goodness-of-fit estimation is not yet available for censored data.")
        
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
        res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                    gof=gof, hessian = opt$hessian, optim.function=custom.optim,
                    loglik=loglik(opt$par, fix.arg, data, ddistname), fix.arg = fix.arg,
                    optim.method=NULL, fix.arg.fun = fix.arg.fun)
    }   
   
    return(res)                
}
