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
#        distname <- substring(as.character(match.call()$distr), 2)
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
   
    if (!is.null(fix.arg) & is.null(start))
        stop("Starting values must be defined when some distribution parameters are fixed")
    
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
    # definition of starting values if not previously defined
    if (is.null(start)) {
        if (distname == "norm") {
            n <- length(data)
            sd0 <- sqrt((n - 1)/n) * sd(data)
            mx <- mean(data)
            start <- list(mean=mx, sd=sd0)
        }
        if (distname == "lnorm") {
            if (any(data <= 0)) 
                stop("values must be positive to fit a lognormal distribution")
            n <- length(data)
            ldata <- log(data)
            sd0 <- sqrt((n - 1)/n) * sd(ldata)
            ml <- mean(ldata)
            start <- list(meanlog=ml, sdlog=sd0)
        }
        if (distname == "exp") {
            start <- list(rate=1/mean(data))
        }
        if (distname == "gamma") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            start <- list(shape=m^2/v,rate=m/v)
        }
        if (distname == "beta") {
            if (any(data < 0) | any(data > 1)) 
                stop("values must be in [0-1] to fit a beta distribution")
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            aux <- m*(1-m)/v - 1
            start <- list(shape1=m*aux,shape2=(1-m)*aux)
        }
        if (distname == "weibull") {
            m <- mean(log(data))
            v <- var(log(data))
            shape <- 1.2/sqrt(v)
            scale <- exp(m + 0.572/shape)
            start <- list(shape = shape, scale = scale)
        }
        if (distname == "logis") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            start <- list(location=m,scale=sqrt(3*v)/pi)
        }
        if (distname == "cauchy") {
            start <- list(location=median(data),scale=IQR(data)/2)
        }
        if (distname == "unif") {
            start <- list(min=min(data),max=max(data))
        }
        if (!is.list(start)) 
            stop("'start' must be defined as a named list for this distribution") 
   } # end of the definition of starting values 
   
   ############# MGE fit using optim or custom.optim ##########
    vstart <- unlist(start)
    vfix.arg <- unlist(fix.arg)
    # check of the names of the arguments of the density function
    argddistname <- names(formals(ddistname))   
    m <- match(names(start), argddistname)
    mfix <- match(names(vfix.arg), argddistname)
    if (any(is.na(m)))
        stop("'start' must specify names which are arguments to 'distr'")
    if (any(is.na(mfix)))
        stop("'fix.arg' must specify names which are arguments to 'distr'")
    # check that some parameters are not both in fix.arg and start
    minter <- match(names(start), names(fix.arg))
    if (any(!is.na(minter)))
        stop("a distribution parameter cannot be specified both in 'start' and 'fix.arg'")

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
                1/(12*n) + sum( ( theop - (2 * seq(1:n) - 1)/(2 * n) )^2 )
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
                - n - mean( (2 * seq(1:n) - 1) * (log(theop) + log(1 - rev(theop))) ) 
            }
        else
        if (gof == "ADR")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                n/2 - 2 * sum(theop) - mean ( (2 * seq(1:n) - 1) * log(1 - rev(theop)) )
            }
        else
        if (gof == "ADL")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                -3*n/2 + 2 * sum(theop) - mean ( (2 * seq(1:n) - 1) * log(theop) )
            }
        else  
        if (gof == "AD2R")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                2 * sum(log(1 - theop)) + mean ( (2 * seq(1:n) - 1) / (1 - rev(theop)) )
            }
        else  
        if (gof == "AD2L")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                2 * sum(log(theop)) + mean ( (2 * seq(1:n) - 1) / theop )
            }
         else  
        if (gof == "AD2")
            fnobj <- function(par, fix.arg, obs, pdistnam)
            { 
                n <- length(obs)
                s <- sort(obs)
                theop <- do.call(pdistnam,c(list(q=s),as.list(par),as.list(fix.arg)))
                2 * sum(log(theop) + log(1 - theop) ) + 
                mean ( ((2 * seq(1:n) - 1) / theop) + ((2 * seq(1:n) - 1) / (1 - rev(theop))) )
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
            warnings("The function optim encountered an error and stopped")
            return(list(estimate = rep(NA,length(vstart)), convergence = 100, loglik = NA, 
                        hessian = NA))
        }
        
        if (opt$convergence>0) {
            warnings("The function optim failed to converge, with the error code ",
                     opt$convergence)
            return(list(estimate = rep(NA,length(vstart)), convergence = opt$convergence, 
                        value = NA, hessian = NA))
        }
        
        return(list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                    hessian = opt$hessian, gof=gof, optim.function="optim",
                    loglik=loglik(opt$par, fix.arg, data, ddistname) ))  
        
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
            print(opttryerror)
            warnings("The customized optimization function encountered an error and stopped")
            return(list(estimate = rep(NA,length(vstart)), convergence = 100, value = NA, 
                        hessian = NA))
        }
        
        if (opt$convergence>0) {
            warnings("The customized optimization function failed to converge, with the error code ",
                     opt$convergence)
            return(list(estimate = rep(NA,length(vstart)), convergence = opt$convergence, 
                        value = NA, hessian = NA))
        }
        
        return(list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                    gof=gof, hessian = opt$hessian, optim.function=custom.optim,
                    loglik=loglik(opt$par, fix.arg, data, ddistname) ))  

    }   
        
     
}

## old function with previous name for censored data
mledistcens<-function (censdata, distr, start=NULL,optim.method="default",lower=-Inf,upper=Inf)
{
    stop("The function \"mledistcens\" is no more used. Now the same function \"mledist\" must be used for censored and non censored data.")
}
