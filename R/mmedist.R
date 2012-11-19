#############################################################################
#   Copyright (c) 2010 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis, Christophe Dutang                                                                                                  
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
### Matching moment estimation for non-censored data
###
###         R functions
### 

mmedist <- function (data, distr, order, memp, start=NULL, fix.arg=NULL,
    optim.method="default", lower=-Inf, upper=Inf, custom.optim=NULL, ...) 
{
    if (!is.character(distr)) 
#        distname <- substring(as.character(match.call()$distr), 2)
    stop("distr must be a character string naming a distribution")
    else 
        distname <- distr
    
    if (is.element(distname, c("norm", "lnorm", "pois", "exp", "gamma", "nbinom", "geom",
       "beta", "unif", "logis")))
        meth <- "closed formula"
    else
        meth <- optim.method
    
    mdistname <- paste("m", distname, sep="")
    ddistname <- paste("d", distname, sep="")

    if(meth != "closed formula")
    {
        if (!exists(mdistname, mode="function"))
            stop(paste("The moment function must be defined."))     
    # mdistname contains the good name of the theoretical moment function    
    }
	if (!(is.numeric(data) & length(data)>1)) 
		stop("data must be a numeric vector of length greater than 1.")
	
    
    if(meth == "closed formula")
    {
        if (!is.null(fix.arg))
            warnings("argument fix.arg cannot be used when a closed formula is used")
        # Fitting by matching moments
        if (!(is.vector(data) & is.numeric(data) & length(data)>1))
            stop("data must be a numeric vector of length greater than 1")
        if (distname == "norm") {
            n <- length(data)
            sd0 <- sqrt((n - 1)/n) * sd(data)
            mx <- mean(data)
            estimate <- c(mean=mx, sd=sd0)
            order <- 1:2
#           names(estimate) <- c("mean", "sd")   
        }
        if (distname == "lnorm") {
            if (any(data <= 0)) 
                stop("values must be positive to fit a lognormal distribution")
            n <- length(data)
#biased estimator by Jensen inequality
#            ldata <- log(data)
#            sd0 <- sqrt((n - 1)/n) * sd(ldata)
#            ml <- mean(ldata)
#
			sd2 <- log(1+var(data)/mean(data)^2)
			estimate <- c(meanlog=log(mean(data)) - sd2/2, sdlog=sqrt(sd2))
            order <- 1:2            
#           names(estimate) <- c("meanlog", "sdlog")
        }
        if (distname == "pois") {
            estimate <- c(lambda=mean(data))
            order <- 1          
#           names(estimate) <- "lambda" 
        }
        if (distname == "exp") {
            estimate <- c(rate=1/mean(data))
            order <- 1          
#           names(estimate) <- "rate" 
        }
        if (distname == "gamma" ) {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            shape <- m^2/v
            rate <- m/v
            estimate<-c(shape=shape, rate=rate)
            order <- 1:2            
#           names(estimate) <- c("shape", "rate")
       }
       if (distname == "nbinom" ) {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            size <- if (v > m) m^2/(v - m)
                    else NaN
            estimate<-c(size=size, mu=m)
            order <- 1:2           
#           names(estimate)<-c("size","mu")
       }
       if (distname == "geom" ) {
            m <- mean(data)
            prob<-if (m>0) 1/(1+m)
                    else NaN
            estimate<-c(prob=prob)
            order <- 1         
#           names(estimate)<-"prob"
       }
        if (distname == "beta" ) {
            if (any(data < 0) | any(data > 1)) 
                stop("values must be in [0-1] to fit a beta distribution")
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            aux<-m*(1-m)/v - 1
            shape1 <- m*aux
            shape2 <- (1-m)*aux
            estimate<-c(shape1=shape1, shape2=shape2)
            order <- 1:2            
#           names(estimate) <- c("shape1", "shape2")
       }
        if (distname == "unif" ) {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            min1 <- m-sqrt(3*v)
            max1 <- m+sqrt(3*v)
            estimate<-c(min1,max1)
            order <- 1:2            
#           names(estimate) <- c("min", "max")
       }
        if (distname == "logis" ) {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            scale <- sqrt(3*v)/pi
            estimate<-c(location=m, scale=scale)
            order <- 1:2            
#           names(estimate) <- c("location", "scale")
       }
        res <- list(estimate=estimate, convergence=0, order=order, memp=NULL)
    }else #an optimimisation has to be done
    {
        
        if(length(start) != length(order))
            stop("wrong dimension for the moment order to match.")
        if(!exists("memp", mode="function")) 
            stop("the empirical moment function must be defined.")

        
        ############# MME fit using optim or custom.optim ##########
        vstart <- unlist(start)
        vfix.arg <- unlist(fix.arg)
        # check of the names of the arguments of the density function
        argmdistname <- names(formals(mdistname))   
        m <- match(names(start), argmdistname)
        mfix <- match(names(vfix.arg), argmdistname)
        if (any(is.na(m)) || length(m) == 0)
			stop("'start' must specify names which are arguments to 'distr'")
        if (any(is.na(mfix)))
			stop("'fix.arg' must specify names which are arguments to 'distr'")
        # check that some parameters are not both in fix.arg and start
        minter <- match(names(start), names(fix.arg))
        if (any(!is.na(minter)))
			stop("a distribution parameter cannot be specified both in 'start' and 'fix.arg'")
        
        # definition of the function to minimize : least square
        #Cramer - von Mises
        DIFF2 <- function(par, fix.arg, order, obs, mdistnam, memp)
        {
            momtheo <- do.call(mdistnam, c(as.list(order), as.list(par), as.list(fix.arg)) )
            momemp <- as.numeric(memp(obs, order))
            (momemp - momtheo)^2
        }
        fnobj <- function(par, fix.arg, obs, mdistnam, memp)
            sum( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp)) )
        
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
        
        cens <- FALSE
        
        # Try to minimize the stat distance using the base R optim function
        if(is.null(custom.optim))
        {
            if (!cens)
                opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, mdistnam=mdistname, memp=memp,
                                            hessian=TRUE, method=meth, lower=lower, upper=upper, ...), silent=FALSE)        
            else 
                stop("Moment matching estimation for censored data is not yet available.")
            
            if (inherits(opttryerror,"try-error"))
            {
                warnings("The function optim encountered an error and stopped")
                print(opttryerror)				
                return(list(estimate = rep(NA,length(vstart)), convergence = 100, value = NA, 
                            hessian = NA))
            }
            
            if (opt$convergence>0) {
                warnings("The function optim failed to converge, with the error code ",
                         opt$convergence)
                return(list(estimate = rep(NA,length(vstart)), convergence = opt$convergence, 
                            value = NA, hessian = NA))
            }
            
            res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, hessian = opt$hessian, 
                        order=order, optim.function="optim", memp=memp)  
            
        }else # Try to minimize the stat distance using a user-supplied optim function 
        {
            if (!cens)
                opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, mdistnam=mdistname, memp=memp,
                                                   par=vstart, ...), silent=TRUE)
            else
                stop("Moment matching estimation for censored data is not yet available.")
            
            if (inherits(opttryerror,"try-error"))
            {
                warnings("The customized optimization function encountered an error and stopped")
                print(opttryerror)				
                return(list(estimate = rep(NA,length(vstart)), convergence = 100, value = NA, 
                            hessian = NA))
            }
            
            if (opt$convergence>0) {
                warnings("The customized optimization function failed to converge, with the error code ",
                         opt$convergence)
                return(list(estimate = rep(NA,length(vstart)), convergence = opt$convergence, 
                            value = NA, hessian = NA))
            }
            
            res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, hessian = opt$hessian, 
                        order=order, optim.function=custom.optim, memp=memp)  
            
        }   
        
    }
    
    loglik <- function(par, fix.arg, obs, ddistnam) {
        sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }
    if(exists(ddistname))
        loglik <- loglik(res$estimate, fix.arg, data, ddistname)
    else
        loglik <- NULL
	res <- c(res, fix.arg=fix.arg)
    
    return( c(res, list(loglik=loglik, method=meth)) )
    
 
}

## old function with previous name 
momdist<-function (data, distr) 
{
    stop("the name \"momdist\" for matching moments function is NO MORE used and is replaced by \"mmedist\".")
}
