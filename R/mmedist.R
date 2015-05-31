#############################################################################
#   Copyright (c) 2010 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
		n <- length(data)
		m <- mean(data)
		v <- (n - 1)/n*var(data)
		
        if (!is.null(fix.arg))
            warnings("argument fix.arg cannot be used when a closed formula is used")
        # Fitting by matching moments
        if (!(is.vector(data) & is.numeric(data) & length(data)>1))
            stop("data must be a numeric vector of length greater than 1")
        if (distname == "norm") {
            estimate <- c(mean=m, sd=sqrt(v))
            order <- 1:2
        }
        if (distname == "lnorm") {
            if (any(data <= 0)) 
                stop("values must be positive to fit a lognormal distribution")
            sd2 <- log(1+v/m^2)
            estimate <- c(meanlog=log(m) - sd2/2, sdlog=sqrt(sd2))
            order <- 1:2            
        }
        if (distname == "pois") {
            estimate <- c(lambda=m)
            order <- 1          
        }
        if (distname == "exp") {
            estimate <- c(rate=1/m)
            order <- 1          
        }
        if (distname == "gamma" ) {
            shape <- m^2/v
            rate <- m/v
            estimate<-c(shape=shape, rate=rate)
            order <- 1:2            
       }
       if (distname == "nbinom" ) {
            size <- if (v > m) m^2/(v - m)
                    else NaN
            estimate<-c(size=size, mu=m)
            order <- 1:2           
       }
       if (distname == "geom" ) {
            prob<-if (m>0) 1/(1+m)
                    else NaN
            estimate<-c(prob=prob)
            order <- 1         
       }
        if (distname == "beta" ) {
            if (any(data < 0) | any(data > 1)) 
                stop("values must be in [0-1] to fit a beta distribution")
            aux<-m*(1-m)/v - 1
            shape1 <- m*aux
            shape2 <- (1-m)*aux
            estimate<-c(shape1=shape1, shape2=shape2)
            order <- 1:2            
       }
        if (distname == "unif" ) {
            min1 <- m-sqrt(3*v)
            max1 <- m+sqrt(3*v)
            estimate<-c(min1,max1)
            order <- 1:2            
       }
        if (distname == "logis" ) {
            scale <- sqrt(3*v)/pi
            estimate<-c(location=m, scale=scale)
            order <- 1:2            
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
                  opt.meth <- "Nelder-Mead"
               else 
                 opt.meth <- "BFGS"
            }else
              opt.meth <- "L-BFGS-B"
               
        }else
          opt.meth <- optim.method
        
        cens <- FALSE
        
        # Try to minimize the stat distance using the base R optim function
        if(is.null(custom.optim))
        {
            if (!cens)
                opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, mdistnam=mdistname, memp=memp,
                                            hessian=TRUE, method=opt.meth, lower=lower, upper=upper, ...), silent=TRUE)        
            else 
                stop("Moment matching estimation for censored data is not yet available.")
            
            if (inherits(opttryerror,"try-error"))
            {
                warnings("The function optim encountered an error and stopped.")
                if(getOption("show.error.messages")) print(attr(opttryerror, "condition"))              
                return(list(estimate = rep(NA,length(vstart)), convergence = 100, value = NA, 
                            hessian = NA))
            }
            
            if (opt$convergence>0) {
                warnings("The function optim failed to converge, with the error code ",
                         opt$convergence)
            }
            
            res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                        hessian = opt$hessian, optim.function="optim", order=order, memp=memp,
                        optim.method=opt.meth)  
            
        }else # Try to minimize the stat distance using a user-supplied optim function 
        {
            if (!cens)
                opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, mdistnam=mdistname, memp=memp,
                                                   par=vstart, ...), silent=TRUE)
            else
                stop("Moment matching estimation for censored data is not yet available.")
            
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
            
            res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                        hessian = opt$hessian, optim.function=custom.optim, order=order, memp=memp)  
            
        }   
        
    }
    
    loglik <- function(par, fix.arg, obs, ddistnam) {
        sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }
    if(exists(ddistname))
        loglik <- loglik(res$estimate, fix.arg, data, ddistname)
    else
        loglik <- NULL
    res <- c(res, fix.arg=fix.arg, loglik=loglik, method=meth)
    
    return(res)
}

## old function with previous name 
momdist<-function (data, distr) 
{
    stop("the name \"momdist\" for matching moments function is NO MORE used and is replaced by \"mmedist\".")
}
