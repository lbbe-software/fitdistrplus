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
    optim.method="default", lower=-Inf, upper=Inf, custom.optim=NULL, 
    weights=NULL, silent=TRUE, gradient=NULL, ...) 
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
    if(is.null(custom.optim))
      optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
    
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

    if(meth != "closed formula")
    {
        if (!exists(mdistname, mode="function"))
            stop(paste0("The moment ", mdistname, " function must be defined."))     
    # mdistname contains the good name of the theoretical moment function    
    }
    if (!(is.numeric(data) & length(data)>1)) 
        stop("data must be a numeric vector of length greater than 1")
    
    if(meth == "closed formula")
    {
		  n <- length(data)
      if(is.null(weights))
      {
        m <- mean(data)
        v <- (n - 1)/n*var(data)
      }else #weighted version from util-wtdstat.R
      {
        m <- wtd.mean(data, weights=weights)
        v <- wtd.var(data, weights=weights)
      }
		  
        if (!is.null(fix.arg))
            stop("argument fix.arg cannot be used when a closed formula is used.")
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
        res <- list(estimate=estimate, convergence=0, value=NULL, hessian=NULL,
                    optim.function=NULL, order=order, memp=NULL, counts=NULL,
                    optim.message=NULL)
		    opt.meth <- fix.arg.fun <- NULL
		    
    }else #an optimimisation has to be done, where fix.arg and start can be a function
    {
        start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
        if(is.vector(start.arg)) #backward compatibility
          start.arg <- as.list(start.arg)
        
        # definition of starting/fixed values values
        argmdistname <- names(formals(mdistname))
        chfixstt <- checkparam(start.arg=start.arg, fix.arg=fix.arg, argdistname=argmdistname, 
                               errtxt=NULL, data10=head(data, 10), distname=distname)
        if(!chfixstt$ok)
          stop(chfixstt$txt)
        #unlist starting values as needed in optim()
        if(is.function(chfixstt$start.arg))
          vstart <- unlist(chfixstt$start.arg(data))
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
        
        if(length(vstart) != length(order))
            stop("wrong dimension for the moment order to match")
        if(missing(memp)) 
            stop("the empirical moment function must be defined")
        #backward compatibility when memp is the name of the function and not the function itself
        if(is.character(memp))
            memp <- get0(memp, envir=pos.to.env(1))
        
        
        #check the memp function
        if(!is.function(memp)) 
          stop("the empirical moment must be defined as a function")
        if(is.null(weights))
        {
          txt <- "the empirical moment function must be a two-argument function of 'x', 'order'"
          if(length(formals(memp)) != 2)
            stop(txt)
          if(any(names(formals(memp)) != c("x", "order")))
            stop(txt)
        }else
        {
          txt <- "the empirical moment function must be a three-argument function of 'x', 'order', 'weights'"
          if(length(formals(memp)) != 3)
            stop(txt)
          if(any(names(formals(memp)) != c("x", "order", "weights")))
            stop(txt)
        }

        
        ############# MME fit using optim or custom.optim ##########
        
        # definition of the function to minimize : least square (Cramer - von Mises type)
        if(is.null(weights))
        {
          DIFF2 <- function(par, fix.arg, order, obs, mdistnam, memp, weights)
          {
            momtheo <- do.call(mdistnam, c(as.list(order), as.list(par), as.list(fix.arg)) )
            momemp <- as.numeric(memp(obs, order))
            (momemp - momtheo)^2
          }
          fnobj <- function(par, fix.arg, obs, mdistnam, memp, weights)
            sum( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp)) )
        }else
        {
          DIFF2 <- function(par, fix.arg, order, obs, mdistnam, memp, weights)
          {
            momtheo <- do.call(mdistnam, c(as.list(order), as.list(par), as.list(fix.arg)) )
            momemp <- as.numeric(memp(obs, order, weights))
            (momemp - momtheo)^2
          }
          fnobj <- function(par, fix.arg, obs, mdistnam, memp, weights)
            sum( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp, weights)) )
        }
        
        cens <- FALSE
        if(cens)
          stop("Moment matching estimation for censored data is not yet available.")
        
        
        owarn <- getOption("warn")
        # Try to minimize the stat distance using the base R optim function
        if(is.null(custom.optim))
        {
            hasbound <- any(is.finite(lower) | is.finite(upper))
          
            # Choice of the optimization method  
            if (optim.method == "default")
            {
              opt.meth <- ifelse(length(vstart) > 1, "Nelder-Mead", "BFGS") 
            }else
              opt.meth <- optim.method
          
            if(opt.meth == "BFGS" && hasbound && is.null(gradient))
            {
              opt.meth <- "L-BFGS-B"
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
                if(opt.meth == "Nelder-Mead")
                  opt.fun <- "constrOptim"
                else if(opt.meth %in% c("L-BFGS-B", "Brent"))
                  opt.fun <- "optim"
                else
                {
                  txt1 <- paste("The method", opt.meth, "cannot be used by constrOptim() nor optim() without gradient and bounds.")
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
                                                      fix.arg=fix.arg, obs=data, mdistnam=mdistname, memp=memp, 
                                                      hessian=!is.null(gradient), method=opt.meth, weights=weights, 
                                                      ...), silent=TRUE)
                if(!inherits(opttryerror, "try-error"))
                  if(length(opt$counts) == 1) #appears when the initial point is a solution
                    opt$counts <- c(opt$counts, NA)
                
              }else #opt.fun == "optim"
              {
                opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                                mdistnam=mdistname, memp=memp, hessian=TRUE, method=opt.meth, 
                                                lower=lower, upper=upper, weights=weights, ...), silent=TRUE)       
              }
              
            }else #hasbound == FALSE
            {
              opt.fun <- "optim"
              opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                                              mdistnam=mdistname, memp=memp, hessian=TRUE, method=opt.meth, lower=lower, 
                                              upper=upper, ...), silent=TRUE)       
            }
            options(warn=owarn)
            
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
            if(is.null(names(opt$par)))
              names(opt$par) <- names(vstart)
            res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                        hessian = opt$hessian, optim.function=opt.fun, order=order, 
                        memp=memp, counts=opt$counts, optim.message=opt$message)  
            
        }else # Try to minimize the stat distance using a user-supplied optim function 
        {
            opt.meth <- NULL
            if (!cens)
            {
                options(warn=ifelse(silent, -1, 0))
                opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, mdistnam=mdistname, 
                                                  memp=memp, par=vstart, weights=weights, ...), silent=TRUE)
                options(warn=owarn)
            }else
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
                        hessian = opt$hessian, optim.function=custom.optim, order=order, 
                        memp=memp, counts=opt$counts, optim.message=opt$message)  
            
        }   
        
    }
    
    if(is.null(weights))
    {  
      loglik <- function(par, fix.arg, obs, ddistnam) 
        sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }else
    {
      loglik <- function(par, fix.arg, obs, ddistnam) 
        sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }
    if(exists(ddistname))
        loglik <- loglik(res$estimate, fix.arg, data, ddistname)
    else
        loglik <- NULL
    
    res <- c(res, fix.arg=fix.arg, loglik=loglik, method=meth, optim.method=opt.meth, 
             fix.arg.fun=fix.arg.fun, list(weights = weights))
    
    return(res)
}

## old function with previous name 
momdist<-function (data, distr) 
{
    stop("the name \"momdist\" for matching moments function is NO MORE used and is replaced by \"mmedist\".")
}
