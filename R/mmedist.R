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
                     weights=NULL, silent=TRUE, gradient=NULL, checkstartfix=FALSE, 
                     calcvcov=FALSE, ...) 
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
  argddistname <- names(formals(ddistname))
  
  if(is.null(custom.optim))
    optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))
  
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
  
  if(meth != "closed formula")
  {
    if (!exists(mdistname, mode="function"))
      stop(paste0("The moment ", mdistname, " function must be defined."))     
    # mdistname contains the good name of the theoretical moment function    
  }
  if (!(is.numeric(data) & length(data)>1)) 
    stop("data must be a numeric vector of length greater than 1")
  checkUncensoredNAInfNan(data)
  
  if(is.null(weights))
  {  
    loglik <- function(par, fix.arg, obs, ddistnam) 
      sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
  }else
  {
    loglik <- function(par, fix.arg, obs, ddistnam) 
      sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
  }
  
  if(meth == "closed formula")
  {
    n <- length(data)
    if(is.null(weights))
    {
      m <- mean(data)
      v <- (n - 1)/n*var(data)
    }else #weighted version from util-wtdstat.R
    {
      m <- wtdmean(data, weights=weights)
      v <- wtdvar(data, weights=weights)
    }
    
    if (!is.null(fix.arg))
      stop("argument fix.arg cannot be used when a closed formula is used.")
    # Fitting by matching moments
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
      stop("data must be a numeric vector of length greater than 1")
    if (distname == "norm") 
    {
      estimate <- c(mean=m, sd=sqrt(v))
      mnorm <- function(order, mean, sd)
      {
        if(order == 1)
          return(mean)
        if(order == 2)
          return(sd^2+mean^2)
        if(order == 3)
          return(3*mean*sd^2+mean^3)
        if(order == 4)
          return(3*sd^4+6*mean^2*sd^2+mean^4)
      }
      order <- 1:2
    }
    if (distname == "lnorm") {
      if (any(data <= 0)) 
        stop("values must be positive to fit a lognormal distribution")
      sd2 <- log(1+v/m^2)
      estimate <- c(meanlog=log(m) - sd2/2, sdlog=sqrt(sd2))
      mlnorm <- function(order, meanlog, sdlog)
      {
        exp(order*meanlog + order^2*sdlog^2/2)
      }
      order <- 1:2            
    }
    if (distname == "pois") {
      estimate <- c(lambda=m)
      mpois <- function(order, lambda) 
      {
        ifelse(order == 1, lambda, lambda+lambda^2)
      }
      order <- 1          
    }
    if (distname == "exp") {
      estimate <- c(rate=1/m)
      mexp <- function(order, rate)
      {
        ifelse(order == 1, 1/rate, 2/rate^2)
      }
      order <- 1          
    }
    if (distname == "gamma" ) {
      shape <- m^2/v
      rate <- m/v
      estimate <- c(shape=shape, rate=rate)
      mgamma <- function(order, shape, rate)
      {
        res <- shape/rate
        if(order == 1)
          return(res)
        res <- res * (shape+1)/rate
        if(order == 2)
          return(res)
        res <- res * (shape+2)/rate
        if(order == 3)
          return(res)
        res <- res * (shape+3)/rate
        res
      }
      order <- 1:2            
    }
    if (distname == "nbinom" ) {
      size <- if (v > m) m^2/(v - m)
      else NaN
      estimate <- c(size=size, mu=m)
      mnbinom <- function(order, size, mu)
      {
        if(order == 1)
          return(mu)
        if(order == 2)
          return(mu+mu^2/size+mu^2)
        prob <- size/(size+mu)
        if(order == 3)
        {
          res <- size*prob^2 + 3*prob*(1-prob)*(size) + (1-prob)^2*(size+size^2)
          res <- res * (1-prob)/prob^3
        }else if(order == 4)
        {
          res <- size*prob^3 + 7*prob^2*(1-prob)*(size) + 6*prob*(1-prob)^2*(size+size^2) + (1-prob)^3*(2*size+3*size^2+size^3)
          res <- res * (1-prob)/prob^4
        }else
          stop("not implemented")
        res
      }
      order <- 1:2           
    }
    if (distname == "geom" ) {
      prob <- if (m>0) 1/(1+m)
      else NaN
      estimate <- c(prob=prob)
      mcumulgeom <- function(order, prob) #E[N(N-1)..(N-j+1)] = (1-p)^j j! / p^j
      {
        (1-prob)^order * factorial(order) / prob^order
      }
      mgeom <- function(order, prob)
      {
        if(order == 0)
          return(1)
        m1 <- mcumulgeom(1, prob)
        if(order == 1)
          return(m1)
        m2 <- mcumulgeom(2, prob) + m1
        if(order == 2)
          return(m2)
        m3 <- mcumulgeom(3, prob) + 3*m2 + 2*m1
        if(order == 3)
          return(m3)
        m4 <- mcumulgeom(4, prob) + 6*m3 + 11*m2 - 6*m1
        if(order == 4)
          return(m4)
        else
          stop("not yet implemented")
      }
      order <- 1         
    }
    if (distname == "beta" ) {
      if (any(data < 0) | any(data > 1)) 
        stop("values must be in [0-1] to fit a beta distribution")
      aux <- m*(1-m)/v - 1
      shape1 <- m*aux
      shape2 <- (1-m)*aux
      estimate <- c(shape1=shape1, shape2=shape2)
      mbeta <- function(order, shape1, shape2)
      {
        stot <- shape1+shape2
        ifelse(order == 1, shape1/stot,
               (shape1^3+2*shape1*shape2+shape1^2)/stot^2/(stot+1))
      }
      order <- 1:2            
    }
    if (distname == "unif" ) {
      min1 <- m-sqrt(3*v)
      max1 <- m+sqrt(3*v)
      estimate <- c(min1,max1)
      munif <- function(order, min, max)
      {
        ifelse(order == 1, (min+max)/2, (max^2 + min^2 + max*min)/3)
      }
      order <- 1:2            
    }
    if (distname == "logis" ) {
      scale <- sqrt(3*v)/pi
      estimate <- c(location=m, scale=scale)
      mlogis <- function(order, location, scale)
      {
        ifelse(order == 1, location, scale^2*pi^2/3 + location^2)
      }
      order <- 1:2            
    }
    if (exists(ddistname)) 
      loglikval <- loglik(estimate, fix.arg, data, ddistname)  
    else 
      loglikval <- NULL
    
    #create memp() for closed formula solution
    if(meth == "closed formula")
    {
      if(is.null(weights))
        memp <- function(x, order, weights) mean(x^order)
      else
        memp <- function(x, order, weights) sum(x^order * weights)/sum(weights)
    }
    #compute asymptotic covariance matrix proposed by
    #I. Ibragimov and R. Has'minskii. Statistical Estimation - Asymptotic Theory. Springer-Verlag, 1981, p368
    #see R/util-mmedist-vcov.R
    if(calcvcov)
    {
      #manually look for the appropriate theoretical raw moment function
      if (distname == "norm") 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mnorm, memp, weights)
      }else if (distname == "lnorm") 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mlnorm, memp, weights)          
      }else if (distname == "pois") 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mpois, memp, weights)                  
      }else if (distname == "exp") 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mexp, memp, weights)
      }else if (distname == "gamma" ) 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mgamma, memp, weights)
      }else if (distname == "nbinom" ) 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mnbinom, memp, weights)          
      }else if (distname == "geom" ) 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mgeom, memp, weights)    
      }else if (distname == "beta" ) 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mbeta, memp, weights)
      }else if (distname == "unif" ) 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, munif, memp, weights)      
      }else if (distname == "logis" ) 
      {
        varcovar <- mme.vcov(estimate, NULL, order, data, mlogis, memp, weights)                 
      }else
        varcovar <- NULL
    }else
      varcovar <- NULL
    #add names
    if(!is.null(varcovar))
      colnames(varcovar) <- rownames(varcovar) <- names(estimate)
    res <- list(estimate=estimate, convergence=0, value=NULL, hessian=NULL,
                optim.function=NULL, opt.meth=NULL, fix.arg=NULL, fix.arg.fun=NULL,
                weights=weights, counts=NULL, optim.message=NULL, 
                loglik= loglikval, method=meth, order=order, memp=NULL,
                vcov=varcovar)
  }else #an optimimisation has to be done, where fix.arg and start can be a function
  {
    if(is.vector(start)) #backward compatibility
      start <- as.list(start)
    
    if(!checkstartfix) #pre-check has not been done by fitdist() or bootdist()
    {
      #cat("checkstartfix is carried out\n")
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
        mean( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp)) )
    }else
    {
      DIFF2 <- function(par, fix.arg, order, obs, mdistnam, memp, weights)
      {
        momtheo <- do.call(mdistnam, c(as.list(order), as.list(par), as.list(fix.arg)) )
        momemp <- as.numeric(memp(obs, order, weights))
        (momemp - momtheo)^2
      }
      fnobj <- function(par, fix.arg, obs, mdistnam, memp, weights)
        mean( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp, weights)) )
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
      #compute asymptotic covariance matrix proposed by
      #I. Ibragimov and R. Has'minskii. Statistical Estimation - Asymptotic Theory. Springer-Verlag, 1981, p368
      #see R/util-mmedist-vcov.R
      if(calcvcov)
      {
        varcovar <- mme.vcov(opt$par, fix.arg, order, data, mdistname, memp, weights)
        #add names
        if(!is.null(varcovar))
          colnames(varcovar) <- rownames(varcovar) <- names(opt$par)
      }else
        varcovar <- NULL
      
      res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                  hessian = opt$hessian, optim.function=opt.fun, optim.method=opt.meth,
                  fix.arg=fix.arg, fix.arg.fun=fix.arg.fun, weights=weights, 
                  counts=opt$counts, optim.message=opt$message, 
                  loglik=ifelse(exists(ddistname), loglik(opt$par, fix.arg, data, ddistname), NULL),
                  method=meth, order=order, memp=memp, vcov=varcovar)  
      
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
      if(is.null(names(opt$par)))
        names(opt$par) <- names(vstart)
      argdot <- list(...)
      method.cust <- argdot$method
      #compute asymptotic covariance matrix proposed by
      #I. Ibragimov and R. Has'minskii. Statistical Estimation - Asymptotic Theory. Springer-Verlag, 1981, p368
      #see R/util-mmedist-vcov.R
      if(calcvcov)
      {
        varcovar <- mme.vcov(opt$par, fix.arg, order, data, mdistname, memp, weights)
        #add names
        if(!is.null(varcovar))
          colnames(varcovar) <- rownames(varcovar) <- names(opt$par)
      }else
        varcovar <- NULL
      
      res <- list(estimate = opt$par, convergence = opt$convergence, value = opt$value, 
                  hessian = opt$hessian, optim.function=custom.optim, optim.method=method.cust,
                  fix.arg=fix.arg, fix.arg.fun=fix.arg.fun, weights=weights, 
                  counts=opt$counts, optim.message=opt$message, 
                  loglik=ifelse(exists(ddistname), loglik(opt$par, fix.arg, data, ddistname), NULL),
                  method=meth, order=order, memp=memp, vcov=varcovar)  
    }   
  }
  return(res)
}

## old function with previous name 
momdist <- function (data, distr) 
{
  stop("the name \"momdist\" for matching moments function is deprecated and is replaced by \"mmedist\".")
}
