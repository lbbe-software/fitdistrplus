#############################################################################
#   Copyright (c) 2010 Christophe Dutang and Marie Laure Delignette-Muller
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
### quantile matching estimation for censored or non-censored data
###
###         R functions
### 

qmedist <- function (data, distr, probs, start=NULL, fix.arg=NULL, 
    qtype=7, optim.method="default", lower=-Inf, upper=Inf, custom.optim=NULL, 
    weights=NULL, silent=TRUE, gradient=NULL, ...)
    # data may correspond to a vector for non censored data or to
    # a dataframe of two columns named left and right for censored data 
{
    if (!is.character(distr)) 
        # distname <- substring(as.character(match.call()$distr), 2)
        stop("distr must be a character string naming a distribution")
    else 
        distname <- distr
    qdistname <- paste("q",distname,sep="")
    ddistname <- paste("d",distname,sep="")
    
    if (!exists(qdistname, mode="function"))
        stop(paste("The ", qdistname, " function must be defined"))
    if (!exists(ddistname, mode="function"))
        stop(paste("The ", ddistname, " function must be defined"))

    if (missing(probs))
        stop("missing probs argument for quantile matching estimation")
    if(is.null(custom.optim))
      optim.method <- match.arg(optim.method, c("default", "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN", "Brent"))

    start.arg <- start #to avoid confusion with the start() function of stats pkg (check is done lines 87-100)
    if(is.vector(start.arg)) #backward compatibility
      start.arg <- as.list(start.arg)
    
    if(qtype < 1 || qtype > 9)
        stop("wrong type for the R quantile function")
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
            stop("data must be a numeric vector of length greater than 1 for non censored data
            or a dataframe with two columns named left and right and more than one line for censored data")
    } else 
    {
        cens <- TRUE
        stop("Quantile matching estimation is not yet available for censored data.")
    }
    

    # QME fit 
    # definition of starting/fixed values values
    argddistname <- names(formals(ddistname))
    chfixstt <- checkparam(start.arg=start.arg, fix.arg=fix.arg, argdistname=argddistname, 
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

    if(length(vstart) != length(probs))
        stop("wrong dimension for the quantiles to match.")
    
   ############# QME fit using optim or custom.optim ##########

    # definition of the function to minimize : 
    # for non censored data
    if (!cens && is.null(weights)) 
    {
        # the argument names are:
        # - par for parameters (like in optim function)
        # - fix.arg for optional fixed parameters
        # - obs for observations (previously dat but conflicts with genoud data.type.int argument)
        # - qdistnam for distribution name
        
        DIFF2Q <- function(par, fix.arg, prob, obs, qdistnam, qtype)
        {
            qtheo <- do.call(qdistnam, c(as.list(prob), as.list(par), as.list(fix.arg)) )
            qemp <- as.numeric(quantile(obs, probs=prob, type=qtype))
            (qemp - qtheo)^2
        }
        
        fnobj <- function(par, fix.arg, obs, qdistnam, qtype)
          sum( sapply(probs, function(p) DIFF2Q(par, fix.arg, p, obs, qdistnam, qtype)) )
        
    }else if (!cens && !is.null(weights)) 
    {
      DIFF2Q <- function(par, fix.arg, prob, obs, qdistnam, qtype)
      {
        qtheo <- do.call(qdistnam, c(as.list(prob), as.list(par), as.list(fix.arg)) )
        qemp <- as.numeric(wtd.quantile(x=obs, weights=weights, probs=prob))
        (qemp - qtheo)^2
      }
      fnobj <- function(par, fix.arg, obs, qdistnam, qtype)
        sum( sapply(probs, function(p) DIFF2Q(par, fix.arg, p, obs, qdistnam, qtype)) )
    }
    
    # Function to calculate the loglikelihood to return
    if(is.null(weights))
    {  
      loglik <- function(par, fix.arg, obs, ddistnam) 
        sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }else
    {
      loglik <- function(par, fix.arg, obs, ddistnam) 
        sum(weights * log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
    }
    
    
    owarn <- getOption("warn")
    # Try to minimize the stat distance using the base R optim function
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
                    fix.arg=fix.arg, obs=data, qdistnam=qdistname, qtype=qtype, hessian=!is.null(gradient), 
                    method=meth, ...), silent=TRUE)
            
            if(!inherits(opttryerror, "try-error"))
              if(length(opt$counts) == 1) #appears when the initial point is a solution
                opt$counts <- c(opt$counts, NA)
            
            
          }else #opt.fun == "optim"
          {
            opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                  qdistnam=qdistname, qtype=qtype, hessian=TRUE, method=meth, lower=lower, upper=upper, ...), silent=TRUE)       
          }
          
        }else #hasbound == FALSE
        {
          opt.fun <- "optim"
          opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, gr=gradient,
                  qdistnam=qdistname, qtype=qtype, hessian=TRUE, method=meth, lower=lower, upper=upper, ...), silent=TRUE)       
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
					hessian = opt$hessian, probs=probs, optim.function=opt.fun, 
					loglik=loglik(opt$par, fix.arg, data, ddistname), fix.arg=fix.arg,
          optim.method=meth, weights = weights, counts=opt$counts, optim.message=opt$message)		
    }
    else # Try to minimize the stat distance using a user-supplied optim function 
    {
        if (!cens)
            opttryerror <- try(opt <- custom.optim(fn=fnobj, fix.arg=fix.arg, obs=data, qdistnam=qdistname, 
                qtype=qtype, par=vstart, ...), silent=TRUE)
        else
            stop("Quantile matching estimation is not yet available for censored data.")
        
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
					hessian = opt$hessian, probs=probs, optim.function=custom.optim, 
					loglik=loglik(opt$par, fix.arg, data, ddistname), fix.arg=fix.arg, 
          optim.method=NULL, weights = weights, counts=opt$counts, optim.message=opt$message)		
    }   
    return(res)    
     
}

