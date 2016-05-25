#############################################################################
#   Copyright (c) 2016 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### prefit for maximum likelihood estimation
###
###         R functions
### 

#search good starting values
prefit <- function(data, distr, method = c("mle", "mme", "qme", "mge"), feasible.par, memp=NULL, order=NULL,
                   probs=NULL, qtype=7, gof=NULL, fix.arg=NULL, lower, upper, weights=NULL, silent=TRUE, ...)
{
  if (!is.character(distr)) 
    distname <- substring(as.character(match.call()$distr), 2)
  else 
    distname <- distr
  
  method <- match.arg(method, c("mle", "mme", "qme", "mge"))
  if(method != "qme" && !is.null(probs))
    stop("probs is not needed")
  if(method != "mme" && (!is.null(memp) || !is.null(order)))
    stop("memp, order are not needed")
  if(method != "mge" && !is.null(gof))
    stop("gof is not needed")
  
  
  ddistname <- paste0("d", distname)
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  pdistname <- paste0("p", distname)
  if (!exists(pdistname, mode="function") && method == "mge")
    stop(paste("The ", pdistname, " function must be defined"))
  qdistname <- paste0("q",distname)
  if (!exists(qdistname, mode="function") && method == "qme")
    stop(paste("The ", qdistname, " function must be defined"))
  mdistname <- paste0("m",distname)
  if (!exists(mdistname, mode="function") && method == "mme")
    stop(paste("The ", mdistname, " function must be defined"))
  
  if(is.null(probs) && method == "qme")
    stop("probs must be provided")
  if(is.null(order) && method == "mme")
    stop("order must be provided")
  
  if(missing(feasible.par))
    stop("feasible values must be provided")
  if(missing(lower) || missing(upper))
    stop("bounds (yet infinite) must be provided")
  
  #recycle parameters
  if(is.list(feasible.par))
    feasible.par <- unlist(feasible.par)
  npar <- length(feasible.par) #as in optim() line 34
  lower <- as.double(rep_len(lower, npar)) #as in optim() line 64
  upper <- as.double(rep_len(upper, npar))
  
  if(is.infinite(lower) && is.infinite(upper))
  {
    bnd <- detectbound(distname, feasible.par, data, fix.arg=fix.arg)
  }else
  {
    bnd <- rbind(lower, upper)
    colnames(bnd) <- names(feasible.par)
    rownames(bnd) <- c("lowb", "uppb")
  }
  if(!silent)
    print(bnd)
  
  translist <- invlist <- NULL
  for(i in 1:NCOL(bnd))
  {
    if(bnd["lowb", i] == -Inf && bnd["uppb", i] == Inf)
    {  
      translist <- c(translist, list(function(x) x))
      invlist <- c(invlist, list(function(x) x))
    }else if(bnd["lowb", i] == 0 && bnd["uppb", i] == Inf)
    {  
      translist <- c(translist, list(T0Inf))
      invlist <- c(invlist, list(iT0Inf))
    }else if(bnd["lowb", i] == 1 && bnd["uppb", i] == Inf)
    {  
      translist <- c(translist, list(T1Inf))
      invlist <- c(invlist, list(iT1Inf))
    }else if(bnd["lowb", i] == 0 && bnd["uppb", i] == 1)
    {  
      translist <- c(translist, list(T01))
      invlist <- c(invlist, list(iT01))
    }else if(bnd["lowb", i] == -1 && bnd["uppb", i] == 0)
    {  
      translist <- c(translist, list(Tm10))
      invlist <- c(invlist, list(iTm10))
    }else
    {
      print(bnd)
      stop("unknown parameter domain")
    }
  }
  if(!silent)
    print(translist)
  
  if(!is.null(weights))
  {
    if(any(weights < 0))
      stop("weights should be a vector of numerics greater than 0")
    if(length(weights) != NROW(data))
      stop("weights should be a vector with a length equal to the observation number")
    if(method == "mge")
      stop("weights is not allowed for maximum GOF estimation")
  }
  
  
  #maximum likelihood
  if(method == "mle")
  {
    if(is.null(weights))
      weights <- rep(1, NROW(data))
    fnobj <- function(par, fix.arg, obs, ddistnam, qdistnam, pdistnam, 
                      mdistnam, qtype, memp, gof) 
    {
      if(!is.list(par))
        par <- as.list(par)
      lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
      -sum( weights * log(do.call(ddistnam, c(list(obs), lpar, as.list(fix.arg)) ) ) )
    }
  }
  #quantile matching
  if(method == "qme" && is.null(weights))
  {
    DIFF2Q <- function(par, fix.arg, prob, obs, qdistnam, qtype)
    {
      if(!is.list(par))
        par <- as.list(par)
      lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
      qtheo <- do.call(qdistnam, c(list(prob), lpar, as.list(fix.arg)) )
      qemp <- as.numeric(quantile(obs, probs=prob, type=qtype))
      (qemp - qtheo)^2
    }
    
    fnobj <- function(par, fix.arg, obs, ddistnam, qdistnam, pdistnam, 
                      mdistnam, qtype, memp, gof) 
      sum( sapply(probs, function(p) DIFF2Q(par, fix.arg, p, obs, qdistnam, qtype)) )
  }
  if(method == "qme" && !is.null(weights))
  {
    DIFF2Q <- function(par, fix.arg, prob, obs, qdistnam, qtype)
    {
      if(!is.list(par))
        par <- as.list(par)
      lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
      qtheo <- do.call(qdistnam, c(list(prob), lpar, as.list(fix.arg)) )
      qemp <- as.numeric(wtd.quantile(x=obs, weights=weights, probs=prob))
      (qemp - qtheo)^2
    }
    fnobj <- function(par, fix.arg, obs, ddistnam, qdistnam, pdistnam, 
                      mdistnam, qtype, memp, gof) 
      sum( sapply(probs, function(p) DIFF2Q(par, fix.arg, p, obs, qdistnam, qtype)) )
  }  
  #moment matching
  if(method == "mme" && is.null(weights))
  {
    DIFF2 <- function(par, fix.arg, order, obs, mdistnam, memp, weights)
    {
      if(!is.list(par))
        par <- as.list(par)
      lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
      momtheo <- do.call(mdistnam, c(list(order), lpar, as.list(fix.arg)) )
      momemp <- as.numeric(memp(obs, order))
      (momemp - momtheo)^2
    }
    fnobj <- function(par, fix.arg, obs, ddistnam, qdistnam, pdistnam, 
                      mdistnam, qtype, memp, gof) 
      sum( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp)) )
  }
  if(method == "mme" && !is.null(weights))
  {
    DIFF2 <- function(par, fix.arg, order, obs, mdistnam, memp, weights)
    {
      if(!is.list(par))
        par <- as.list(par)
      lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
      momtheo <- do.call(mdistnam, c(list(order), lpar, as.list(fix.arg)) )
      momemp <- as.numeric(memp(obs, order, weights))
      (momemp - momtheo)^2
    }
    fnobj <- function(par, fix.arg, obs, ddistnam, qdistnam, pdistnam, 
                      mdistnam, qtype, memp, gof) 
      sum( sapply(order, function(o) DIFF2(par, fix.arg, o, obs, mdistnam, memp, weights)) )
  } 
  #gof matching
  if(method == "mge")
  {
    fnobj <- function(par, fix.arg, obs, ddistnam, qdistnam, pdistnam, 
                      mdistnam, qtype, memp, gof) 
    {
      if(!is.list(par))
        par <- as.list(par)
      lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
      n <- length(obs)
      s <- sort(obs)
      theop <- do.call(pdistnam, c(list(s), lpar, as.list(fix.arg)) )
      obspu <- seq(1,n)/n
      obspl <- seq(0,n-1)/n
      
      if (gof == "CvM")
        1/(12*n) + sum( ( theop - (2 * 1:n - 1)/(2 * n) )^2 )
      else if (gof == "KS")
        max(pmax(abs(theop-obspu),abs(theop-obspl)))
      else if (gof == "AD")
        - n - mean( (2 * 1:n - 1) * (log(theop) + log(1 - rev(theop))) ) 
      else if (gof == "ADR")
        n/2 - 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(1 - rev(theop)) )
      else if (gof == "ADL")
        -3*n/2 + 2 * sum(theop) - mean ( (2 * 1:n - 1) * log(theop) )
      else if (gof == "AD2R")
        2 * sum(log(1 - theop)) + mean ( (2 * 1:n - 1) / (1 - rev(theop)) )
      else if (gof == "AD2L")
        2 * sum(log(theop)) + mean ( (2 * 1:n - 1) / theop )
      else if (gof == "AD2")
        2*sum(log(theop) + log(1 - theop)) + mean(((2*1:n - 1) / theop) + ((2*1:n - 1) / (1 - rev(theop))))
          
    }
  }
  
  ltrans.par <- sapply(1:length(feasible.par), function(i) invlist[[i]](feasible.par[[i]]))
  if(!silent)
  {
    cat("before transform\n")
    print(unlist(feasible.par))
    cat("after transform\n")
    print(unlist(ltrans.par))
  }
  
  if(method == "mle")
    test1 <- try(fnobj(par=ltrans.par, fix.arg = fix.arg, obs=data, ddistnam = ddistname), silent=silent)
  if(method == "qme")
    test1 <- try(fnobj(par=ltrans.par, fix.arg = fix.arg, obs=data, qdistnam=qdistname, qtype=qtype), silent=silent)
  if(method == "mme")
    test1 <- try(fnobj(par=ltrans.par, fix.arg = fix.arg, obs=data, mdistnam=mdistname, memp=memp), silent=silent)
  if(method == "mge")
    test1 <- try(fnobj(par=ltrans.par, fix.arg = fix.arg, obs=data, pdistnam=pdistname, gof=gof), silent=silent)
  
  if(class(test1) == "try-error" || silent == FALSE)
    print(test1)
  
  #get old warning value and set it
  owarn <- options(warn=ifelse(silent, -1, 0))
  if(method == "mle")
    opttryerror <- try(opt <- optim(par=ltrans.par, fn=fnobj, fix.arg=fix.arg, obs=data, ddistnam=ddistname, 
                                  hessian=FALSE, method="BFGS", ...), silent=silent)       
  if(method == "qme")
    opttryerror <- try(opt <- optim(par=ltrans.par, fn=fnobj, fix.arg=fix.arg, obs=data, qdistnam=qdistname, 
                                    qtype=qtype, hessian=FALSE, method="BFGS", ...), silent=silent)
  if(method == "mme")
    opttryerror <- try(opt <- optim(par=ltrans.par, fn=fnobj, fix.arg=fix.arg, obs=data, mdistnam=mdistname, 
                                    memp=memp, hessian=FALSE, method="BFGS", ...), silent=silent) 
  if(method == "mge")
    opttryerror <- try(opt <- optim(par=ltrans.par, fn=fnobj, fix.arg=fix.arg, obs=data, pdistnam=pdistname, 
                                    gof=gof, hessian=FALSE, method="BFGS", ...), silent=silent) 
  
  #get back to old warning value
  on.exit(options(owarn), add=TRUE)
  
  
  if(class(opttryerror) == "try-error")
    stop("unsuccessful pre-fitting process")
  if(!silent)
    print(opt)
  
  if(opt$convergence %in% 0:1) #either successful or reached the iteration limit (see ?optim)
  {
    prefitpar <- unlist(sapply(1:length(opt$par), function(i) translist[[i]](opt$par[i])))
  }else
  {
    prefitpar <- rep(NA, length(opt$par))
  }
  names(prefitpar) <- names(feasible.par)
  
  as.list(prefitpar)
}


