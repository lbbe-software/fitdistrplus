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
prefitmle <- function(data, distr, method = c("mle", "mme", "qme", "mge"), feasible.par, 
                      fix.arg=NULL, lower, upper, weights=NULL, silent=TRUE, ...)
{
  if (!is.character(distr)) 
    distname <- substring(as.character(match.call()$distr), 2)
  else 
    distname <- distr
  ddistname <- paste("d", distname, sep="")
  if (!exists(ddistname, mode="function"))
    stop(paste("The ", ddistname, " function must be defined"))
  
  pdistname <- paste("p", distname, sep="")
  if (!exists(pdistname, mode="function"))
    stop(paste("The ", pdistname, " function must be defined"))
  
  method <- match.arg(method, c("mle", "mme", "qme", "mge"))
  if(method != "mle")
    stop("not yet implemented")
  
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
  }else
    weights <- rep(1, NROW(data))
  
  
  #log likelihood (weighted)
  fnobj <- function(par, fix.arg, obs, ddistnam) 
  {
    if(!is.list(par))
      par <- as.list(par)
    #print(unlist(par))
    lpar <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
    #print(lpar)
    -sum( weights * log(do.call(ddistnam, c(list(obs), lpar, as.list(fix.arg)) ) ) )
  }
   
  ltrans.par <- sapply(1:length(feasible.par), function(i) invlist[[i]](feasible.par[[i]]))
  if(!silent)
  {
    cat("before transform\n")
    print(unlist(feasible.par))
    cat("after transform\n")
    print(unlist(ltrans.par))
  }
  test1 <- try(fnobj(par=ltrans.par, fix.arg = fix.arg, obs=data, ddistnam = ddistname), silent=silent)
  if(class(test1) == "try-error" || silent == FALSE)
    print(test1)
  
  #get old warning value and set it
  owarn <- options(warn=ifelse(silent, -1, 0))
  opttryerror <- try(opt <- optim(par=ltrans.par, fn=fnobj, fix.arg=fix.arg, obs=data, ddistnam=ddistname, 
                                  hessian=FALSE, method="BFGS", ...), silent=silent)       
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


