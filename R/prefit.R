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
prefitmle <- function(data, ddistnam, vstart, fix.arg=NULL, lower, upper, custom.optim=NULL, weights=NULL, silent=TRUE, ...)
{
  if(missing(vstart))
    stop("starting values must be provided")
  if(missing(lower) || missing(upper))
    stop("bounds (yet infinite) must be provided")
  
  #recycle parameters
  npar <- length(vstart) #as in optim() line 34
  lower <- as.double(rep_len(lower, npar)) #as in optim() line 64
  upper <- as.double(rep_len(upper, npar))
  
  if(is.infinite(lower) && is.infinite(upper))
  {
    bnd <- detectbound(substr(ddistnam, 2, nchar(ddistnam)), vstart, data, fix.arg=fix.arg)
  }else
  {
    bnd <- rbind(lower, upper)
    colnames(bnd) <- names(vstart)
    rownames(bnd) <- c("lowb", "uppb")
  }
  
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
      stop("unknown parameter domain")
  }
  
  fnobj <- function(par, fix.arg, obs, ddistnam) 
  {
    if(!is.list(par))
      par <- as.list(par)
    #print(par)
    par <- lapply(1:length(par), function(i) translist[[i]](par[[i]]))
    #print(par)
    -sum(log(do.call(ddistnam, c(list(obs), as.list(par), as.list(fix.arg)) ) ) )
  }
  
  test1 <- try(fnobj(par=vstart, fix.arg = fix.arg, obs=data, ddistnam = ddistnam), silent=silent)
  if(class(test1) == "try-error")
    print(test1)
  
  
  
  opttryerror <- try(opt <- optim(par=vstart, fn=fnobj, fix.arg=fix.arg, obs=data, ddistnam=ddistnam, 
                                  hessian=FALSE, method="BFGS", ...), silent=silent)       
  
  if(opt$convergence %in% 0:1) #either successful or reached the iteration limit (see ?optim)
  {
    prefitpar <- unlist(sapply(1:length(opt$par), function(i) translist[[i]](opt$par[i])))
  }else
  {
    prefitpar <- rep(NA, length(par))
  }
  names(prefitpar) <- names(par)
  
  if(!silent)
    print(opt)
  prefitpar
}


