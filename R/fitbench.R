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
### Benchmark of optimization algorithms to compute estimate
###
###         R functions
### 

fitbench <- function(data, distr, method, grad=NULL, control=list(trace=0, REPORT=1, maxit=1000), lower=-Inf, upper=+Inf, ..., npar=2) 
{
  if(method != "mle")
    stop("not supported")
  
  hasbound <- any(is.finite(lower) | is.finite(upper))
  hasgrad <- !is.null(grad)
  
  reslist <- NULL
  # 1 - methods without gradient without constraint (always possible)
  for(meth in c("BFGS", "Nelder", "CG")) #CG with FR update
  {
    res1fit$time <- system.time(res1fit <- mledist(data, distr=distr, optim.method=meth, control=control, ...))[3]
    reslist <- c(reslist, list(res1fit))
  }
  for(type in 2:3) #CG with PR or BS updates
  {
    res1fit$time <- system.time(res1fit <- mledist(data, distr=distr, optim.method="CG", control=c(control, type=type), ...))[3] 
    reslist <- c(reslist, list(res1fit)) 
  }
  fullname <- c("BFGS", "NM", paste0("CG", c("FR", "PR", "BS")))
  
  # 2 - methods without gradient with constraints
  if(hasbound)
  {
    for(meth in c("L-BFGS-B", "Nelder"))
    {
      res1fit$time <- system.time(res1fit <- mledist(data, distr=distr, optim.method=meth, control=control, lower=lower, upper=upper, ...))[3]
      reslist <- c(reslist, list(res1fit))
    }
    fullname <- c(fullname, "L-BFGS-B", "NM-B")
  }
  
  # 3 - methods with gradient without constraint
  if(hasgrad)
  {
    for(meth in c("BFGS", "CG")) #CG with FR update
    {
      res1fit$time <- system.time(res1fit <- mledist(data, distr=distr, gradient=grad, optim.method=meth, control=control, ...))[3]
      reslist <- c(reslist, list(res1fit))
    }
    for(type in 2:3) #CG with PR or BS updates
    {
      res1fit$time <- system.time(res1fit <- mledist(data, distr=distr, gradient=grad, optim.method="CG", control=c(control, type=type), ...))[3] 
      reslist <- c(reslist, list(res1fit)) 
    }
    fullname <- c(fullname, paste0("G-",c("BFGS", paste0("CG", c("FR", "PR", "BS")))) )
  }
  
  # 4 - methods with gradient with constraints
  if(hasbound && hasgrad)
  {
    for(meth in c("BFGS", "Nelder", "CG")) #CG with FR update
    {
      res1fit$time <- system.time(res1fit <- mledist(x, distr=distr, optim.method=meth, control=control, lower=lower, upper=upper, gradient=grad, ...))[3]
      reslist <- c(reslist, list(res1fit))
    }
    for(type in 2:3) #CG with PR or BS updates
    {
      res1fit$time <- system.time(res1fit <- mledist(x, distr=distr, optim.method="CG", control=c(control, type=type), lower=lower, upper=upper, gradient=grad, ...))[3] 
      reslist <- c(reslist, list(res1fit)) 
    }
    fullname <- c(fullname, paste0("G-", c("BFGS", "NM", paste0("CG", c("FR", "PR", "BS"))), "-B") )
  }
  names(reslist) <- fullname
  
  getval <- function(x)
    c(x$estimate, loglik=x$loglik, x$counts, x$time)
  
  #suspect behavior
  if(any(sapply(reslist, length) != 12))
  {
    print(sapply(reslist, length))
  }
  
  resmat <- sapply(reslist, getval)
  
#   resNM <- sapply(reslist[grep("NM", names(reslist))], getval)
#   resCGFR <- sapply(reslist[grep("CGFR", names(reslist))], getval)
#   resCGPR <- sapply(reslist[grep("CGPR", names(reslist))], getval)
#   resCGBS <- sapply(reslist[grep("CGBS", names(reslist))], getval)  
#   resBFGS <- sapply(reslist[grep("BFGS", names(reslist))], getval)  
  
  allname <- c(paste("fitted", names(reslist[[1]]$estimate)), "fitted loglik", "func. eval. nb.", "grad. eval. nb.", "time (sec)")
  
#  rownames(resNM) <- rownames(resCGFR) <- rownames(resCGPR) <- rownames(resCGBS) <- rownames(resBFGS) <- allname
  rownames(resmat) <- allname
  
  resmat
}