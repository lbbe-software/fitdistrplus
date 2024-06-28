#############################################################################
#   Copyright (c) 2024 Christophe Dutang, Marie Laure Delignette-Muller
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
### Estimate J and A matrices for MME asymptotic covariance matrix proposed by
### I. Ibragimov and R. Has'minskii. Statistical Estimation - Asymptotic Theory. Springer-Verlag, 1981, p368
###
###         R functions
### 


#compute J^{-1} A J^{-T}, see below
mme.vcov <- function(par, fix.arg, order, obs, mdistnam, memp, weights, 
                     epsilon = sqrt(.Machine$double.eps), echo=FALSE)
{
  stopifnot(length(par) == length(order))
  if(!is.null(weights))
    stopifnot(length(weights) == length(obs))
  stopifnot(length(epsilon) == 1)
  stopifnot(epsilon > 0)
  
  Jmat <- mme.Jhat(par, fix.arg, order, mdistnam, epsilon)
  if(all(!is.na(Jmat)) && qr(Jmat)$rank == NCOL(Jmat))
  {  
    Jinv <- solve(Jmat)
    Amat <- mme.Ahat(par, fix.arg, order, obs, mdistnam, memp, weights)
    res <- Jinv %*% Amat %*% t(Jinv)
  }else 
    res <- NULL
  if(echo)
  {
    cat("Amat\n")
    print(Amat)
    cat("Jmat\n")
    print(Jmat)
    cat("Jinv\n")
    print(Jinv)
  }
  res
}


#compute E[X^(i+j)] - \overline{x^i}_n \overline{x^j}_n
mme.Ahat <- function(par, fix.arg, order, obs, mdistnam, memp, weights)
{
  stopifnot(length(par) == length(order))
  if(!is.null(weights))
    stopifnot(length(weights) == length(obs))
  
  mtheo <- function(order)
  {
    do.call(mdistnam, c(as.list(order), as.list(par), as.list(fix.arg)) )
  }
  if(is.null(weights))
  {
    memp2 <- function(order)
      as.numeric(memp(obs, order))
  }else
  {
    memp2 <- function(order)
      as.numeric(memp(obs, order, weights))
  }
  res <- matrix(nrow=length(order), ncol=length(order))
  for(i in 1:length(order))
  {
    for(j in i:length(order))
    {
      res[i, j] <- mtheo(i+j) - memp2(i)*memp2(j)
      if(i != j)
        res[j, i] <- res[i, j]
    }
  }
  res
}
  

#compute (\partial m_k(theta)/\partial theta_i)_{k,i}
mme.Jhat <- function(par, fix.arg, order, mdistnam, epsilon = sqrt(.Machine$double.eps))
{
  stopifnot(length(par) == length(order))
  stopifnot(length(epsilon) == 1)
  stopifnot(epsilon > 0)
  
  mtheo <- function(x, order, fix.arg)
  {
    do.call(mdistnam, c(as.list(order), as.list(x), as.list(fix.arg)) )
  }
  
  f <- function(i, k)
  {
    dfapprox(mtheo, i=i, x=par, order=order[k], fix.arg=fix.arg, epsilon=epsilon)
  }
  res <- sapply(1:length(order), function(i)
      sapply(1:length(order), function(k) f(i, k)))
  res
}

#compute \partial f/\partial x_i(x)
#see e.g. https://en.wikipedia.org/wiki/Numerical_differentiation
dfapprox <- function(f, i, x, ..., epsilon = sqrt(.Machine$double.eps))
{
  
  if(epsilon <= 0)
    epsilon <- 1e-6
  xiplus <- ximinus <- x
  xiplus[i] <- xiplus[i] + epsilon
  ximinus[i] <- ximinus[i] - epsilon
  
  fx <- try(f(x, ...))
  if(inherits(fx, "try-error"))
    stop("the function cannot be evaluated at x")
  fplus <- try(f(xiplus, ...))
  if(inherits(fplus, "try-error"))
    stop("the function cannot be evaluated at x+epsilon")
  fminus <- try(f(ximinus, ...))
  if(inherits(fminus, "try-error"))
    stop("the function cannot be evaluated at x-epsilon")
  
  #start with centered approximation
  res <- (fplus - fminus) / (2 * epsilon)
  if(any(is.infinite(res) | is.nan(res)))
  {
    idx <- is.infinite(res) | is.nan(res)
    #switch to uncentered
    res[idx] <- (fplus[idx] - fx[idx]) / epsilon
  }
  if(any(is.infinite(res) | is.nan(res)))
  {
    idx <- is.infinite(res) | is.nan(res)
    #switch to uncentered
    res[idx] <- (fx[idx] - fminus[idx]) / epsilon
  }
  res
}