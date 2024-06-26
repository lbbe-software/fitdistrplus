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
# testdpqfun function returns a vector of TRUE when the d,p,q functions exist
# and behave like in R (e.g. d,p,qexp()), otherwise a list of messages. 

# INPUTS 
#distr : the distribution name
#fun: a character vector with letters among d, p, q
#start.arg: an initial value list
#fix.arg: a fixed value list
#discrete: a logical whether the distribution is discrete

# OUTPUTS
# a vector of logical TRUE or a vector of text messages
testdpqfun <- function(distr, fun=c("d","p","q"), start.arg, 
                       fix.arg=NULL, discrete=FALSE)
{
  stopifnot(all(is.character(fun)))
  fun <- fun[fun %in% c("d","p","q")]
  stopifnot(length(fun) > 0)
  
  if(is.vector(start.arg)) 
    start.arg <- as.list(start.arg)
  if(is.function(fix.arg))
    stop("fix.arg should be either a named list or NULL but not a function")
  
  op <- options() #get current options
  #print(getOption("warn"))
  options(warn=-1)
  
  res <- NULL
  if("d" %in% fun)
    res <- rbind(res, test1fun(paste0("d", distr), start.arg, fix.arg))
  if("p" %in% fun)
    res <- rbind(res, test1fun(paste0("p", distr), start.arg, fix.arg))
  if("q" %in% fun)
    res <- rbind(res, test1fun(paste0("q", distr), start.arg, fix.arg))
  
  options(op)     # reset (all) initial options
  res
} 

test1fun <- function(fn, start.arg, fix.arg, dpqr)
{
  res <- data.frame(ok=FALSE, txt="")
  stopifnot(is.list(start.arg))
  if(!is.null(fix.arg))
    stopifnot(is.list(fix.arg))
  
  #does the function exist?
  if(!exists(fn, mode="function"))
  {
    res$txt <- paste("The", fn, "function must be defined")
    return(res)
  }
  
  #naming convention
  if(missing(dpqr))
    dpqr <- substr(fn, 1, 1)
  firstarg_theo <- switch(dpqr, "d"="x", "p"="q", "q"="p", "r"="n")
  firstarg_found <- names(formals(fn))[1]
  if(firstarg_found != firstarg_theo)
  {    
    t0 <- paste("The", fn, "function should have its first argument named:", firstarg_theo)
    res$txt <- paste(t0, "as in base R")
    return(res)
  }
  
  #zero-component vector
  res0 <- try(do.call(fn, c(list(numeric(0)), start.arg, fix.arg)), silent=TRUE)
  t0 <- paste("The", fn, "function should return a zero-length vector when input has length zero and not raise an error")
  t1 <- paste("The", fn, "function should return a zero-length vector when input has length zero")
  if(inherits(res0, "try-error"))
  {
    res$txt <- t0
    return(res)
  }
  if(length(res0) != 0)
  {
    res$txt <- t1
    return(res)
  }
  
  #inconsistent value
  x <- c(0, 1, Inf, NaN, -1)
  res1 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t2 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent values and not raise an error")
  if(inherits(res1, "try-error"))
  {
    res$txt <- t2
    return(res)
  }
  
  #missing value 
  x <- c(0, 1, NA)
  res2 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t4 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not raise an error")
  t5 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not remove missing values")
  if(inherits(res2, "try-error"))
  {
    res$txt <- t4
    return(res)
  }
  if(length(res2) != length(x))
  {
    res$txt <- t5
    return(res)
  }
  
  #inconsistent parameter
  x <- 0:1
  start.arg <- lapply(start.arg, function(x) -x)
  res3 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t6 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent parameters and not raise an error")
  if(inherits(res3, "try-error"))
  {
    res$txt <- t6
    return(res)
  }
  
  #wrong parameter name
  x <- 0:1
  names(start.arg) <- paste0(names(start.arg), "_")
  res4 <- try(do.call(fn, c(list(x), start.arg, fix.arg)), silent=TRUE)
  t8 <- paste("The", fn, "function should raise an error when names are incorrectly named")
  if(!inherits(res4, "try-error"))
  {
    res$txt <- t8
    return(res)
  }  
  return(data.frame(ok=TRUE, txt=""))
}
