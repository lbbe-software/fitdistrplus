# testdpqfun function returns a vector of TRUE when the d,p,q functions exist
# and behave like in R (e.g. d,p,qexp()), otherwise a list of messages. 

# INPUTS 
#distr : the distribution name
#fun: a character vector with letters among d, p, q
#start.arg: an initial value vector
#discrete: a logical whether the distribution is discrete

# OUTPUTS
# a vector of logical TRUE or a vector of text messages
testdpqfun <- function(distr, fun=c("d","p","q"), start.arg, 
                       fix.arg=NULL, discrete=FALSE)
{
  stopifnot(all(is.character(fun)))
  fun <- fun[fun %in% c("d","p","q")]
  stopifnot(length(fun) > 0)
  
  if(is.list(start.arg)) 
    start.arg <- as.vector(start.arg)
  if(is.function(fix.arg))
    stop("fix.arg should be either a named vector or NULL but not a function")
  
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

test1fun <- function(fn, start.arg, fix.arg)
{
  res <- data.frame(ok=FALSE, txt="")
  
  #does the function exist?
  if(!exists(fn, mode="function"))
  {
    res$txt <- paste("The", fn, "function must be defined")
    return(res)
  }
  
  #zero-component vector
  res0 <- try(do.call(fn, c(list(numeric(0)), as.list(start.arg), as.list(fix.arg))), silent=TRUE)
  t0 <- paste("The", fn, "function should return a zero-length vector when input has length zero and not raise an error")
  t1 <- paste("The", fn, "function should return a zero-length vector when input has length zero")
  if(class(res0) == "try-error")
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
  res1 <- try(do.call(fn, c(list(x), as.list(start.arg), as.list(fix.arg))), silent=TRUE)
  t2 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent values and not raise an error")
  if(class(res1) == "try-error")
  {
    res$txt <- t2
    return(res)
  }
  
  #missing value 
  x <- c(0, 1, NA)
  res2 <- try(do.call(fn, c(list(x), as.list(start.arg), as.list(fix.arg))), silent=TRUE)
  t4 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not raise an error")
  t5 <- paste("The", fn, "function should return a vector of with NA values when input has missing values and not remove missing values")
  if(class(res2) == "try-error")
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
  start.arg <- -1*start.arg
  res3 <- try(do.call(fn, c(list(x), as.list(start.arg), as.list(fix.arg))), silent=TRUE)
  t6 <- paste("The", fn, "function should return a vector of with NaN values when input has inconsistent parameters and not raise an error")
  if(class(res3) == "try-error")
  {
    res$txt <- t6
    return(res)
  }
  
  #wrong parameter name
  x <- 0:1
  names(start.arg) <- paste0(names(start.arg), "_")
  res4 <- try(do.call(fn, c(list(x), as.list(start.arg), as.list(fix.arg))), silent=TRUE)
  t8 <- paste("The", fn, "function should raise an error when names are incorrectly named")
  if(class(res4) != "try-error")
  {
    res$txt <- t8
    return(res)
  }  
  return(data.frame(ok=TRUE, txt=""))
}
