#start : starting values for optimization
#fix.arg : fixed values of paramaters
#argdistname : parameter names of the distribution
#errtxt : error text messages
#data10 : the first ten values of data
#distname : name of the distribution
checkparam <- function(start.arg, fix.arg, argdistname, errtxt=NULL, data10, distname)
{
  if(is.null(errtxt))
    errtxt <- list(t0="Fixed values must be either a named list or a function returning a named list.",
          t1="Starting values must be either a named list or a function returning a named list.",
          t2="Starting and fixed values must be either a named list or a function returning a named list.",
          t3="'start' must specify names which are arguments to 'distr'.",
          t4="'fix.arg' must specify names which are arguments to 'distr'.",
          t5="A distribution parameter cannot be specified both in 'start' and 'fix.arg'.",
          t6="Starting values must be defined when some distribution parameters are fixed.",
          t7="Missing starting values.")
  
  start.arg.class <- class(start.arg)
  
  if(is.null(start.arg) && is.null(fix.arg))
  {
    start.arg <- start.arg.default(data10, distr=distname)
  }else if(is.null(start.arg) && !is.null(fix.arg))
  {
    if(!is.list(fix.arg) && !is.function(fix.arg))
      return(list(ok=FALSE, txt=errtxt$t0))
    
    if(is.function(fix.arg))
      namarg <- names(fix.arg(data10))
    else 
      namarg <- names(fix.arg)
    start.arg <- start.arg.default(data10, distr=distname)
    start.arg <- start.arg[!names(start.arg) %in% namarg]
    
  }else if(!is.null(start.arg) && is.null(fix.arg))
  {
    #start should be a function or a named list
    if(!is.list(start.arg) && !is.function(start.arg))
      return(list(ok=FALSE, txt=errtxt$t1)) 
    
  }else if(!is.null(start.arg) && !is.null(fix.arg))
  {
    #fix.arg and start should be a function or a named list
    if(!is.list(fix.arg) && !is.function(fix.arg) && !is.list(start.arg) && !is.function(start.arg))
      return(list(ok=FALSE, txt=errtxt$t2)) 
    
  }else
    stop("wrong implementation")
  
  #check start 
  if(is.function(start.arg))
  {
    start2 <- start.arg(data10)
    if(!is.list(start2) && is.null(names(start2)))
      stop(errtxt$t3)
    vstart <- unlist(start2)
  }else #a list
    vstart <- unlist(start.arg)
  m <- match(names(vstart), argdistname)
  if (any(is.na(m)) || length(m) == 0)
    return(list(ok=FALSE, txt=errtxt$t3))
  
  #check fix.arg when non-NULL
  if(is.function(fix.arg))
  {
    fix.arg2 <- fix.arg(data10)
    if(!is.list(fix.arg2) && is.null(names(fix.arg2)))
      stop(errtxt$t4)
    vfix.arg <- unlist(fix.arg2)
  }else if(is.list(fix.arg)) #a list
    vfix.arg <- unlist(fix.arg)
  else
    vfix.arg <- NULL
  
  mfix <- match(names(vfix.arg), argdistname)
  if (any(is.na(mfix)))
    return(list(ok=FALSE, txt=errtxt$t4))
  
  # check that some parameters are not both in fix.arg and start
  minter <- match(names(vstart), names(vfix.arg))
  if (any(!is.na(minter)))
    return(list(ok=FALSE, txt=errtxt$t5))
  
  if(start.arg.class == "NULL" && is.null(fix.arg))
    start.arg <- function(x) start.arg.default(x, distr=distname)
  else if(start.arg.class == "NULL" && !is.null(fix.arg))
    start.arg <- function(x){
      namarg <- names(fix.arg)
      start.arg <- start.arg.default(x, distr=distname)
      start.arg[!names(start.arg) %in% namarg]
    } 
  #otherwise start.arg is a named list or a function
  
  return(list(ok=TRUE, txt=NULL, start.arg=start.arg))
}