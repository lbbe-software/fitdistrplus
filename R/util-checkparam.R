# checkparam function checks start.arg and fix.arg that parameters are named correctly

# INPUTS 
# start.arg : starting values for optimization or the function to compute them from data
# fix.arg : fixed values of paramaters or the function to compute them from data
# argdistname : parameter names of the distribution
# errtxt : error text messages
# data10 : the first ten values of data
# distname : name of the distribution

# OUTPUTS 
# a named list with components: ok (TRUE or FALSE), txt (NULL or the error message), 
# start.arg : a named list of starting values for optimization 
# or a function to compute them from data
checkparam <- function(start.arg, fix.arg, argdistname, errtxt=NULL, data10, distname)
{
  if(is.null(errtxt))
    errtxt <- list(t0="Fixed values must be either a named list or a function returning a named list.",
          t1="Starting values must be either a named list or a function returning a named list.",
          t2="Starting and fixed values must be either a named list or a function returning a named list.",
          t3="'start' must specify names which are arguments to 'distr'.",
          t4="'fix.arg' must specify names which are arguments to 'distr'.",
          t5="A distribution parameter cannot be specified both in 'start' and 'fix.arg'.")
          #t6 = "Unknown starting values..."
  
  #before any treatment
  start.arg.was.null <- is.null(start.arg)
  
  #if clause with 4 different cases:
  #start.arg \ fix.arg | NULL | non NULL
  # NULL               | 1    | 2 
  # non NULL           | 3    | 4
  if(is.null(start.arg) && is.null(fix.arg)) #1
  { #default case from fitdist, mledist,...
    start.arg <- start.arg.default(data10, distr=distname)
  }else if(is.null(start.arg) && !is.null(fix.arg)) #2
  { #fix.arg should be a function or a named list
    if(!is.list(fix.arg) && !is.function(fix.arg))
      return(list(ok=FALSE, txt=errtxt$t0))
    
    #get param names
    if(is.function(fix.arg))
      namarg <- names(fix.arg(data10))
    else 
      namarg <- names(fix.arg)
    start.arg <- start.arg.default(data10, distr=distname) #could return "Unknown starting values..."
    start.arg <- start.arg[!names(start.arg) %in% namarg]
    
  }else if(!is.null(start.arg) && is.null(fix.arg)) #3
  {  #start should be a function or a named list
    if(!is.list(start.arg) && !is.function(start.arg))
      return(list(ok=FALSE, txt=errtxt$t1)) 
    
  }else if(!is.null(start.arg) && !is.null(fix.arg)) #4 
  {
    #fix.arg and start should be a function or a named list
    if( (!is.list(fix.arg) && !is.function(fix.arg)) || (!is.list(start.arg) && !is.function(start.arg)) )
      return(list(ok=FALSE, txt=errtxt$t2)) 
    
  }else
    stop("wrong implementation")
  
  #check start 
  #start.arg : function() | list()
  #start.arg cannot be null because set to a named list (by start.arg.default) when NULL 
  if(is.function(start.arg)) #a function
  {
    start2 <- start.arg(data10)
    if(!is.list(start2) && is.null(names(start2))) #check a named list
      return(list(ok=FALSE, txt=errtxt$t3))
    vstart <- unlist(start2)
  }else #a list
    vstart <- unlist(start.arg)
  m <- match(names(vstart), argdistname)
  if (any(is.na(m))) #check unexpected names
    return(list(ok=FALSE, txt=errtxt$t3))
  
  #check fix.arg
  #fix.arg : function() | list() | NULL
  if(is.function(fix.arg)) #a function
  {
    fix.arg2 <- fix.arg(data10)
    if(!is.list(fix.arg2) && is.null(names(fix.arg2))) #check a named list
      return(list(ok=FALSE, txt=errtxt$t4))
    vfix.arg <- unlist(fix.arg2)
  }else if(is.list(fix.arg)) #a list
    vfix.arg <- unlist(fix.arg)
  else
    vfix.arg <- NULL
  
  mfix <- match(names(vfix.arg), argdistname)
  if (any(is.na(mfix))) #check unexpected names
    return(list(ok=FALSE, txt=errtxt$t4))
  
  # check that some parameters are not both in fix.arg and start
  minter <- match(names(vstart), names(vfix.arg))
  if (any(!is.na(minter)))
    return(list(ok=FALSE, txt=errtxt$t5))
  
  #prepare the starg.arg for outputs, i.e. when start.arg=NULL, 
  # returns start.arg.default if not fixed param
  # returns a subset of start.arg.default if fixed param
  if(start.arg.was.null && is.null(fix.arg))
    start.arg <- function(x) start.arg.default(x, distr=distname) #could return "Unknown starting values..."
  else if(start.arg.was.null && !is.null(fix.arg))
  {
    if(is.function(fix.arg))
      namarg <- names(fix.arg(data10))
    else 
      namarg <- names(fix.arg)
    start.arg <- function(x){
      start.arg <- start.arg.default(x, distr=distname) #could return "Unknown starting values..."
      start.arg[!names(start.arg) %in% namarg]
    }
  }
  #otherwise start.arg is a named list or a function
  
  return(list(ok=TRUE, txt=NULL, start.arg=start.arg))
}