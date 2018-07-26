# checkparam function checks start.arg and fix.arg that parameters are named correctly

# INPUTS 
# start.arg : a named list
# fix.arg : NULL or a named list
# argdistname : argument names of the distribution
# hasnodefaultval : vector of logical indicating no default value of argument

# OUTPUTS 
# a named list with components: ok (TRUE or FALSE), txt (NULL or the error message), 
# start.arg : a named list of starting values for optimization 
# or a function to compute them from data
checkparamlist <- function(start.arg, fix.arg, argdistname, hasnodefaultval)
{
  errtxt <- list(t3="'start' must specify names which are arguments to 'distr'.",
          t4="'fix.arg' must specify names which are arguments to 'distr'.",
          t5="A distribution parameter cannot be specified both in 'start' and 'fix.arg'.",
          t6="'start' should not have NA or NaN values.",
          t7="'fix.arg' should not have NA or NaN values.",
          t8="Some parameter names have no starting/fixed value.",
          t9="Some parameter names have no starting/fixed value but have a default value.")
  
  vstart <- unlist(start.arg)
  #check unexpected names
  m <- match(names(vstart), argdistname)
  if (any(is.na(m))) 
    stop(errtxt$t3)
  #check NA/NaN values
  if(any(is.na(vstart) || is.nan(vstart)))
    stop(errtxt$t6)
  if(!is.null(fix.arg))
  {
    vfix <- unlist(fix.arg)
    #check unexpected names
    mfix <- match(names(vfix), argdistname)
    if (any(is.na(mfix))) 
      stop(errtxt$t4)
    
    # check that some parameters are not both in fix.arg and start
    minter <- match(names(vstart), names(vfix))
    if (any(!is.na(minter)))
      stop(errtxt$t5)
    
    #check NA/NaN values
    if(any(is.na(vfix) || is.nan(vfix)))
      stop(errtxt$t7)
    allparname <- names(c(vstart, vfix))
  }else
    allparname <- names(vstart)
  
  theoparam <- computegetparam(argdistname)
  #special case where both scale and rate are allowed, see ?dgamma
  if("scale" %in% theoparam && "rate" %in% theoparam)
  {
    errt8 <- any(!allparname %in% theoparam) || length(allparname) != length(theoparam)-1
  #special case where both prob and mu are allowed, see ?dnbinom
  }else if(length(theoparam) == 3 && all(c("size", "prob", "mu") %in% theoparam))
  {
    errt8 <- any(!allparname %in% theoparam) || length(allparname) != length(theoparam)-1
  }else
    errt8 <- any(!theoparam %in% allparname)
  #only make a warning if unset arguments have a default value
  if(errt8)
  {
    unsetarg <- theoparam[!theoparam %in% allparname] 
    if(any(hasnodefaultval[unsetarg]))
      stop(errtxt$t8)
    else
      warning(errtxt$t9)
  }
  
  list("start.arg"=start.arg, "fix.arg"=fix.arg)
}