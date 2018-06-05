# compute pseudo data from a censored dataset

# INPUTS 
# censdata : a two-column matrix with left, right and/or interval censored data

# OUTPUTS 
# a vector of pseudo data (needed to compute starting values)
# a vector of right censored data
# a vector of left censored data
# a vector of non censored data
# a vector of interval censored data

cens2pseudo <- function(censdata)
{
  # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector, 
  #   icens (interval censored) = dataframe with left and right 
  # and ncens (not censored) = vector
  
  
  irow.lcens <- is.na(censdata$left) # rows corresponding to lcens data
  lcens <- censdata[irow.lcens, ]$right
  irow.rcens <- is.na(censdata$right)  # rows corresponding to rcens data
  rcens <- censdata[irow.rcens, ]$left
  if (any(is.na(lcens)) )
    stop("An observation cannot be both right and left censored, coded with two NA values")
  
  irow.ncens <- censdata$left==censdata$right & !is.na(censdata$left) & 
    !is.na(censdata$right)  # rows corresponding to ncens data
  ncens<-censdata[irow.ncens, ]$left
  irow.icens <- censdata$left!=censdata$right & !is.na(censdata$left) & 
    !is.na(censdata$right)  # rows corresponding to icens data
  icens<-censdata[irow.icens, ]
  
  pseudo <- c(rcens, lcens, ncens, (icens$left+icens$right)/2)
  list(pseudo=pseudo, rcens=rcens, lcens=lcens, ncens=ncens, icens=icens)
}


# compute row indexes from a censored dataset

cens2idxrow <- function(censdata)
{
  # Definition of datasets lcens (left censored)=vector, rcens (right censored)= vector, 
  #   icens (interval censored) = dataframe with left and right 
  # and ncens (not censored) = vector
  
  irow.lcens <- is.na(censdata$left) # rows corresponding to lcens data
  irow.rcens <- is.na(censdata$right)  # rows corresponding to rcens data
  irow.ncens <- censdata$left==censdata$right & !is.na(censdata$left) & 
    !is.na(censdata$right)  # rows corresponding to ncens data
  irow.icens <- censdata$left!=censdata$right & !is.na(censdata$left) & 
    !is.na(censdata$right)  # rows corresponding to icens data
  list(lcens=irow.lcens, rcens=irow.rcens, ncens=irow.ncens, icens=irow.icens)
}