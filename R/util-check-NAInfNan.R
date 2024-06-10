

checkUncensoredNAInfNan <- function(x)
{
  if(!is.numeric(x))
    stop("data is not a numeric vector.")
  if(any(is.nan(x)))
    stop("data contain NaN (not a numeric) values.")
  if(any(is.infinite(x)))
    stop("data contain Inf (infinite) values.")
  if(any(is.na(x)))
    stop("data contain NA values.")
  invisible(NULL)
}


checkCensoredDataFrameNAInfNan <- function(x)
{
  if(!is.data.frame(x))
    stop("censdata is not a dataframe with two columns named left and right.")
  if(NCOL(x) != 2)
    stop("censdata is not a dataframe with two columns named left and right.")
  if(!"left" %in% colnames(x) || !"right" %in% colnames(x))
    stop("censdata is not a dataframe with two columns named left and right.")
  if(any(!is.numeric(x$left) | !is.numeric(x$right)))
    stop("censdata contain NaN (not a numeric) values.")
  if(any(is.nan(x$left) | is.nan(x$right)))
    stop("censdata contain NaN (not a numeric) values.")
  if(any(is.infinite(x$left) | is.infinite(x$right)))
    stop("censdata contain Inf (infinite) values.")
  if(any(is.na(x$left) & is.na(x$right)))
    stop("censdata contain two NA values on the same line.")
  leftsupright <- x$left > x$right
  leftsupright <- leftsupright[!is.na(leftsupright)]
  if (any(leftsupright))
    stop("each censdata$left value must be less or equal to the corresponding censdata$right value")
  invisible(NULL)
}
