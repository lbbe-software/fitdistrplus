
#three functions from Hmisc also under GPL

#From wtd.stats.s (line 43) of the Hmisc package
wtd.quantile <- function(x, weights=NULL, probs=c(0, .25, .5, .75, 1), 
                         type='quantile', normwt=FALSE, na.rm=TRUE)
{
  if(!length(weights))
    return(quantile(x, probs=probs, na.rm=na.rm))
  
  type <- match.arg(type)
  if(any(probs < 0 | probs > 1))
    stop("Probabilities must be between 0 and 1 inclusive")
  
  nams <- paste(format(round(probs * 100, if(length(probs) > 1) 
    2 - log10(diff(range(probs))) else 2)), 
    "%", sep = "")
  
  w <- wtd.table(x, weights, na.rm=na.rm, normwt=normwt, type='list')
  x     <- w$x
  wts   <- w$sum.of.weights
  n     <- sum(wts)
  order <- 1 + (n - 1) * probs
  low   <- pmax(floor(order), 1)
  high  <- pmin(low + 1, n)
  order <- order %% 1
  ## Find low and high order statistics
  ## These are minimum values of x such that the cum. freqs >= c(low,high)
  allq <- approx(cumsum(wts), x, xout=c(low,high), 
                 method='constant', f=1, rule=2)$y
  k <- length(probs)
  quantiles <- (1 - order)*allq[1:k] + order*allq[-(1:k)]
  names(quantiles) <- nams
  return(quantiles)
  
}

#From wtd.stats.s (line 119) of the Hmisc package
wtd.table <- function(x, weights=NULL, type=c('list','table'), 
                      normwt=FALSE, na.rm=TRUE)
{
  type <- match.arg(type)
  if(!length(weights))
    weights <- rep(1, length(x))
  
  #isdate <- testDateTime(x)  ## 31aug02 + next 2
  ax <- attributes(x)
  ax$names <- NULL
  
  if(is.character(x)) x <- as.factor(x)
  lev <- levels(x)
  x <- unclass(x)
  
  if(na.rm) {
    s <- !is.na(x + weights)
    x <- x[s, drop=FALSE]    ## drop is for factor class
    weights <- weights[s]
  }
  
  n <- length(x)
  if(normwt)
    weights <- weights * length(x) / sum(weights)
  
  i <- order(x)  # R does not preserve levels here
  x <- x[i]; weights <- weights[i]
  
  if(anyDuplicated(x)) {  ## diff(x) == 0 faster but doesn't handle Inf
    weights <- tapply(weights, x, sum)
    if(length(lev)) {
      levused <- lev[sort(unique(x))]
      if((length(weights) > length(levused)) &&
           any(is.na(weights)))
        weights <- weights[!is.na(weights)]
      
      if(length(weights) != length(levused))
        stop('program logic error')
      
      names(weights) <- levused
    }
    
    if(!length(names(weights)))
      stop('program logic error')
    
    if(type=='table')
      return(weights)
    
    x <- all.is.numeric(names(weights), 'vector')
    #if(isdate)
    #  attributes(x) <- c(attributes(x),ax)
    
    names(weights) <- NULL
    return(list(x=x, sum.of.weights=weights))
  }
  
  xx <- x
  #if(isdate)
  #  attributes(xx) <- c(attributes(xx),ax)
  
  if(type=='list')
    list(x=if(length(lev))lev[x]
         else xx, 
         sum.of.weights=weights)
  else {
    names(weights) <- if(length(lev)) lev[x]
    else xx
    weights
  }
}

#From Misc.s (line 241) of the Hmisc package
all.is.numeric <- function(x, what=c('test','vector'),
                           extras=c('.','NA'))
{
  what <- match.arg(what)
  x <- sub('[[:space:]]+$', '', x)
  x <- sub('^[[:space:]]+', '', x)
  xs <- x[!x %in% c('',extras)] #originally %nin%
  isnum <- suppressWarnings(!any(is.na(as.numeric(xs))))
  if(what=='test')
    isnum
  else if(isnum)
    as.numeric(x)
  else x
}
