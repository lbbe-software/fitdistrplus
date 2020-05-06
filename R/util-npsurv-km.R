# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Internal kaplan meier function                                          #
# ----------------------------------------------------------------------- #
# Original code from Yong Wang, 2020                                      #
# ----------------------------------------------------------------------- #


# Kaplan-Meier estimate of the survival function for right-censored data

km = function(data, w=1) {
  x = icendata(data, w)
  if(any(x$o[,2] != Inf))
    stop("Not all observations are exact or right-censored")
  if(nrow(x$o) == 0) {              # no right-censored observations
    f = idf(x$t, x$t, x$wt)
    ll = sum(x$wt * log(f$p))
    return(list(f=f, ll=ll))
  }
  c = colSums(x$wo * outer(x$o[,1], x$t, "<"))
  n = sum(x$wt, x$wo)                            # number of observations
  r = n - c - c(0,cumsum(x$wt))[1:length(x$t)]   # no. at risk
  S = cumprod(1 - x$wt/r)                        # survival prob.
  # tab = cbind(x$t, x$wt, c, r, S)
  p = rev(diff(rev(c(1,S,0))))
  dc = x$wt + c
  if(max(x$t) > max(x$o[,1])) {
    f = idf(x$t, x$t, p[-length(p)])
    ll = sum( x$wt * log(f$p) )
  }
  else {
    f = idf(c(x$t,max(x$o[,1])), c(x$t,Inf), p)
    ll = sum(c(x$wt, n - sum(x$wt)) * log(f$p))
  }
  list(f=f, ll=ll)
}

