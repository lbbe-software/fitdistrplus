# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Main function :       npsurvminimal                                    #
# ----------------------------------------------------------------------- #
# Original code from Yong Wang, 2020                                      #
# ----------------------------------------------------------------------- #

## Arguments:
##   data: Data
##   w:  Weights
##   maxit: Maximum number of iterations to perform
##   tol: Tolerance for the stopping condition (in log-likelihood value)
##   verb: For internal use only: depth of recursion (up to 4)
##   pkg: package used in NNLS_constrSum()
##   probtol: Tolerance for keeping output probabilities (default value has
##   been chosen so that the output is similar to original npsurv())
  
npsurvminimal <- function(data, w=1, maxit=100, tol=1e-6, verb=0,
                           pkg="stats", probtol=2e-4, ...) 
{
  #sanity checks
  pkg <- match.arg(pkg, c("stats"))
  if(length(maxit) > 1) stop("maxit should be a positive scalar")
  if(maxit <= 0) stop("maxit should be a positive scalar")
  if(length(tol) > 1) stop("maxit should be a positive scalar")
  if(tol <= 0) stop("maxit should be a positive scalar")
  if(length(verb) > 1) stop("maxit should be a non-negative scalar")
  if(verb < 0) stop("maxit should be a non-negative scalar")
  if(length(probtol) > 1) stop("maxit should be a positive probability")
  if(probtol <= 0 || probtol >= 1) stop("maxit should be a positive probability")
  if(sum(is.na(data))) stop("data should have no NA values")
  
  x2 = icendata(data, w) #see npsurv-intercens.R
  
  #sanity checks
  if(sum(apply(x2$o, 1, function(x) x[1] > x[2])) > 0)
    stop("some intervals have left bounds strictly greater than right bounds")
  
  # exact or right-censored only
  if(nrow(x2$o) == 0 || all(x2$o[,2] == Inf)) { 
    if(verb > 0)
      cat("call to km()\n")
    r0 = km(x2) #see npsurv-km.R
    r = list(f=r0$f, upper=max(x2$t, x2$o[,1]), convergence=TRUE, ll=r0$ll,
             maxgrad=0, numiter=1, method="km")
    return(r)
  }
  #get left/right interval bounds
  x = rbind(cbind(x2$t, x2$t), x2$o)
  
  nx = nrow(x)
  w = c(x2$wt, x2$wo)
  wr = sqrt(w)
  n = sum(w)
  upper = x2$upper
  # compute incidence matrix of maximal intersection intervals
  # see section 2.2, p3 of Wang & Taylor : Delta matrix is denoted S
  dmat = Deltamatrix(x) #see npsurv-intercens.R
  
  left = dmat$left
  right = dmat$right
  D = dmat$Delta
  m = length(left)
  p = double(m)
  if(verb > 0)
    cat("nb row orig data=", nx, "\nnb of maximal intersection intervals=", m, "\n")
  
  i = rowSums(D) != 1
  j = colSums(D[!i,,drop=FALSE]) > 0
  j[c(1,m)] = TRUE
  # Initial p must ensure P > 0
  repeat {                 
    jm = which.max(colSums(D[i,,drop=FALSE]))
    j[jm] = TRUE
    i[D[,jm]] = FALSE
    if( sum(i) == 0 ) break
  }
  p = colSums(w * D) * j
  p = p / sum(p)
  
  # interval censored and large dataset
  if(m >= 30) {                     ## Turn to HCNM
    if(verb > 0)
      cat("call to hcnm()\n")
    r = hcnm(w=w, D=D, p0=p, maxit=maxit, tol=tol, verb=verb) #see npsurv-hcnm.R
    
    j = r$pf > probtol
    
    f = idf(left=left[j], right=right[j], p=r$pf[j]) 
    #normalize prob vector
    f$p <- f$p/sum(f$p)
    r = list(f=f, upper=upper, convergence=r$convergence, method="hcnm", ll=r$ll,
             maxgrad=r$maxgrad, numiter=r$numiter, m=m)
    return(r)
  }
  
  # interval censored and small dataset
  if(verb > 0)
    cat("body of npsurvminimal()\n\n")
  P = drop(D %*% p)
  ll = sum( w * log(P) )
  converge = FALSE
  for(i in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / P
    d = colSums(w * S)
    dmax = max(d) - n
    if(verb > 0) {
      cat("\n##### Iteration", i, "#####\n")
      cat("Log-likelihood: ", signif(ll, 6), "\n")
    }
    if(verb > 1) cat("Maximum gradient: ", signif(dmax, 6), "\n")
    if(verb > 2) {cat("Probability vector:\n"); print(p)} 
    j[which(j)-1 + aggregate(d, by=list(group=cumsum(j)), which.max)[,2]] = TRUE
    #orignal call
    #pj = pnnls(wr * S[,j,drop=FALSE], wr * 2, sum=1)$x
    #new call
    resNNLS <- NNLS_constrSum(a=wr * S[,j,drop=FALSE], b=wr * 2, pkg=pkg, 
                              sumtotal=1, control=list(trace=verb), ...) #see npsurv-NNLS.R
    if(resNNLS$convergence != 0)
      break
    else
      pj <- resNNLS$prob
    p[j] = pj / sum(pj)
    # line search
    alpha = 1                
    pd = p - p.old
    lld = sum(d * pd)
    p.alpha = p
    repeat {
      P.alpha = drop(D %*% p.alpha)
      ll.alpha = sum(w * log(P.alpha))
      if(ll.alpha >= ll + alpha * lld * .33)
      { p = p.alpha; P = P.alpha; ll = ll.alpha; break }
      if((alpha <- alpha * .5) < 1e-10) break
      p.alpha = p.old + alpha * pd
    }
    j = p > 0
    if( ll <= ll.old + tol ) #convergence attained
    {
      converge=TRUE; break
    }
  }
  f <- idf(left[j], right[j], p[j]) #see npsurv-intercens.R
  res <- list(f=f, upper=upper, convergence=converge, method="cnm", ll=ll,
           maxgrad=max(crossprod(w/P, D))-n, numiter=i, m=m)
  res
}
