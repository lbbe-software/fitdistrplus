# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Main function :       npsurv.minimal                                    #
# ----------------------------------------------------------------------- #
# Original code from Yong Wang, 2020                                      #
# ----------------------------------------------------------------------- #

  
npsurv.minimal <- function(data, w=1, maxit=100, tol=1e-6, verb=0,
                           pkg="stats", ...) 
{
  pkg <- match.arg(pkg, c("limSolve", "stats"))
  
  x2 = icendata(data, w)
  # exact or right-censored only
  if(nrow(x2$o) == 0 || all(x2$o[,2] == Inf)) { 
    r0 = km(x2)
    r = list(f=r0$f, upper=max(x2$t, x2$o[,1]), convergence=TRUE, ll=r0$ll,
             maxgrad=0, numiter=1)
    return(r)
  }
  x = rbind(cbind(x2$t, x2$t), x2$o)
  
  nx = nrow(x)
  w = c(x2$wt, x2$wo)
  wr = sqrt(w)
  n = sum(w)
  upper = x2$upper
  dmat = Deltamatrix(x)
  
  left = dmat$left
  right = dmat$right
  D = dmat$Delta
  m = length(left)
  p = double(m)
  
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
    r = hcnm(w=w, D=D, p0=p, maxit=maxit, tol=tol, verb=verb)
    j = r$pf > 0
    f = idf(left[j], right[j], r$pf[j]) 
    r = list(f=f, upper=upper, convergence=r$convergence, method="hcnm", ll=r$ll,
             maxgrad=r$maxgrad, numiter=r$numiter)
    return(r)
  }
  
  # interval censored and small dataset
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
      cat("##### Iteration", i, "#####\n")
      cat("Log-likelihood: ", signif(ll, 6), "\n")
    }
    if(verb > 1) cat("Maximum gradient: ", signif(dmax, 6), "\n")
    if(verb > 2) {cat("Probability vector:\n"); print(p)} 
    j[which(j)-1 + aggregate(d, by=list(group=cumsum(j)), which.max)[,2]] = TRUE
    
    pj <- NNLS_constrSum(wr * S[,j,drop=FALSE], wr * 2, pkg=pkg, ...)
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
    if( ll <= ll.old + tol ) {converge=TRUE; break}
  }
  f = idf(left[j], right[j], p[j])
  r = list(f=f, upper=upper, convergence=converge, method="cnm", ll=ll,
           maxgrad=max(crossprod(w/P, D))-n, numiter=i)
  r
}
