# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Internal optimi functions for Hierarchical constrained Newton Method    #
# ----------------------------------------------------------------------- #
# Original code from Yong Wang, 2020                                      #
# ----------------------------------------------------------------------- #





## ==========================================================================
## Hierarchical CNM: a variant of the Constrained Newton Method for finding
## the NPMLE survival function of a data set containing interval censoring.
## This is a new method to build on those in the Icens and MLEcens
## packages.  It uses the idea of block subsets of the S matrix to move
## probability mass among blocks of candidate support intervals.
##
## Usage (parameters and return value) is similar to the methods in package
## Icens, although note the transposed clique matrix.
##
## Arguments:
##   data: Data
##   w:  Weights
##   D: Clique matrix, n*m (note, transposed c.f. Icens::EMICM,
##      MLEcens::reduc).  The clique matrix may contain conditional
##      probabilities rather than just membership flags, for use in HCNM
##      recursively calling itself.
##   p0: Vector (length m) of initial estimates for the probabilities of
##      the support intervals.
##   maxit: Maximum number of iterations to perform
##   tol: Tolerance for the stopping condition (in log-likelihood value)
##   blockpar:
##     NA or NULL  means choose a value based on the data (using n and r)
##     ==0  means same as cnm (don't do blocks)
##      <1  means nblocks is this power of sj, e.g. 0.5 for sqrt
##      >1  means exactly this block size (e.g. 40)
##   recurs.maxit: For internal use only: maximum number of iterations in
##      recursive calls
##   depth: For internal use only: depth of recursion
##   verb: For internal use only: depth of recursion
##   pkg: package used in NNLS_constrSum()

## Author: Stephen S. Taylor and Yong Wang

## Reference: Wang, Y. and Taylor, S. M. (2013). Efficient computation of
## nonparametric survival functions via a hierarchical mixture
## formulation. Statistics and Computing, 23, 713-725.
## ==========================================================================

hcnm = function(data, w=1, D=NULL, p0=NULL, maxit=100, tol=1e-6,
                blockpar=NULL, recurs.maxit=2, depth=1, verb=0,
                pkg="stats", ...) {
  if(missing(D)) {
    x2 = icendata(data, w) #see npsurv-intercens.R
    if(nrow(x2$o) == 0 || all(x2$o[,2] == Inf)) { # exact or right-censored only
      r0 = km(x2)
      r = list(f=r0$f, convergence=TRUE, ll=r0$ll, maxgrad=0, numiter=1)
      class(r) = "npsurv"
      return(r)
    }
    x = rbind(cbind(x2$t, x2$t), x2$o)
    nx = nrow(x)
    w = c(x2$wt, x2$wo)
    dmat = Deltamatrix(x) #see npsurv-intercens.R
    left = dmat$left
    right = dmat$right
    intervals = cbind(left, right)
    D = dmat$Delta
  }
  else {
    if (missing(p0)) stop("Must provide 'p0' with D.")
    if (!missing(data)) warning("D and data both provided.  LR ignored!")
    nx = nrow(D)
    w = rep(w, length=nx)
    intervals = NULL
  }
  n = sum(w)
  wr = sqrt(w)
  converge = FALSE
  m = ncol(D)
  m1 = 1:m
  nblocks = 1
  maxdepth = depth
  i = rowSums(D) == 1
  r = mean(i)         # Proportion of exact observations
  if(is.null(p0)) {
    ## Derive an initial p vector.
    j = colSums(D[i,,drop=FALSE]) > 0
    while(any(c(FALSE,(i <- rowSums(D[,j,drop=FALSE])==0)))) {
      j[which.max(colSums(D[i,,drop=FALSE]))] = TRUE
    }
    p = colSums(w * D) * j
  }
  else { if(length(p <- p0) != m) stop("Argument 'p0' is the wrong length.") }
  p = p / sum(p)
  P = drop(D %*% p)
  ll = sum(w * log(P))
  evenstep = FALSE
  
  for(iter in 1:maxit) {
    p.old = p
    ll.old = ll
    S = D / P
    g = colSums(w * S)
    dmax = max(g) - n
    if(verb > 0) {
      cat("\n##### Iteration", iter, "#####\n")
      cat("Log-likelihood: ", signif(ll, 6), "\n")
    }
    if(verb > 1) cat("Maximum gradient: ", signif(dmax, 6), "\n")
    if(verb > 2) 
    {
      cat("Probability vector of length ", length(p), " :\n")
      if(length(p) > 12)
        cat(head(as.numeric(p)), "\t ...\t", tail(as.numeric(p)), "\n")
      else
        cat(as.numeric(p), "\n")
    } 
    j = p > 0
    if(depth==1) {
      s = unique(c(1,m1[j],m))
      if (length(s) > 1) for (l in 2:length(s)) {
        j[s[l-1] + which.max(g[s[l-1]:s[l]]) - 1] = TRUE
      }
    }
    sj = sum(j)
    ## BW: matrix of block weights: sj rows, nblocks columns
    if(is.null(blockpar) || is.na(blockpar))
      ## Default blockpar based on log(sj), Equation (14) p6 of Wang & Taylor
      iter.blockpar = ifelse(sj < 30, 0,
                             1 - log(max(20,10*log(sj/100)))/log(sj))
    else iter.blockpar = blockpar
    if(iter.blockpar==0 | sj < 30) {
      nblocks = 1
      BW = matrix(1, nrow=sj, ncol=1)
    }
    else {
      nblocks = max(1, if(iter.blockpar>1) round(sj/iter.blockpar)
                       else floor(min(sj/2, sj^iter.blockpar)))
      i = seq(0, nblocks, length=sj+1)[-1]
      if(evenstep) {
        nblocks = nblocks + 1
        BW = outer(round(i)+1, 1:nblocks, "==")
      }
      else BW = outer(ceiling(i), 1:nblocks, "==")
      storage.mode(BW) = "numeric"
    }

    for(block in 1:nblocks) {
      jj = logical(m)
      jj[j] = BW[,block] > 0
      sjj = sum(jj)
      if (sjj > 1 && (delta <- sum(p.old[jj])) > 0) {
        Sj = S[,jj]
        #original call
        #res = pnnls(wr * Sj, wr * drop(Sj %*% p.old[jj]) + wr, sum=delta) 
        #new call
        resNNLS <- NNLS_constrSum(a=wr * Sj, b=wr * drop(Sj %*% p.old[jj]) + wr, 
                              pkg=pkg, sumtotal=delta, control=list(trace=verb, maxit=5000), ...) #see npsurv-NNLS.R
        if(resNNLS$convergence != 0)
          warning("Problem NNLS_constrSum(a,b)")
        xj <- resNNLS$prob
        
        if(verb > 3) {cat("Block:", block, "\t delta", delta, "\t Optimized vector by NNLS:\n"); print(xj)} 
        
        p[jj] = p[jj] +  BW[jj[j],block] * (xj * (delta / sum(xj)) - p.old[jj])
      }
    }
    
    ## Maximise likelihood along the line between p and p.old
    p.gap = p - p.old              # vector from old to new estimate
    ## extrapolated rise in ll, based on gradient at old estimate
    ll.rise.gap = sum(g * p.gap) 
    alpha = 1
    p.alpha = p
    ll.rise.alpha = ll.rise.gap
    repeat {
      P = drop(D %*% p.alpha)
      ll = sum(w * log(P))
      if(ll >= ll.old && ll + ll.rise.alpha <= ll.old) {
        p = p.alpha               # flat land reached
        converge = TRUE
        break
      }
      if(ll > ll.old && ll >= ll.old + ll.rise.alpha * .33) {
        p = p.alpha               # Normal situation:  new ll is higher
        break
      }
      if((alpha <- alpha * 0.5) < 1e-10) {
        p = p.old
        P = drop(D %*% p)
        ll = ll.old
        converge = TRUE
        break
      }
      p.alpha = p.old + alpha * p.gap
      ll.rise.alpha = alpha * ll.rise.gap
    }
    if(converge) break

    if (nblocks > 1) {
      ## Now jiggle p around among the blocks
      Q = sweep(BW,1,p[j],"*")  # Matrix of weighted probabilities: [sj,nblocks]
      q = colSums(Q)            # its column sums (total in each block)
      ## Now Q is n*nblocks Matrix of probabilities for mixture components
      Q = sweep(D[,j] %*% Q, 2, q, "/")  
      if (any(q == 0)) {
        warning("A block has zero probability!")
      }
      else {
        if(verb >= 4)
          cat("recursion\n")
        ## Recursively call HCNM to allocate probability among the blocks 
        res = hcnm(w=w, D=Q, p0=q, blockpar=iter.blockpar,
                   maxit=recurs.maxit, recurs.maxit=recurs.maxit,
                   depth=depth+1, pkg=pkg, ...)
        maxdepth = max(maxdepth, res$maxdepth)
        if (res$ll > ll) {
          p[j] = p[j] * (BW %*% (res$pf / q))
          P = drop(D %*% p)
          ll = sum(w * log(P))  # should match res$lval
        }
      }
    }
    if(iter > 2) if( ll <= ll.old + tol ) {converge=TRUE; break}
    evenstep = !evenstep
  }
  list(pf=p, intervals=intervals, convergence=converge, method="hcnm", ll=ll,
       maxgrad=max(crossprod(w/P, D))-n, numiter=iter)
}
