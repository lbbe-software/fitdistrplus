# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Internal functions for handling interval censored data                  #
# ----------------------------------------------------------------------- #
# Original code from Yong Wang, 2020                                      #
# ----------------------------------------------------------------------- #



# Delta matrix = incidence matrix of maximal intersection intervals
# Delta[i,j] = 1 indicates that ith observation is within jth interval
# 
# see section 2.2, p3 of Wang & Taylor : Delta matrix is denoted S

# An interval is either (Li, Ri] if Li < Ri, or [Li, Ri] if Li = Ri. 
Deltamatrix = function(LR) {
  #remove NAs
  id_noNA <- rowSums(is.na(LR)) == 0
  LR <- LR[id_noNA,]
  
  L = LR[,1]
  R = LR[,2]
  ic = L != R             # inverval-censored
  nc = sum(ic)
  # tol = max(R[R!=Inf]) * 1e-8
  if(nc > 0) {
    L1 = L[ic] + max(R[R!=Inf]) * 1e-8       # open left endpoints
    LRc = cbind(c(L1, R[ic]), c(rep(0,nc), rep(1,nc)), rep(1:nc, 2))
    LRc.o = LRc[order(LRc[,1]),]
    j = which(diff(LRc.o[,2]) == 1)
    left = L[ic][LRc.o[j,3]]
    right = R[ic][LRc.o[j+1,3]]
  }
  else left = right = numeric(0)
  if(nrow(LR) - nc > 0) {
    ut = unique(L[!ic])
    jin = colSums(outer(ut, left, ">") & outer(ut, right, "<=")) > 0
    left = c(ut, left[!jin])     # remove those that contain exact obs.
    right = c(ut, right[!jin])
    o = order(left, right)
    left = left[o]
    right = right[o]
  }
  ## D = outer(L, left, "<=") & outer(R, right, ">=") 
  D = outer(L, left, "<=") & outer(R, right, ">=") &
    (outer(L, right, "<") | outer(R, left, "=="))  
  colnames(D) <- paste0("left=", round(left, 1), "-right=", round(right, 1))
  rownames(D) <- paste0("obs", 1:length(L))

  names(left) = names(right) = NULL
  list(left=left, right=right, Delta=D)
}

# interval distribution function, i.e., a distribution function defined on
# a set of intervals.

# left      Left endpoints of the intervals
# right     Right endpoints of the intervals
# p         Probability masses allocated to the intervals

idf = function(left, right, p) {
  if(length(left) != length(right)) stop("length(left) != length(right)")
  names(left) = names(right) = names(p) = NULL
  p = rep(p, length=length(left))
  f = list(left=left, right=right, p=p/sum(p))
  f
}

idf2data.frame <- function(object)
{
  if(inherits(object, "idf"))
     data.frame(left=object$left, right=object$right, p=object$p)
  else
    as.data.frame(object)
}

icendata = function(x, w=1) {
  if(is.null(x)) return(NULL)
  if(is.icendata(x)) {
    if(all(w == 1)) return(x)
    w = rep(w, length = length(x$t) + nrow(x$o))
    if(length(x$t) > 0) x$wt = x$wt * w[1:length(x$wt)]
    if(nrow(x$o) > 0) x$wo = x$wo * w[length(x$wt)+1:nrow(x$o)]
    return(x)
  }
  z = vector("list", 7)
  names(z) = c("t", "wt", "o", "wo", "i1", "upper", "u")
  if(is.vector(x)) x = cbind(x, x)
  if(!is.matrix(x)) x = as.matrix(x)
  if(ncol(x) == 3) {w = w * x[,3]; x = x[,1:2]}
  if(length(w) != nrow(x)) w = rep(w, len=nrow(x))
  iw = w > 0
  w = w[iw]
  x = x[iw,,drop=FALSE]
  o = order(x[,1], x[,2])
  x = x[o,]
  w = w[o]
  id = c(TRUE, diff(x[,1]) > 0 | diff(x[,2]) > 0)
  id[is.na(id)] = FALSE            # for Inf's
  w = aggregate(w, by=list(group=cumsum(id)), sum)[,2]
  x = x[id,]
  i = x[,1] == x[,2]
  z$t = x[i,1]
  names(z$t) = NULL
  z$wt = w[i]
  z$o = x[!i,1:2,drop=FALSE]
  dimnames(z$o) = list(NULL, c("L","R"))
  z$wo = w[!i]
  z$upper = max(x[,1])
  z$i1 = z$t != z$upper
  z$u = sort(unique(c(0, pmin(c(x[,1], x[,2]), z$upper))))
  class(z) = "icendata"
  z
}

is.icendata = function(x) "icendata" %in% class(x)

