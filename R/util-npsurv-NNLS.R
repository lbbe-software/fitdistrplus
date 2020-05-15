# ----------------------------------------------------------------------- #
# Nonparametric maximum likelihood estimation from interval-censored data #
# ----------------------------------------------------------------------- #
# Internal optimization functions for nonnegative least square method     #
# ----------------------------------------------------------------------- #

## ==========================================================================
## Nonnegative least squares (NNLS): 
##                                   
##        Minimize     ||ax - b||    
##        subject to non-negativity for all components  
##           x >= 0 i.e. diag(ncol(x)) x - 0 >= 0
##        sum up to one 
##         sum(x) = 1, ie. <x,1> - 1 >= 0, <x,-1> + 1 >= 0           
## ==========================================================================

NNLS_constrSum <- function(a, b, control=list(), pkg="stats", sumtotal=1, ...)
{
  pkg <- match.arg(pkg, c("stats"))
  
  controlnames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", 
                    "reltol", "alpha", "beta", "gamma", "REPORT", "warn.1d.NelderMead", "type", 
                    "lmm", "factr", "pgtol", "tmax", "temp")
  if(length(control) >= 1)
    control <- control[names(control) %in% controlnames]
  else
    control <- list(maxit=500)
  if (pkg == "stats")  #control parameter for BFGS used below
  {
    control$maxit <- max(control$maxit, 1000)
    control$reltol <- 1e-6
  }
  
  #sanity check
  if(!is.vector(b)) b = drop(b)
  if(!is.matrix(a)) stop("a not matrix")
  if(length(sumtotal) > 1) stop("sumtotal must be a positive scalar")
  if(!is.numeric(sumtotal)) stop("sumtotal must be a positive scalar")
  if(sumtotal <= 0) stop("sumtotal must be a positive scalar")
  m <- NROW(a)
  n <- NCOL(a)
  
  if (pkg == "stats")  
  {
    if(control$trace >= 2)
      cat("nb of optimized variables", n, "\n")
    if(n > 2)
    {
      
      #residual least square sum with theta=x[1:(n-1)] ; x[n] = 1-sum(theta)
      RLS <- function(theta) 
      {
        x <- c(theta, sumtotal-sum(theta))
        y <- a %*% x - b
        sum(y^2)
      }
      gradRLS <- function(theta)
        {
          x <- c(theta, sumtotal-sum(theta))
          diffa <- a[, 1:(NCOL(a)-1)] - a[, NCOL(a)]
          y <- a %*% x - b
          2*crossprod(diffa, y)
        }
      
      # non negativity constraint
      one_n <- rep(1, n-1)
      ui <- diag(n-1)
      rownames(ui) <- c(paste0("theta", 1:(n-1)))
      ci <- rep(0, n-1)
      
      #initial guess
      x0 <- rep(1/(n-1), n-1)
      #call to constrOptim
      control$trace <- control$trace >=5
      res <- constrOptim(theta=x0, f=RLS, grad=gradRLS, ui=ui, 
                         ci=ci, method="BFGS", control=control, ...)
      #warning : convergence is not reached if(res$convergence != 0)
      res$prob <- c(res$par, sumtotal-sum(res$par))
      #sanity check
      res$prob <- pmax(res$prob, 0)
      res$prob <- res$prob/sum(res$prob)*sumtotal
    }else if(n == 2)
    {
      #residual least square sum with theta=x[1:(n-1)] ; x[n] = 1-sum(theta)
      RLS <- function(theta) 
      {
        x <- c(theta, sumtotal-sum(theta))
        y <- a %*% x - b
        sum(y^2)
      }
      gradRLS <- function(theta)
      {
        x <- c(theta, sumtotal-sum(theta))
        diffa <- a[,1]-a[,2]
        2*crossprod(diffa, a %*% x - b)
      }
      
      # non negativity constraint
      one_n <- rep(1, n-1)
      ui <- diag(n-1)
      rownames(ui) <- c(paste0("theta", 1:(n-1)))
      ci <- rep(0, n-1)
      
      #initial guess
      x0 <- rep(1/(n-1), n-1)
      #call to constrOptim
      
      res <- constrOptim(theta=x0, f=RLS, grad=gradRLS, ui=ui, 
                         ci=ci, method="BFGS", control=control, ...)
      #warning : convergence is not reached if(res$convergence != 0)
      res$prob <- c(res$par, sumtotal-sum(res$par))
      #sanity check
      res$prob <- pmax(res$prob, 0)
      res$prob <- res$prob/sum(res$prob)*sumtotal
    }else if(n == 1)
    {
      xstar <- sumtotal
      res <- list(prob=xstar, value=(a*xstar-b)^2, counts=0, par=NA,
                  convergence=1*(xstar < 0), message="no optimization, because fully constrained problem")
    }else
      stop("wrong argument a for NNLS_constrSum()")
  }else 
    #   if (pkg == "limSolve")  
    # {
    #   #TODO : pass control argument to limSolve::lsei
    #   require(limSolve)
    #   res <- limSolve::lsei(A=a, B=b, E=rep(1,n), F=1, G=diag(n), H=rep(0, n), ...)
    #   if(res$IsError)
    #     stop("error in limSolve::lsei when computing NNLS")
    #   xstar <- res$X
    # }else
    stop("wrong package")
  
  res
}

