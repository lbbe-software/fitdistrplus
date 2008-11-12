mledist<-function (data, distr, start=NULL) 
{
    if (!is.character(distr)) distname<-substring(as.character(match.call()$distr),2)
    else distname<-distr
    ddistname<-paste("d",distname,sep="")
    if (!exists(ddistname,mode="function"))
        stop(paste("The ",ddistname," function must be defined"))
    if (distname == "unif")
        stop("Maximum likelihood estimation is not available for the uniform distribution")
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")
    # MLE fit 
    # definition of starting values if not previously defined
    if (is.null(start)) {
        if (distname == "norm") {
            n<-length(data)
            sd0 <- sqrt((n - 1)/n) * sd(data)
            mx <- mean(data)
            start <- list(mean=mx, sd=sd0)
        }
        if (distname == "lnorm") {
            if (any(data <= 0)) 
                stop("values must be positive to fit a lognormal distribution")
            n<-length(data)
            ldata <- log(data)
            sd0 <- sqrt((n - 1)/n) * sd(ldata)
            ml <- mean(ldata)
            start <- list(meanlog=ml, sdlog=sd0)
        }
        if (distname == "pois") {
            start<-list(lambda=mean(data))
        }
        if (distname == "exp") {
            start <- list(rate=1/mean(data))
        }
        if (distname == "gamma") {
            n<-length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            start<-list(shape=m^2/v,rate=m/v)
        }
        if (distname == "nbinom") {
            n<-length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            size <- if (v > m) m^2/(v - m)
                else 100
            start <- list(size = size, mu = m) 
        }
        if (distname == "geom" ) {
            m <- mean(data)
            prob<-if (m>0) 1/(1+m)
                    else 1
            start<-list(prob=prob)        
        }
        if (distname == "beta") {
            if (any(data < 0) | any(data > 1)) 
                stop("values must be in [0-1] to fit a beta distribution")
            n<-length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            aux<-m*(1-m)/v - 1
            start<-list(shape1=m*aux,shape2=(1-m)*aux)
        }
        if (distname == "weibull") {
            m <- mean(log(data))
            v <- var(log(data))
            shape <- 1.2/sqrt(v)
            scale <- exp(m + 0.572/shape)
            start <- list(shape = shape, scale = scale)
        }
        if (distname == "logis") {
            n<-length(data)
            m <- mean(data)
            v <- (n - 1)/n*var(data)
            start<-list(location=m,scale=sqrt(3*v)/pi)
        }
        if (distname == "cauchy") {
            start<-list(location=median(data),scale=IQR(data)/2)
        }
        if (!is.list(start)) 
            stop("'start' must be defined as a named list for this distribution") 
   } # end of the definition of starting values 
   
   # MLE fit using optim
    vstart<-unlist(start)
    if (length(vstart)>1) meth<-"Nelder-Mead"
    else meth<-"BFGS"
    fnobj<-function(par,dat,ddistnam) 
    -sum(log(do.call(ddistnam,c(list(x=dat),as.list(par)))))
    opt<-try(optim(par=vstart,fn=fnobj,dat=data,ddistnam=ddistname,hessian=TRUE,method=meth),silent=TRUE)
    if (inherits(opt,"try-error"))
    {
        warnings("The function optim encountered an error and stopped")
        return(list(estimate = rep(NA,length(vstart)), convergence = 100, loglik = NA, hessian = NA))
    }
    if (opt$convergence>0) {
        warnings("The function optim failed to converge, with the error code ",opt$convergence)
        return(list(estimate = rep(NA,length(vstart)), convergence = opt$convergence, loglik = NA, hessian = NA))
    }
    return(list(estimate = opt$par, convergence = opt$convergence, loglik = -opt$value, hessian = opt$hessian))       
}
