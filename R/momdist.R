momdist<-function (data, distr) 
{
    if (!is.character(distr)) distname<-substring(as.character(match.call()$distr),2)
    else distname<-distr
    if (!is.element(distname,c("norm","lnorm","pois","exp","gamma","nbinom","geom",
        "beta","unif","logis")))
        stop(paste("Method of matching moments is not available for ",distname,
        " distribution"))
        # Ajouter d'autres distributions éventuellement
    # Fitting by matching moments
    if (!(is.vector(data) & is.numeric(data) & length(data)>1))
        stop("data must be a numeric vector of length greater than 1")
    if (distname == "norm") {
        n<-length(data)
        sd0 <- sqrt((n - 1)/n) * sd(data)
        mx <- mean(data)
        estimate <- c(mx, sd0)
        names(estimate) <- c("mean", "sd")
        return(estimate)        
    }
    if (distname == "lnorm") {
        if (any(data <= 0)) 
            stop("values must be positive to fit a lognormal distribution")
        n<-length(data)
        ldata <- log(data)
        sd0 <- sqrt((n - 1)/n) * sd(ldata)
        ml <- mean(ldata)
        estimate <- c(ml, sd0)
        names(estimate) <- c("meanlog", "sdlog")
        return(estimate)        
    }
    if (distname == "pois") {
        estimate <- mean(data)
        names(estimate) <- "lambda"
        return(estimate)        
    }
    if (distname == "exp") {
        estimate <- 1/mean(data)
        names(estimate) <- "rate"
        return(estimate)        
    }
    if (distname == "gamma" ) {
        n<-length(data)
        m <- mean(data)
        v <- (n - 1)/n*var(data)
        shape <- m^2/v
        rate <- m/v
        estimate<-c(shape,rate)
        names(estimate) <- c("shape", "rate")
        return(estimate)        
   }
   if (distname == "nbinom" ) {
        n<-length(data)
        m <- mean(data)
        v <- (n - 1)/n*var(data)
        size <- if (v > m) m^2/(v - m)
                else NaN
        estimate<-c(size,m)
        names(estimate)<-c("size","mu")
        return(estimate)        
   }
   if (distname == "geom" ) {
        m <- mean(data)
        prob<-if (m>0) 1/(1+m)
                else NaN
        estimate<-prob
        names(estimate)<-"prob"
        return(estimate)        
   }
    if (distname == "beta" ) {
        if (any(data < 0) | any(data > 1)) 
            stop("values must be in [0-1] to fit a beta distribution")
        n<-length(data)
        m <- mean(data)
        v <- (n - 1)/n*var(data)
        aux<-m*(1-m)/v - 1
        shape1 <- m*aux
        shape2 <- (1-m)*aux
        estimate<-c(shape1,shape2)
        names(estimate) <- c("shape1", "shape2")
        return(estimate)        
   }
    if (distname == "unif" ) {
        n<-length(data)
        m <- mean(data)
        v <- (n - 1)/n*var(data)
        min1 <- m-sqrt(3*v)
        max1 <- m+sqrt(3*v)
        estimate<-c(min1,max1)
        names(estimate) <- c("min", "max")
        return(estimate)        
   }
    if (distname == "logis" ) {
        n<-length(data)
        m <- mean(data)
        v <- (n - 1)/n*var(data)
        scale <- sqrt(3*v)/pi
        estimate<-c(m,scale)
        names(estimate) <- c("location", "scale")
        return(estimate)        
   }
 
}
