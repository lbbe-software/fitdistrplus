bootdistcens<-function (f, niter=999)
{ 
    if (niter<10) 
        stop("niter must be an integer above 10")
    if (!inherits(f, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    # non parametric bootstrap
    n<-length(f$censdata[,1])
    numrow<-seq(1,n)
    rnumrow<-sample(numrow,size=niter*n,replace=TRUE)
    dim(rnumrow)<-c(n,niter)
    start<-f$estimate
    funcmle<-function(iter) {
    mle<-mledistcens(data.frame(left=f$censdata[rnumrow[,iter],]$left,
            right=f$censdata[rnumrow[,iter],]$right),f$distname,start)
        return(c(mle$estimate,mle$convergence))
    }
    resboot<-sapply(1:niter,funcmle)
    rownames(resboot)<-c(names(start),"convergence")
    if (length(resboot[,1])>2) {
        estim<-data.frame(t(resboot)[,-length(resboot[,1])])
        bootCI <- cbind(apply(resboot[-length(resboot[,1]),],1,median,na.rm=TRUE),
            apply(resboot[-length(resboot[,1]),],1,quantile,0.025,na.rm=TRUE),
            apply(resboot[-length(resboot[,1]),],1,quantile,0.975,na.rm=TRUE))
        colnames(bootCI) <- c("Median","2.5%","97.5%")
    }
    else {
        estim<-as.data.frame(resboot[,-length(resboot[,1])])
        names(estim)<-names(f$estimate)
        bootCI <- c(median(resboot[-length(resboot[,1]),],na.rm=TRUE),
            quantile(resboot[-length(resboot[,1]),],0.025,na.rm=TRUE),
            quantile(resboot[-length(resboot[,1]),],0.975,na.rm=TRUE)) 
        names(bootCI) <- c("Median","2.5%","97.5%")
    }       
    return(structure(list(estim=estim,
        converg=t(resboot)[,length(resboot[,1])],CI=bootCI),
        class="bootdistcens"))
}

print.bootdistcens <- function(x,...){
    if (!inherits(x, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    cat("Parameter values obtained with nonparametric bootstrap \n")
    op<-options()
    options(digits=3)
    print(x$estim,...)    
    if (!is.null(x$converg)) { 
        nconverg<-length(x$converg[x$converg==0])
        cat("\n")
        cat("Maximum likelihood method converged for ",nconverg," among ",
        length(x$converg)," iterations \n")
    }
    options(op)

}

plot.bootdistcens <- function(x,...){
    if (!inherits(x, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    if (dim(x$estim)[1]==1) {
        stripchart(x$estim,method="jitter",
        xlab="Scatterplot of the boostrapped values of the parameter",...)
    }
    else {
        if (dim(x$estim)[2]==2)
            plot(x$estim,
            main="Scatterplot of the boostrapped values of the two parameters",...)
        else 
            plot(x$estim,
            main="Scatterplots of the boostrapped values of parameters",...)
    }
}

summary.bootdistcens <- function(object,...){
    if (!inherits(object, "bootdistcens"))
        stop("Use only with 'bootdistcens' objects")
    op<-options()
    options(digits=3)
    cat("Nonparametric bootstrap medians and 95% CI \n")
    print(object$CI)
    
     if (!is.null(object$converg)) { 
        nconverg<-length(object$converg[object$converg==0])
        cat(" \n")
        cat("Maximum likelihood method converged for ",nconverg," among ",
        length(object$converg)," iterations \n")
    }
   options(op)
}
