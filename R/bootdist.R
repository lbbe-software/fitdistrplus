bootdist<-function (f, bootmethod="param", niter=999)
{ 
    if (niter<10) 
        stop("niter must be an integer above 10")
    if (!is.element(bootmethod,c("param","nonparam")))
        stop("bootmethod must be affected to 'param' or 'nonparam'") 
    if (!inherits(f, "fitdist"))
        stop("Use only with 'fitdist' objects")
    rdistname<-paste("r",f$distname,sep="")
    if (!exists(rdistname,mode="function"))
        stop(paste("The ",rdistname," function must be defined"))
    if (bootmethod=="param") { # parametric bootstrap
        rdata<-do.call(rdistname,c(list(n=niter*f$n),as.list(f$estimate)))
        dim(rdata)<-c(f$n,niter)
    }
    else { # non parametric bootstrap
        rdata<-sample(f$data,size=niter*f$n,replace=TRUE)
        dim(rdata)<-c(f$n,niter)
    }
    if (f$method=="mle") {
        start<-f$estimate
        funcmle<-function(iter) {
            mle<-mledist(rdata[,iter],f$distname,start)
            return(c(mle$estimate,mle$convergence))
        }
        resboot<-sapply(1:niter,funcmle)
        rownames(resboot)<-c(names(start),"convergence")
        if (length(resboot[,1])>2) 
            bootCI <- cbind(apply(resboot[-length(resboot[,1]),],1,median,na.rm=TRUE),
            apply(resboot[-length(resboot[,1]),],1,quantile,0.025,na.rm=TRUE),
            apply(resboot[-length(resboot[,1]),],1,quantile,0.975,na.rm=TRUE))
        else
            bootCI <- cbind(median(resboot[-length(resboot[,1]),],na.rm=TRUE),
            quantile(resboot[-length(resboot[,1]),],0.025,na.rm=TRUE),
            quantile(resboot[-length(resboot[,1]),],0.975,na.rm=TRUE))        
        colnames(bootCI) <- c("Median","2.5%","97.5%")
        return(structure(list(bootestim=t(resboot)[,-length(resboot[,1])],
        bootconverg=t(resboot)[,length(resboot[,1])],bootmethod=bootmethod, bootCI=bootCI),
        class="bootdist"))
    }
    else { # f$method=="mom"
        funcmom<-function(iter) {
            mom<-momdist(rdata[,iter],f$distname)
        }
        resboot<-sapply(1:niter,funcmom)
        if (is.vector(resboot)) {
            bootCI <- cbind(median(resboot,na.rm=TRUE),quantile(resboot,0.025,na.rm=TRUE),
            quantile(resboot,0.975,na.rm=TRUE))
            colnames(bootCI)<-c("Median","2.5%","97.5%") # A MODIFIER POUR BON AFFICHAGE
        }           
         else {
            bootCI <- cbind(apply(resboot,1,median,na.rm=TRUE),apply(resboot,1,quantile,0.025,na.rm=TRUE),
            apply(resboot,1,quantile,0.975,na.rm=TRUE))
            colnames(bootCI) <- c("Median","2.5%","97.5%")
         }
        return(structure(list(bootestim=t(resboot),
        bootconverg=NULL,bootmethod=bootmethod,bootCI=bootCI),
        class="bootdist"))
    }
}

print.bootdist <- function(x,...){
    if (!inherits(x, "bootdist"))
        stop("Use only with 'bootdist' objects")
    if (x$bootmethod=="param") 
        cat("Parameter values obtained with parametric bootstrap \n")
    else
       cat("Parameter values obtained with nonparametric bootstrap \n")
    op<-options()
    options(digits=3)
    print(x$bootestim,...)    
    if (!is.null(x$bootconverg)) { 
        nconverg<-length(x$bootconverg[x$bootconverg==0])
        cat("\n")
        cat("Maximum likelihood method converged for ",nconverg," among ",length(x$bootconverg)," iterations \n")
    }
    options(op)

}

plot.bootdist <- function(x,...){
    if (!inherits(x, "bootdist"))
        stop("Use only with 'bootdist' objects")
    if (dim(x$bootestim)[1]==1) {
        stripchart(x$bootestim,method="jitter",xlab="Scatterplot of the boostrapped values of the parameter",...)
    }
    else {
        if (dim(x$bootestim)[2]==2)
            plot(x$bootestim[,1],x$bootestim[,2],xlab=colnames(x$bootestim)[1],ylab=colnames(x$bootestim)[2],
            main="Scatterplot of the boostrapped values of the two parameters",...)
        else 
            pairs(x$bootestim,main="Scatterplots of the boostrapped values of parameters",...)
    }
}

summary.bootdist <- function(object,...){
    if (!inherits(object, "bootdist"))
        stop("Use only with 'bootdist' objects")
    op<-options()
    options(digits=3)
    if (object$bootmethod=="param") 
        cat("Parametric bootstrap medians and 95% CI \n")
    else
       cat("Nonparametric bootstrap medians and 95% CI \n")
    print(object$bootCI)
    
     if (!is.null(object$bootconverg)) { 
        nconverg<-length(object$bootconverg[object$bootconverg==0])
        cat(" \n")
        cat("Maximum likelihood method converged for ",nconverg," among ",length(object$bootconverg)," iterations \n")
    }
   options(op)
}
