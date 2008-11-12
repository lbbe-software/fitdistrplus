fitdistcens<-function (censdata, distr, start) 
{
    if (!is.character(distr)) distname<-substring(as.character(match.call()$distr),2)
    else distname<-distr
    ddistname<-paste("d",distname,sep="")
    if (!exists(ddistname,mode="function"))
        stop(paste("The ",ddistname," function must be defined"))
    pdistname<-paste("p",distname,sep="")
    if (!exists(pdistname,mode="function"))
        stop(paste("The ",pdistname," function must be defined"))
    # MLE fit with mledistcens 
    if (missing(start))
        mle<-mledistcens(censdata,distname) 
    else 
    mle<-mledistcens(censdata,distname,start)
    if (mle$convergence>0) 
        stop("the function mle failed to estimate the parameters, with the error code ",mle$convergence) 
    estimate<-mle$estimate
    sd<-sqrt(diag(solve(mle$hessian)))
    loglik<-mle$loglik
         
    return(structure(list(estimate = estimate, sd = sd, loglik = loglik, 
        censdata=censdata, distname=distname), class = "fitdistcens"))
        
}

print.fitdistcens <- function(x,...){
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    cat("Fitting of the distribution '",x$distname,"' on censored data by maximum likelihood \n")
    cat("Parameters:\n")
    op<-options()
    options(digits=3)
    print(data.frame("estimate" = x$estimate),...)
    options(op)

}

plot.fitdistcens <- function(x,...){
    if (!inherits(x, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    plotdistcens(censdata=x$censdata,distr=x$distname,
    para=as.list(x$estimate),...)
}

summary.fitdistcens <- function(object,...){
    if (!inherits(object, "fitdistcens"))
        stop("Use only with 'fitdistcens' objects")
    ddistname<-paste("d",object$distname,sep="")
    pdistname<-paste("p",object$distname,sep="")
    
    op<-options()
    options(digits=3)
    cat("FITTING OF THE DISTRIBUTION '",object$distname,"' BY MAXIMUM LIKELIHOOD ON CENSORED DATA \n")
    cat("PARAMETERS\n")
    print(cbind.data.frame("estimate" = object$estimate, "Std. Error" = object$sd))
    cat("Loglikelihood: ",object$loglik,"\n")
    options(op)
}
