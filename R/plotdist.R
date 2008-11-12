plotdist <- function(data,distr,para,breaks="default",discrete=FALSE,...){
    def.par <- par(no.readonly = TRUE)
    if (missing(data) || !is.vector(data,mode="numeric"))
        stop("data must be a numeric vector")
    if ((missing(distr) & !missing(para)) || 
    (missing(distr) & !missing(para)))
        stop("distr and para must defined")
    xlim<-c(min(data),max(data)) # for plot of discrete distributions
    if (missing(distr)) {
        par(mfrow=c(2,1))
        s<-sort(data)
        n<-length(data)
        obsp<-ecdf(s)(s)
        if (!discrete) {
            # plot for continuous data
            if (breaks=="default") 
                h<-hist(data,freq=FALSE,xlab="data",main=paste("Histogram"),...)
            else 
                h<-hist(data,freq=FALSE,xlab="data",main=paste("Histogram"),breaks=breaks,...)
            plot(s,obsp,main=paste("Cumulative distribution plot"),xlab="data",
            xlim=c(h$breaks[1],h$breaks[length(h$breaks)]),ylab="CDF",type='l',...)
        }
        else {
            # plot for discrete data
            if (breaks!="default") 
            warning("Breaks are not taken into account for discrete data")
            # plot of empirical distribution
            t<-table(data)
            xval<-as.numeric(names(t))
            xvalfin<-seq(min(xval),max(xval))
            ydobs<-as.vector(t)/n
            ydmax<-max(ydobs)
            plot(xval,ydobs,type='h',lwd=5,xlim=xlim,ylim=c(0,ydmax),
            main="Empirical distribution",xlab="data",ylab="Density",...)
            # plot of the cumulative probability distributions
            ycdfobs<-ecdf(data)(xvalfin)
            plot(xvalfin,ycdfobs,type='h',lwd=5,xlim=xlim,ylim=c(0,1),
            main="Empirical CDFs",xlab="data",ylab="CDF",...)
        }
    } #end of if (missing(distr))
    else {
        if (!is.character(distr)) distname<-substring(as.character(match.call()$distr),2)
            else distname<-distr
        if (!missing(discrete))
        warning("the argument discrete is not taken into account when distr is defined")
        if (!is.list(para)) 
        stop("'para' must be a named list")
        ddistname<-paste("d",distname,sep="")
        if (!exists(ddistname,mode="function"))
            stop(paste("The ",ddistname," function must be defined"))
        pdistname<-paste("p",distname,sep="")
        if (!exists(pdistname,mode="function"))
            stop(paste("The ",pdistname," function must be defined"))
        qdistname<-paste("q",distname,sep="")
        if (!exists(qdistname,mode="function"))
            stop(paste("The ",qdistname," function must be defined"))
        densfun <- get(ddistname,mode="function")    
        nm <- names(para)
        f <- formals(densfun)
        args <- names(f)
        m <- match(nm, args)
        if (any(is.na(m))) 
            stop(paste("'para' specifies names which are not arguments to ",ddistname))

        n<-length(data) 
        if (is.element(distname,c("binom","nbinom","geom","hyper","pois"))) 
            discrete<-TRUE
        else
            discrete<-FALSE
        if (!discrete) {
        # plot of continuous data with theoretical distribution
            par(mfrow=c(2,2))
            s<-sort(data)
            obsp<-ecdf(s)(s)
            theop<-do.call(pdistname,c(list(q=s),as.list(para)))
            # plot of the histogram with theoretical density
            if (breaks=="default")
                h<-hist(data,freq=FALSE,xlab="data",
                main=paste("Empirical and theoretical distr."),...)
            else
                h<-hist(data,freq=FALSE,xlab="data",
                main=paste("Empirical and theoretical distr."),breaks=breaks,...)           
            xhist<-seq(min(h$breaks),max(h$breaks),length=1000)
            yhist<-do.call(ddistname,c(list(x=xhist),as.list(para)))
            lines(xhist,yhist)
            # plot of the qqplot
            theoq<-do.call(qdistname,c(list(p=obsp),as.list(para)))
            plot(theoq,s,main=" QQ-plot",xlab="theoretical quantiles",ylab="sample quantiles",...)
            abline(0,1)
            # plot of the cumulative probability distributions
            plot(s,obsp,main=paste("Empirical and theoretical CDFs"),xlab="data",
            ylab="CDF",xlim=c(h$breaks[1],h$breaks[length(h$breaks)]),type='l',...)
            lines(s,theop,lty=2)
            legend(s[1],max(obsp),lty=c(1,2),legend=c("empirical",paste("theoretical")),bty='n')
            # plot of the ppplot
            plot(theop,obsp,main="PP-plot",xlab="theoretical probabilities",
            ylab="sample probabilities",...)
            abline(0,1)
        }
        else {
        # plot of discrete data with theoretical distribution
            par(mfrow=c(2,1))
            if (breaks!="default") 
            warning("Breaks are not taken into account for discrete distributions")
            # plot of empirical and theoretical distributions
            t<-table(data)
            xval<-as.numeric(names(t))
            xvalfin<-seq(min(xval),max(xval))
            yd<-do.call(ddistname,c(list(x=xvalfin),as.list(para)))
            ydobs<-as.vector(t)/n
            ydmax<-max(yd,ydobs)
            plot(xvalfin+0.1,yd,type='h',lwd=5,xlim=c(min(xval),max(xval)+0.1),ylim=c(0,ydmax),lty=2,
            main="Empirical and theoretical distr.",xlab="data",ylab="Density",col='red',...)
            points(xval,ydobs,type='h',lwd=5,lty=1)
            legend(xval[1]+0.8*(max(xval)-min(xval)),ydmax,lty=c(1,2),legend=c("empirical",paste("theoretical")),
            col=c('black','red'),bty='n',cex=0.8)
            # plot of the cumulative probability distributions
            ycdfobs<-ecdf(data)(xvalfin)
            ycdf<-do.call(pdistname,c(list(q=xvalfin),as.list(para)))
            plot(xvalfin+0.1,ycdf,type='h',lwd=5,xlim=c(min(xval),max(xval)+0.1),ylim=c(0,1),lty=2,
            main="Empirical and theoretical CDFs",xlab="data",ylab="CDF",col='red',...)
            points(xvalfin,ycdfobs,type='h',lwd=5,lty=1)
            legend(xval[1],1,lty=c(1,2),legend=c("empirical",paste("theoretical")),
            col=c('black','red'),bty='n',cex=0.8)
        }
    }
    par(def.par)    
}
