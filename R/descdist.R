descdist <- function(data,discrete=FALSE,boot=NULL,graph=TRUE)
{
    #if(is.mcnode(data)) data <- as.vector(data)
    if (missing(data) || !is.vector(data,mode="numeric"))
        stop("data must be a numeric vector")
    moment <- function(data,k){
        m1<-mean(data)
        return(sum((data-m1)^k)/length(data))
    }
    skewness <- function(data){
        sd<-sqrt(moment(data,2))
        return(moment(data,3)/sd^3)
    }
    kurtosis <- function(data){
        var<-moment(data,2)
        return(moment(data,4)/var^2)
    }
    
    res<-list(min=min(data),max=max(data),median=median(data),
    mean=mean(data),sd=sqrt(moment(data,2)),
    skewness=skewness(data),kurtosis=kurtosis(data))
    
    op<-options()
    options(digits=3)
    cat("summary statistics\n")
    cat("------\n")
    cat("min: ",res$min,"  max: ",res$max,"\n")
    cat("median: ",res$median,"\n")
    cat("mean: ",res$mean,"\n")
    cat("sample sd: ",res$sd,"\n")
    cat("sample skewness: ",res$skewness,"\n")
    cat("sample kurtosis: ",res$kurtosis,"\n")
    options(op)
    
    skewdata<-res$skewness
    kurtdata<-res$kurtosis
        
    # Cullen and Frey graph
    if (graph) {
        # bootstrap sample for observed distribution
        # and computation of kurtmax from this sample
        if (!is.null(boot)) {
            if (!is.numeric(boot) || boot<10) {
                stop("boot must be NULL or a integer above 10")
            }
            n<-length(data)

            databoot<-matrix(sample(data,size=n*boot,replace=TRUE),nrow=n,ncol=boot)
            s2boot<-sapply(1:boot,function(iter) skewness(databoot[,iter])^2)
            kurtboot<-sapply(1:boot,function(iter) kurtosis(databoot[,iter]))
            
            kurtmax<-max(10,ceiling(max(kurtboot)))
            xmax<-max(4,ceiling(max(s2boot)))
        } 
        else{
            kurtmax<-max(10,ceiling(kurtdata))
            xmax<-max(4,ceiling(skewdata^2))
        }

        ymax<-kurtmax-1
        plot(skewdata^2,kurtmax-kurtdata,pch=16,xlim=c(0,xmax),ylim=c(0,ymax),
        yaxt="n",xlab="square of skewness",ylab="kurtosis",main="Cullen and Frey graph")
        yax<-as.character(kurtmax-0:ymax)
        axis(side=2,at=0:ymax,labels=yax)
        if (!discrete) {
            # beta dist 
            p<-exp(-100)
            lq<-seq(-100,100,0.1)
            q<-exp(lq)
            s2a<-(4*(q-p)^2*(p+q+1))/((p+q+2)^2*p*q)
            ya<-kurtmax-(3*(p+q+1)*(p*q*(p+q-6)+2*(p+q)^2)/(p*q*(p+q+2)*(p+q+3)))
            p<-exp(100)
            lq<-seq(-100,100,0.1)
            q<-exp(lq)
            s2b<-(4*(q-p)^2*(p+q+1))/((p+q+2)^2*p*q)
            yb<-kurtmax-(3*(p+q+1)*(p*q*(p+q-6)+2*(p+q)^2)/(p*q*(p+q+2)*(p+q+3)))
            s2<-c(s2a,s2b)
            y<-c(ya,yb)
            polygon(s2,y,col="grey",border="grey80")
            # gamma dist
            lshape<-seq(-100,100,0.1)
            shape<-exp(lshape)
            s2<-4/shape
            y<-kurtmax-(3+6/shape)
            lines(s2,y,lty=2)
            # lnorm dist
            lshape<-seq(-100,100,0.1)
            shape<-exp(lshape)
            es2<-exp(shape^2)
            s2<-(es2+2)^2*(es2-1)
            y<-kurtmax-(es2^4+2*es2^3+3*es2^2-3)
            lines(s2,y,lty=3)
                
            legend(xmax*0.2,ymax*1.03,pch=16,legend="Observation",bty="n",cex=0.8,pt.cex=1.2)
            if (!is.null(boot)) {
            legend(xmax*0.2,ymax*0.98,pch=1,legend="bootstrapped values",
                bty="n",cex=0.8,col="blue")        
            }
            legend(xmax*0.55,ymax*1.03,legend="Theoretical distributions",bty="n",cex=0.8)
            legend(xmax*0.6,0.98*ymax,pch=8,legend="normal",bty="n",cex=0.8)
            legend(xmax*0.6,0.94*ymax,pch=2,legend="uniform",bty="n",cex=0.8)
            legend(xmax*0.6,0.90*ymax,pch=7,legend="exponential",bty="n",cex=0.8)
            legend(xmax*0.6,0.86*ymax,pch=3,legend="logistic",bty="n",cex=0.8)
            legend(xmax*0.6,0.82*ymax,fill="grey80",legend="beta",bty="n",cex=0.8)
            legend(xmax*0.6,0.78*ymax,lty=3,legend="lognormal",bty="n",cex=0.8)
            legend(xmax*0.6,0.74*ymax,lty=2,legend="gamma",bty="n",cex=0.8)
            legend(xmax*0.58,0.69*ymax,legend=c("(Weibull is close to gamma and lognormal)"),
            bty="n",cex=0.6)
        }
        else {         
            # negbin dist
            p<-exp(-10)
            lr<-seq(-100,100,0.1)
            r<-exp(lr)
            s2a<-(2-p)^2/(r*(1-p))
            ya<-kurtmax-(3+6/r+p^2/(r*(1-p)))
            p<-1-exp(-10)
            lr<-seq(100,-100,-0.1)
            r<-exp(lr)
            s2b<-(2-p)^2/(r*(1-p))
            yb<-kurtmax-(3+6/r+p^2/(r*(1-p)))
            s2<-c(s2a,s2b)
            y<-c(ya,yb)
            polygon(s2,y,col="grey80",border="grey80")
            legend(xmax*0.2,ymax*1.03,pch=16,legend="Observation",bty="n",cex=0.8,pt.cex=1.2)
            if (!is.null(boot)) {
            legend(xmax*0.2,ymax*0.98,pch=1,legend="bootstrapped values",
                bty="n",cex=0.8,col="blue")        
            }
            legend(xmax*0.55,ymax*1.03,legend="Theoretical distributions",bty="n",cex=0.8)
            legend(xmax*0.6,0.98*ymax,pch=8,legend="normal",bty="n",cex=0.8)
            legend(xmax*0.6,0.94*ymax,fill="grey80",legend="negative binomial",
                bty="n",cex=0.8)
           legend(xmax*0.6,0.90*ymax,lty=2,legend="Poisson",bty="n",cex=0.8)
            # poisson dist
            llambda<-seq(-100,100,0.1)
            lambda<-exp(llambda)
            s2<-1/lambda
            y<-kurtmax-(3+1/lambda)
            lines(s2,y,lty=2)
        }
        # bootstrap sample for observed distribution
        if (!is.null(boot)) {
            points(s2boot,kurtmax-kurtboot,pch=1,col="blue",cex=0.5)
        } 
        # observed distribution
        points(skewness(data)^2,kurtmax-kurtosis(data),pch=16,cex=2)
        # norm dist
        points(0,kurtmax-3,pch=8,cex=1.5,lwd=2)
        if (!discrete) {
            # unif dist   
            points(0,kurtmax-9/5,pch=2,cex=1.5,lwd=2)
            # exp dist
            points(2^2,kurtmax-9,pch=7,cex=1.5,lwd=2)
            # logistic dist
            points(0,kurtmax-4.2,pch=3,cex=1.5,lwd=2)
        }
    } # end of is (graph)
    
    invisible(res)
}
