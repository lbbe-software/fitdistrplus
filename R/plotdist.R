#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot, Jean-Baptiste Denis                                                                                                  
#                                                                                                                                                                        
#   This program is free software; you can redistribute it and/or modify                                               
#   it under the terms of the GNU General Public License as published by                                         
#   the Free Software Foundation; either version 2 of the License, or                                                   
#   (at your option) any later version.                                                                                                            
#                                                                                                                                                                         
#   This program is distributed in the hope that it will be useful,                                                             
#   but WITHOUT ANY WARRANTY; without even the implied warranty of                                          
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                 
#   GNU General Public License for more details.                                                                                    
#                                                                                                                                                                         
#   You should have received a copy of the GNU General Public License                                           
#   along with this program; if not, write to the                                                                                           
#   Free Software Foundation, Inc.,                                                                                                              
#   59 Temple Place, Suite 330, Boston, MA 02111-1307, USA                                                             
#                                                                                                                                                                         
#############################################################################
### plot functions for non-censored data
###
###         R functions
### 

plotdist <- function(data,distr,para,breaks="default",discrete=FALSE,...){
    def.par <- par(no.readonly = TRUE)
    if (missing(data) || !is.vector(data,mode="numeric"))
        stop("data must be a numeric vector")
    if ((missing(distr) & !missing(para)) || 
    (missing(distr) & !missing(para)))
        stop("distr and para must defined")
    xlim<-c(min(data),max(data)) # for plot of discrete distributions
    if (missing(distr)) { ## Plot of data only
        par(mfrow=c(2,1))
        s<-sort(data)
        n<-length(data)
        if (!discrete) {
            # plot for continuous data
            obsp <- ppoints(s)
            if (breaks=="default") 
                h<-hist(data,freq=FALSE,xlab="data",main=paste("Histogram"),...)
            else 
                h<-hist(data,freq=FALSE,xlab="data",main=paste("Histogram"),breaks=breaks,...)
            plot(s,obsp,main=paste("Cumulative distribution"),xlab="data",
            xlim=c(h$breaks[1],h$breaks[length(h$breaks)]),ylab="CDF",...)
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
            plot(xval,ydobs,type='h',xlim=xlim,ylim=c(0,ydmax),
            main="Empirical distribution",xlab="data",ylab="Density",...)
            # plot of the cumulative probability distributions
            ycdfobs<-ecdf(data)(xvalfin)
            plot(xvalfin,ycdfobs,type='h',xlim=xlim,ylim=c(0,1),
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

        n <- length(data) 
        if (is.element(distname,c("binom","nbinom","geom","hyper","pois"))) 
            discrete <- TRUE
        else
            discrete <- FALSE
        if (!discrete) {
        # plot of continuous data with theoretical distribution
            par(mfrow=c(2,2))
            s <- sort(data)
            obsp <- ppoints(s)
            theop <- do.call(pdistname,c(list(q=s),as.list(para)))
            # plot of the histogram with theoretical density
            # computes densities in order to define limits for y-axis
            if (breaks=="default")
                h <- hist(data,plot=FALSE)
            else
                h <- hist(data,breaks=breaks,plot=FALSE,...)           
            xhist <- seq(min(h$breaks),max(h$breaks),length=1000)
            yhist <- do.call(ddistname,c(list(x=xhist),as.list(para)))
            ymax <- ifelse(is.finite(max(yhist)),max(max(h$density),max(yhist)),max(h$density)) 
            # plot of histograms and theoretical density
            hist(data,freq=FALSE,xlab="data",ylim=c(0,ymax),breaks=h$breaks,
                main=paste("Empirical and theoretical distr."),...)
            lines(xhist,yhist)
           
            # plot of the qqplot
            theoq<-do.call(qdistname,c(list(p=obsp),as.list(para)))
            plot(theoq,s,main=" QQ-plot",xlab="theoretical quantiles",
            ylab="sample quantiles",...)
            abline(0,1)
            # plot of the cumulative probability distributions
            xmin<-h$breaks[1]
            xmax<-h$breaks[length(h$breaks)]
            plot(s,obsp,main=paste("Empirical and theoretical CDFs"),xlab="data",
            ylab="CDF",xlim=c(xmin,xmax),...)
            sfin<-seq(xmin,xmax,by=(xmax-xmin)/100)
            theopfin<-do.call(pdistname,c(list(q=sfin),as.list(para)))
            lines(sfin,theopfin,lty=1)
            #legend(s[1],max(obsp),lty=c(1,2),legend=c("empirical",paste("theoretical")),
            #bty='n')
            
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
            xlinesdec <- min((max(xval)-min(xval))/100,0.4)
            yd<-do.call(ddistname,c(list(x=xvalfin),as.list(para)))
            ydobs<-as.vector(t)/n
            ydmax<-max(yd,ydobs)
            plot(xvalfin+xlinesdec,yd,type='h',xlim=c(min(xval),max(xval)+xlinesdec),
                ylim=c(0,ydmax),lty=3,
                main="Empirical (full line) and theoretical (dotted line) distr.",xlab="data",
                ylab="Density",...)
            points(xval,ydobs,type='h',lty=1,...)
            #legend(xval[1]+0.8*(max(xval)-min(xval)),ydmax,lty=c(1,2),
            #    legend=c("empirical",paste("theoretical")),
            #    col=c('black','red'),bty='n',cex=0.8)
            
            # plot of the cumulative probability distributions
            ycdfobs<-ecdf(data)(xvalfin)
            ycdf<-do.call(pdistname,c(list(q=xvalfin),as.list(para)))
            plot(xvalfin+xlinesdec,ycdf,type='h',xlim=c(min(xval),max(xval)+xlinesdec),
                ylim=c(0,1),lty=3,
                main="Empirical (full line) and theoretical (dotted line) CDFs",xlab="data",
                ylab="CDF",...)
            points(xvalfin,ycdfobs,type='h',lty=1,...)
            #legend(xval[1],1,lty=c(1,2),legend=c("empirical",paste("theoretical")),
            #col=c('black','red'),bty='n',cex=0.8)
        }
    }
    par(def.par)    
}
