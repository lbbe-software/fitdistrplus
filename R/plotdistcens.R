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
### plot functions for censored data
###
###         R functions
### 

plotdistcens <- function(censdata,distr,para,leftNA = -Inf,rightNA = Inf,Turnbull = TRUE,
Turnbull.confint = FALSE, ...){
    def.par <- par(no.readonly = TRUE)
    if (missing(censdata) ||
        !(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[,1])>1))
        stop("datacens must be a dataframe with two columns named left 
            and right and more than one line")
    if ((missing(distr) & !missing(para)) || 
    (missing(distr) & !missing(para)))
        stop("distr and para must defined")
        
    if (Turnbull)
    {
        survdata <- Surv(time = censdata$left, time2 = censdata$right, type="interval2")
        survfitted <- survfit(survdata ~ 1)
        if (Turnbull.confint)
            plot(survfitted,fun="event",xlab="censored data",
            ylab="CDF",main="Cumulative distribution",...)
        else
            plot(survfitted,fun="event",xlab="censored data",
            ylab="CDF",main="Cumulative distribution", conf.int = FALSE, ...)
        xmin <- par("usr")[1]
        xmax <- par("usr")[2]
    }
    
   
    else
    {    
        if (is.finite(leftNA) & any(is.na(censdata$left)))
            censdata[is.na(censdata$left),]$left<-leftNA
        if (is.finite(rightNA) & any(is.na(censdata$right)))
            censdata[is.na(censdata$right),]$right<-rightNA
        lcens<-censdata[is.na(censdata$left),]$right
        if (any(is.na(lcens)) )
            stop("An observation cannot be both right and left censored, coded with two NA values")
        rcens<-censdata[is.na(censdata$right),]$left
        noricens<-censdata[!is.na(censdata$left) & !is.na(censdata$right),]
        # definition of mid point for each observation (if not NA) 
        # in order to have the order of plot of each observation
        # and order of left and rigth bounds for censored observations
        midnoricens<-(noricens$left+noricens$right)/2
        ordmid<-order(midnoricens)
        ordlcens<-order(lcens)
        ordrcens<-order(rcens)
        nlcens<-length(lcens)
        nrcens<-length(rcens)
        nnoricens<-length(noricens$left)
        n<-length(censdata$left)
        # definition of xlim 
        xminright<-min(censdata[!is.na(censdata$right),]$right)
        xminleft<-min(censdata[!is.na(censdata$left),]$left)
        xmin=min(xminright,xminleft)
        xmaxright<-max(censdata[!is.na(censdata$right),]$right)
        xmaxleft<-max(censdata[!is.na(censdata$left),]$left)
        xmax=max(xmaxright,xmaxleft)
        xrange<-xmax-xmin
        xmin<-xmin-0.3*xrange
        xmax<-xmax+0.3*xrange
        xlim<-c(xmin,xmax)
        plot(c(0,0),c(0,0),type="n",xlim=xlim,ylim=c(0,1),xlab="censored data",
        ylab="CDF",main="Cumulative distribution",...)
        # functions to plot one interval or point for each observation for 
        # observation ordered i out of n
        plotlcens<-function(i) {
            y<-i/n
            lines(c(xmin,lcens[ordlcens[i]]),c(y,y))
        }
        if (nlcens>=1)
            toto<-sapply(1:nlcens,plotlcens)
        plotnoricens<-function(i) {
            y<-(i+nlcens)/n
            if (noricens[ordmid[i],]$left!=noricens[ordmid[i],]$right)
                lines(c(noricens[ordmid[i],]$left,noricens[ordmid[i],]$right),c(y,y))
            else
                points(noricens[ordmid[i],]$left,y,pch=4)
        }
        if (nnoricens>=1)
            toto<-sapply(1:nnoricens,plotnoricens)
        plotrcens<-function(i) {
            y<-(i+nlcens+nnoricens)/n
            lines(c(rcens[ordrcens[i]],xmax),c(y,y))
        }
        if (nrcens>=1)
            toto<-sapply(1:nrcens,plotrcens)
    } # else if Turnbull
    
    if (!missing(distr)){ # plot of the theoretical cumulative function
        if (!is.character(distr)) distname<-substring(as.character(match.call()$distr),2)
        else distname<-distr
        if (!is.list(para)) 
            stop("'para' must be a named list")
        ddistname<-paste("d",distname,sep="")
        if (!exists(ddistname,mode="function"))
            stop(paste("The ",ddistname," function must be defined"))
        pdistname<-paste("p",distname,sep="")
        if (!exists(pdistname,mode="function"))
            stop(paste("The ",pdistname," function must be defined"))
        densfun <- get(ddistname,mode="function")    
        nm <- names(para)
        f <- formals(densfun)
        args <- names(f)
        m <- match(nm, args)
        if (any(is.na(m))) 
            stop(paste("'para' specifies names which are not arguments to ",ddistname))
        # plot of continuous data with theoretical distribution
        s<-seq(xmin,xmax,by=(xmax-xmin)/100)
        theop<-do.call(pdistname,c(list(q=s),as.list(para)))
        lines(s,theop,...)
     }
    par(def.par)    
    
}
