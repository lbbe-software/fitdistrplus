#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Regis Pouillot                                                                                                 
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

plotdistcens <- function(censdata, distr, para, leftNA = -Inf,rightNA = Inf, NPMLE = TRUE,
                         Turnbull = NULL, Turnbull.confint = FALSE, NPMLE.method = "Wang", ...){
  if (missing(censdata) ||
      !(is.vector(censdata$left) & is.vector(censdata$right) & length(censdata[,1])>1))
    stop("datacens must be a dataframe with two columns named left 
         and right and more than one line")
  if ((missing(distr) & !missing(para)) || 
      (missing(distr) & !missing(para)))
    stop("distr and para must defined")
  if (missing(distr))
  {
    onlyCDFplot <- TRUE
    titleCDF <- "Cumulative distribution"
  } else 
  {
    if (!is.character(distr)) 
      distname <- substring(as.character(match.call()$distr), 2)
    else 
      distname <- distr
    if (!is.list(para)) 
      stop("'para' must be a named list")
    ddistname <- paste("d", distname, sep="")
    if (!exists(ddistname, mode="function"))
      stop(paste("The ", ddistname, " function must be defined"))
    
    onlyCDFplot <- FALSE
    titleCDF <- "Empirical and theoretical CDFs"
  }
  
  ###################################### Part to check !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (!missing(Turnbull))
  {
    warning("The argument Turnbull is deprecated and should note be used any more. 
            Now use the argument NPMLE to tell if you want to compute a nonparametric
            maximum likelihood estimation of the cumulative distribution, and the argument NPMLE.method
            to define the method chosen for the computation (Turnbull or Wang).") 
    if (missing(NPMLE) & missing(NPMLE.method)) 
    {
      if (Turnbull == TRUE)
      {
        NPMLE <- TRUE
        NPMLE.method <- "Turnbull"
      } else
      {
        NPMLE <- FALSE
      }
    }
  }
  
  if ((Turnbull.confint == TRUE) & (NPMLE.method == "Wang"))
  {
    warning("When Turnbull.confint is TRUE NPMLE.method is forced to Turnbull." )
    NPMLE.method <- "Turnbull"
    # so the second part of the message will be printed in the following if needed
    onlyCDFplot <- TRUE
  } 
  if ((NPMLE.method == "Turnbull") & !missing(distr))
  {
    warning("Q-Q plot and P-P plot are available only using the method implemented in the package npsurv (Wang) 
            with the arguments NPMLE.method at Wang (default recommended arguments)." )
    onlyCDFplot <- TRUE
  }
  if ((NPMLE == FALSE) & !missing(distr))
  {
    warning("When NPMLE is FALSE the nonparametric maximum likelihood estimation 
            of the cumulative distribution function is not computed.
            Q-Q plot and P-P plot are available only using the method implemented in the package npsurv (Wang) 
            with the arguments Turnbull.confint at FALSE and NPMLE.method at Wang (default recommended arguments)." )
    onlyCDFplot <- TRUE
  }
  if (!onlyCDFplot)
  {
    pdistname <- paste("p", distname, sep="")
    if (!exists(pdistname, mode="function"))
      stop(paste("The ", pdistname, " function must be defined"))
    qdistname <- paste("q", distname, sep="")
    if (!exists(qdistname, mode="function"))
      stop(paste("The ", qdistname, " function must be defined"))
    densfun <- get(ddistname, mode="function")    
    nm <- names(para)
    f <- formals(densfun)
    args <- names(f)
    m <- match(nm, args)
    if (any(is.na(m))) 
      stop(paste("'para' specifies names which are not arguments to ", ddistname))
    
    def.par <- par(no.readonly = TRUE)
    par(mfrow = c(2, 2))
  }
  
  # Plot of the empirical distribution as an ECDF       
  if (NPMLE)
  {
    if (NPMLE.method == "Wang") # plot using package npsurv
    {
      db <- censdata
      db$left[is.na(db$left)] <- -Inf
      db$right[is.na(db$right)] <- Inf
      f <- npsurv(db)$f
      bounds <- c(f$right, f$left)
      finitebounds <- bounds[is.finite(bounds)]
      upper <- max(finitebounds)
      lower <- min(finitebounds)
      width <- upper - lower
      xmax <- upper + width * 0.1 # limits for the plot
      xmin <- lower - width * 0.1
      xmaxinf <- xmax + width * 10 # larger limits for the non finite rectangles
      xmininf <- xmin - width * 10
      k <- length(f$left)
      Fnpsurv <- cumsum(f$p) 
      
      ## calul des points points pour Q et P dans les GOF stat et graph
      Fbefore <- c(0, Fnpsurv[-k])
      df <- data.frame(left = f$left, right = f$right)
      dfb <- df[(df$left != -Inf) & (df$right != Inf), ]
      
      # Definition of vertices of each rectangle
      Qi.left <- df$left # dim k
      Qi.left4plot <- Qi.left
      if (Qi.left4plot[1] == - Inf) Qi.left4plot[1] <- xmininf
      Qi.right <- df$right
      Qi.right4plot <- Qi.right
      if (Qi.right4plot[k] == Inf) Qi.right4plot[k] <- xmaxinf
      Pi.low <- Fbefore
      Pi.up <- Fnpsurv
      
      # Plot of the ECDF
      plot(1, 1, type = "n", xlim = c(xmin, xmax), ylim = c(0,1), xlab = "Censored data", 
           ylab = "CDF", main = titleCDF)
      
      # the line at right of the rectangles
      dright <- c(f$left[1], rep(f$right, rep(2,k)), f$right[k]) 
      Fright <- rep(c(0,Fnpsurv), rep(2,k+1))
      lines(dright, Fright)
      ### the line at left of the rectangles
      dleft = rep(c(f$left,f$right[k]), rep(2,k+1))
      Fleft = c(0,rep(Fnpsurv, rep(2,k)),1)
      lines(dleft, Fleft)
      
      # Add of the filled rectangles
      for(i in 1:k) {
        rect(xleft = Qi.left4plot, ybottom = Pi.low, xright = Qi.right4plot, ytop = Pi.up,
             border = "black", col = "lightgrey")
      }
    } else # plot using package survival
    {
      survdata <- Surv(time = censdata$left, time2 = censdata$right, type="interval2")
      survfitted <- survfit(survdata ~ 1)
      if (Turnbull.confint)
        plot(survfitted,fun="event",xlab="Censored data",
             ylab="CDF",main = titleCDF, ...)
      else
        plot(survfitted,fun="event",xlab="Censored data",
             ylab="CDF", conf.int = FALSE, main = titleCDF, ...)
      xmin <- par("usr")[1]
      xmax <- par("usr")[2]
    }
  }else # if !NPMLE
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
    xmin<-min(xminright,xminleft)
    xmaxright<-max(censdata[!is.na(censdata$right),]$right)
    xmaxleft<-max(censdata[!is.na(censdata$left),]$left)
    xmax<-max(xmaxright,xmaxleft)
    xrange<-xmax-xmin
    xmin<-xmin-0.3*xrange
    xmax<-xmax+0.3*xrange
    xlim<-c(xmin,xmax)
    plot(c(0,0),c(0,0),type="n",xlim=xlim,ylim=c(0,1),xlab="Censored data",
         ylab="CDF",main="Cumulative distribution", ...)
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
  } # en of else if NPMLE
  
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
    lines(s,theop,col="red")
  }
  if (!onlyCDFplot)
  {
    ## Q-Q plot
    plot(1, 1, type = "n", main = "Q-Q plot", xlim = c(xmin, xmax), ylim = c(xmin, xmax),
         xlab = "Theoretical quantiles", ylab = "Empirical quantiles")
    
    # plot of rectangles
    Qitheo.left <- do.call(qdistname, c(list(p = Pi.low), as.list(para)))
    Qitheo.right <- do.call(qdistname, c(list(p = Pi.up), as.list(para)))
    Qitheo.left4plot <- Qitheo.left
    if (Qitheo.left4plot[1] == - Inf) Qitheo.left4plot[1] <- xmininf
    Qitheo.right4plot <- Qitheo.right
    if (Qitheo.right4plot[k] == Inf) Qitheo.right4plot[k] <- xmaxinf
    rect(xleft = Qitheo.left4plot, ybottom = Qi.left4plot, xright = Qitheo.right4plot, ytop = Qi.right4plot, 
         border = "black", col = "lightgrey")
    abline(0,1)
    
    
    ## P-P plot
    plot(1, 1, type = "n", main = "P-P plot", xlim = c(0,1), ylim = c(0, 1),
         xlab = "Theoretical probabilities", ylab = "Empirical probabilities")
    
    # plot of rectangles
    Pitheo.low <- do.call(pdistname, c(list(q=Qi.left), as.list(para)))
    Pitheo.up <- do.call(pdistname, c(list(q=Qi.right), as.list(para)))
    rect(xleft = Pitheo.low, ybottom = Pi.low, xright = Pitheo.up, ytop = Pi.up, 
         border = "black", col = "lightgrey")
    abline(0,1)
    
    par(def.par)    
  }
  invisible()
}
