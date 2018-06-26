#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller                                                                                                  
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

plotdist <- function(data, distr, para, histo = TRUE, breaks="default", demp = FALSE, discrete,  ...)
{
    def.par <- par(no.readonly = TRUE)
    if (missing(data) || !is.vector(data, mode="numeric"))
        stop("data must be a numeric vector")
    if ((missing(distr) & !missing(para)) || 
    (missing(distr) & !missing(para)))
        stop("distr and para must defined")
    if (!histo & !demp)
      stop("one the arguments histo and demp must be put to TRUE")
    xlim <- c(min(data), max(data)) # for plot of discrete distributions
    s <- sort(data)
    n <- length(data)
    
    
    if (missing(distr)) 
  { ## Plot of data only
        par(mfrow=c(1, 2))
		if(missing(discrete))
			discrete <- FALSE
      if (!discrete) 
      {
            # plot for continuous data alone
            obsp <- ppoints(s)
            # PLOT 1 ---
            if (histo)
            {
              if(demp)
              {
                if (breaks=="default") 
                  h <- hist(data, freq=FALSE, xlab="Data", main="Empirical density", ...)
                else 
                  h <- hist(data, freq=FALSE, xlab="Data", main="Empirical density", 
                            breaks=breaks, ...)
                lines(density(data)$x, density(data)$y, lty=2, col="black")
              }
              else
              {
                if (breaks=="default") 
                  h <- hist(data, freq=FALSE, xlab="Data", main="Histogram", ...)
                else 
                  h <- hist(data, freq=FALSE, xlab="Data", main="Histogram", 
                            breaks=breaks, ...)
              }
            }
            else
            {
              h <- hist(data, freq=FALSE, xlab="Data", main="Histogram", plot = FALSE, ...)
              plot(density(data)$x, density(data)$y, lty=1, col="black", type = "l",
                   xlab="Data", main=paste("Empirical density"), ylab = "Density",...)
            }
            
            # PLOT 2 ---
            plot(s, obsp, main=paste("Cumulative distribution"), xlab="Data", 
				      xlim=c(h$breaks[1], h$breaks[length(h$breaks)]), ylab="CDF", ...)
      }
      else 
      {
            # plot for discrete data alone
            if (breaks!="default") 
				    warning("Breaks are	not taken into account for discrete data")
            # plot of empirical distribution
            t <- table(data)
            xval <- as.numeric(names(t))
#            xvalfin <- seq(min(xval), max(xval),by=1)
            ydobs <- as.vector(t)/n
            ydmax <- max(ydobs)
            plot(xval, ydobs, type="h", xlim=xlim, ylim=c(0, ydmax), 
				      main="Empirical distribution", xlab="Data", ylab="Density", ...)
			
            # plot of the cumulative probability distributions
            ycdfobs <- cumsum(ydobs)
            plot(xval, ycdfobs, type="p", xlim=xlim, ylim=c(0, 1), 
				        main="Empirical CDFs", xlab="Data", ylab="CDF", ...)
        }
  } #end of if (missing(distr))
  else 
  { 
    # plot of data and distribution
      if (!is.character(distr)) 
			    distname <- substring(as.character(match.call()$distr), 2)
		  else 
			  distname <- distr
        if (!is.list(para)) 
			    stop("'para' must be a named list")
        ddistname <- paste("d", distname, sep="")
        if (!exists(ddistname, mode="function"))
            stop(paste("The ", ddistname, " function must be defined"))
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

        if(missing(discrete))
		    {
			    if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) 
				  discrete <- TRUE
			    else
				  discrete <- FALSE
		    }
		
        if (!discrete) 
        {
        # plot of continuous data with theoretical distribution
            par(mfrow=c(2, 2))
            obsp <- ppoints(s)
            
            # plot of empirical and theoretical density
            # computes densities in order to define limits for y-axis
            if (breaks=="default")
                h <- hist(data, plot=FALSE)
            else
                h <- hist(data, breaks=breaks, plot=FALSE, ...)           
            xhist <- seq(min(h$breaks), max(h$breaks), length=1000)
            yhist <- do.call(ddistname, c(list(xhist), as.list(para)))
            if(length(yhist) != length(xhist))
              stop("problem when computing densities.")
            ymax <- ifelse(is.finite(max(yhist)), max(max(h$density), max(yhist)), max(h$density)) 
            
            # PLOT 1 - plot of empirical and  theoretical density
            # empirical density
            if (histo)
            {
              hist(data, freq=FALSE, xlab="Data", ylim=c(0, ymax), breaks=h$breaks, 
                   main=paste("Empirical and theoretical dens."), ...)
              if(demp)
              {
                lines(density(data)$x, density(data)$y, lty=2, col="black")
              }
            }
            else
              plot(density(data)$x, density(data)$y, lty=2, col="black", type = "l",
                   xlab="Data", main=paste("Empirical and theoretical dens."),
                   ylab = "Density", xlim = c(min(h$breaks), max(h$breaks)), ...)
            if (demp)
              legend("topright",bty="n",lty=c(2,1),col=c("black","red"),
                     legend=c("empirical","theoretical"), bg="white",cex=0.7)
            
            # Add of theoretical density
            lines(xhist, yhist,lty=1,col="red")
           
            # PLOT 2 - plot of the qqplot
            theoq <- do.call(qdistname, c(list(obsp), as.list(para)))
            if(length(theoq) != length(obsp))
              stop("problem when computing quantities.")
            
            plot(theoq, s, main=" Q-Q plot", xlab="Theoretical quantiles", 
				      ylab="Empirical quantiles", ...)
            abline(0, 1)
              
            # PLOT 3 - plot of the cumulative probability distributions
            xmin <- h$breaks[1]
            xmax <- h$breaks[length(h$breaks)]
				    if(length(s) != length(obsp))
				      stop("problem when computing probabilities.")
				      
            plot(s, obsp, main=paste("Empirical and theoretical CDFs"), xlab="Data", 
				      ylab="CDF", xlim=c(xmin, xmax), ...)
            sfin <- seq(xmin, xmax, by=(xmax-xmin)/100)
            theopfin <- do.call(pdistname, c(list(sfin), as.list(para)))
            lines(sfin, theopfin, lty=1,col="red")
            
            # PLOT 4 - plot of the ppplot
            
				    theop <- do.call(pdistname, c(list(s), as.list(para)))
				    if(length(theop) != length(obsp))
				      stop("problem when computing probabilities.")
				    
				    plot(theop, obsp, main="P-P plot", xlab="Theoretical probabilities", 
				      ylab="Empirical probabilities", ...)
            abline(0, 1)
        }
        else 
        {
        # plot of discrete data with theoretical distribution
            par(mfrow=c(1, 2))
            if (breaks!="default") 
            warning("Breaks are not taken into account for discrete distributions")
            # plot of empirical and theoretical distributions
            t <- table(data)
            xval <- as.numeric(names(t))
            xvalfin <- seq(min(xval), max(xval),by=1)
            xlinesdec <- min((max(xval)-min(xval))/30, 0.4)
            yd <- do.call(ddistname, c(list(xvalfin), as.list(para)))
            if(length(yd) != length(xvalfin))
              stop("problem when computing density points.")
            
            ydobs <- as.vector(t)/n
            ydmax <- max(yd, ydobs)
            plot(xvalfin+xlinesdec, yd, type='h', xlim=c(min(xval), max(xval)+xlinesdec), 
                ylim=c(0, ydmax), lty=1, col="red",
                main="Emp. and theo. distr.", xlab="Data", 
                ylab="Density", ...)
            points(xval, ydobs, type='h', lty=1, col="black",...)
            legend("topright", lty=c(1, 1), col=c("black","red"),
                legend=c("empirical", paste("theoretical")), 
                bty="o", bg="white",cex=0.6,...)
            
            # plot of the cumulative probability distributions
            ycdf <- do.call(pdistname, c(list(xvalfin), as.list(para)))
            if(length(ycdf) != length(xvalfin))
              stop("problem when computing probabilities.")
            
			      plot(xvalfin, ycdf, type="s", xlim=c(min(xval), max(xval)+xlinesdec), 
				        ylim=c(0, 1), lty=1, col="red", 
				        main="Emp. and theo. CDFs", xlab="Data", 
				        ylab="CDF", ...)
			
#			plot(xvalfin+xlinesdec, ycdf, type="h", xlim=c(min(xval), max(xval)+xlinesdec), 
#               ylim=c(0, 1), lty=3, col="red", 
#               main="Emp. and theo. CDFs", xlab="Data", 
#               ylab="CDF", ...)
            ycdfobs <- cumsum(ydobs)
            points(xval,ycdfobs, type="p", col="black",...)
            legend("bottomright", lty=c(1, 1), col=c("black","red"), legend=c("empirical", paste("theoretical")), 
             bty="o", bg ="white",cex=0.6,...)
        }
    }
    par(def.par)    
    invisible()
}
