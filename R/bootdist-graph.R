#############################################################################
#   Copyright (c) 2009 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### bootstrap in fitdistrplus
###
###         R functions
### 

pairs4boot <- function(x, trueval, col4ramp = c("green", "yellow", "orange", "red"), 
                       nbgrid = 100, nbcol = 100, enhance=TRUE, ...)
{
  x <- data.matrix(rbind(x, trueval))
  n <- NROW(x)
  if(is.null(trueval))
    id1 <- 1:n
  else
    id1 <- 1:(n-1)
  panel.points <- function(x, y, ...)
  {
    points(x[id1], y[id1], xlim=range(x, na.rm=TRUE), ylim=range(y, na.rm=TRUE))
    if(!is.null(trueval))
      abline(v=x[n], h=y[n], col="red", lwd=2)
  }
  panel.density <- function(x, y, ...)
  {
    id2 <- id1[!is.na(x[id1])]
    #require(MASS)
    k <- kde2d(x[id2], y[id2], n=nbgrid)
    image(k, col=colorRampPalette(col4ramp)(nbcol), add=TRUE, 
          xlim=range(x, na.rm=TRUE), ylim=range(y, na.rm=TRUE))
    if(!is.null(trueval))
      abline(v=x[n], h=y[n], col="black", lty="dashed")
  }
  if(enhance)
    pairs(x, upper.panel=panel.points,
          lower.panel=panel.density, ...)
  else
    pairs(x, upper.panel=panel.points,
          lower.panel=panel.points, ...)
  invisible()
}  
