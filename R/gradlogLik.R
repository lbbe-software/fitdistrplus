#############################################################################
#   Copyright (c) 2016 Marie Laure Delignette-Muller, Christophe Dutang                                                                                                  
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
### gradient of log-likelihood
###
###         R functions (not exported)
### 

#Beta distribution
grdbeta <- function(x, d1, d2) #individual contribution
  c(log(x)-digamma(d1)+digamma(d1+d2), log(1-x)-digamma(d2)+digamma(d1+d2))
grlnlbeta <- function(par, obs, ...) #total grad loglik
  -rowSums(sapply(obs, function(x) grdbeta(x, d1=par[1], d2=par[2])))


#Gamma distribution
grdgamma <- function(x, shape, rate) #individual contribution
  c(log(x)-log(rate)-digamma(shape), x/rate^2-shape/rate)
grlnlgamma <- function(par, obs, ...) #total grad loglik
{
  n <- length(obs)
  res <- grdgamma(obs, shape=par[1], rate=par[2])
  c(-sum(res[1:n]), -sum(res[1:n+n]))
}  
