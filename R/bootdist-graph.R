
pairs4boot <- function(x, trueval, col4ramp = c("green", "yellow", "orange", "red"), 
                   nbgrid = 100, nbcol = 100, ...)
{
  x <- data.matrix(rbind(x, trueval))
  n <- NROW(x)
  if(is.null(trueval))
    id1 <- 1:n
  else
    id1 <- 1:(n-1)
  panel.upper <- function(x, y, ...)
  {
    points(x[id1], y[id1], xlim=range(x, na.rm=TRUE), ylim=range(y, na.rm=TRUE))
    if(!is.null(trueval))
      abline(v=x[n], h=y[n], col="red", lwd=2)
  }
  panel.lower <- function(x, y, ...)
  {
    id2 <- id1[!is.na(x[id1])]
    #require(MASS)
    k <- kde2d(x[id2], y[id2], n=nbgrid)
    image(k, col=colorRampPalette(col4ramp)(nbcol), add=TRUE, 
          xlim=range(x, na.rm=TRUE), ylim=range(y, na.rm=TRUE))
    if(!is.null(trueval))
      abline(v=x[n], h=y[n], col="black", lty="dashed")
  }
  pairs(x, upper.panel=panel.upper,
        lower.panel=panel.lower, ...)
  invisible()
}  
