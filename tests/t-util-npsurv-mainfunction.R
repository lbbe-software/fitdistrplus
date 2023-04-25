library(fitdistrplus)
vizualise <- FALSE
mytrace <- 1
mytrace <- 0

#-------------------------------------------------
# example with right censoring from package npsurv 
ap <- cbind(L=c(1:15, 1:15), 
            R=c(1:15, rep(Inf, 15)),
            count=c(456, 226, 152, 171, 135, 125,  83,  74,  51,  42,  
                    43,  34,  18,   9,   6,  39,  22,  23,  24,  107, 
                    133, 102,  68,  64,  45,  53,  33,  27,  23,  30))
dim(ap)

resap <- fitdistrplus:::npsurvminimal(ap, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 

#------------------------------------------------
# example with left censoring from package npsurv 
ap <- cbind(L=c(1:15, rep(-Inf, 15)), 
            R=c(1:15, 1:15),
            count=c(456, 226, 152, 171, 135, 125,  83,  74,  51,  42,  
                    43,  34,  18,   9,   6,  39,  22,  23,  24,  107, 
                    133, 102,  68,  64,  45,  53,  33,  27,  23,  30))
dim(ap)

resap <- fitdistrplus:::npsurvminimal(ap, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 

if(vizualise)
{
  #-------------------------------------------------------------------------------------------------
  # example with interval censoring from package npsurv leading to 16 maximal intersection intervals  
  ap <- cbind(L=c(0:14,1:15), 
              R=c(1:15, rep(Inf, 15)),
              count=c(456, 226, 152, 171, 135, 125,  83,  74,  51,  42,  
                      43,  34,  18,   9,   6,  39,  22,  23,  24,  107, 
                      133, 102,  68,  64,  45,  53,  33,  27,  23,  30))
  dim(ap)
  
  ap.x2 <- fitdistrplus:::icendata(ap, w=1)
  ap.x <- rbind(cbind(ap.x2$t, ap.x2$t), ap.x2$o)
  ap.Delta <- fitdistrplus:::Deltamatrix(ap.x)
  #cbind(ap.Delta$left, unique(ap.x[,"L"]))
  
  resap <- fitdistrplus:::npsurvminimal(ap, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
  
  # if(FALSE)
  # {
  #   require(npsurv)
  #   rescheck <- npsurv::npsurv(ap, w=1, maxit=100, tol=1e-6, verb=3) 
  #   c(resap$ll, rescheck$ll)
  #   
  #   cbind(resap$f$left, rescheck$f$left)
  #   cbind(resap$f$right, rescheck$f$right)
  #   cbind(resap$f$p, rescheck$f$p)
  #   sum(abs(resap$f$p- rescheck$f$p))
  # }
  
  
  
  #----------------------------------------------------------------------------------
  # example with interval censoring leading to a single maximal intersection interval
  
  LR <- matrix(1:100, ncol=2)
  cnt <- round(1000*pgeom(1:NROW(LR)-1, prob=1/10, lower=FALSE))
  fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)
  
  fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
  fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
  fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x) # a single vector
  
  
  resfk <- fitdistrplus:::npsurvminimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
  
  #---------------------------------------------------------------------------------------
  # geometric example with interval censoring leading to 50 maximal intersection intervals
  
  LR <- matrix(1:100, ncol=2, byrow = TRUE)
  theop <- dgeom(1:NROW(LR)-1, prob=1/10)
  cnt <- round(1000*pgeom(1:NROW(LR)-1, prob=1/10, lower=FALSE))
  fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)
  
  fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
  fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
  fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x) #is diagonal
  
  
  resfk <- fitdistrplus:::npsurvminimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
  
  head(cbind(optimized=resfk$f$p, theoretical=theop))
  
  #----------------------------------------------------------------------------------------
  # multimodal example with interval censoring leading to 50 maximal intersection intervals
  
  LR <- matrix(1:100, ncol=2, byrow = TRUE)
  cnt <- round(450*dgeom(1:NROW(LR)-1, prob=1/10)
               +550*dbinom(1:NROW(LR)-1, size=NROW(LR), prob=4/5))
  theop <- cnt/sum(cnt)
  fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)
  
  fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
  fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
  fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x) #is diagonal
  
  
  resfk <- fitdistrplus:::npsurvminimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
  
  head(cbind(optimized=resfk$f$p, theoretical=theop[1:resfk$m]))
  
  
  #----------------------------------------------------------------------------------------
  # multimodal example with interval censoring leading to 43 maximal intersection intervals
  
  n <- 100
  set.seed(123)
  LR <- sample.int(n)
  LR <- cbind(LR, LR+sample.int(n)/10)
  cnt <- 1+round(450*dgeom(1:NROW(LR)-1, prob=1/10)
                 +550*dbinom(1:NROW(LR)-1, size=NROW(LR), prob=4/5))
  theop <- cnt/sum(cnt)
  fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)
  
  fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
  fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
  fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x)
  str(fakedata.x2.Delta)
  
  resfk <- fitdistrplus:::npsurvminimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
  
  
  # if(FALSE)
  # {
  #   require(npsurv)
  #   rescheck <- npsurv::npsurv(fakedata, w=1, maxit=100, tol=1e-6, verb=3)
  #   c(resfk$ll, rescheck$ll)
  # 
  #   cbind(resfk$f$left, rescheck$f$left)
  #   cbind(resfk$f$right, rescheck$f$right)
  #   cbind(resfk$f$p, rescheck$f$p)
  #   sum(abs(resfk$f$p- rescheck$f$p))
  # }
  
  
  #----------------------------------------------------------------------------------------
  # simulated example with interval censoring leading to 258 maximal intersection intervals
  
  set.seed(1232)
  ns <- 500
  ns <- 100
  r <- rnorm(ns)
  d8 <- data.frame(left = r, right = r)
  delta <- rlnorm(ns)
  icensored <- rbinom(ns, size = 1, prob = 0.2) 
  Lcensored <- rbinom(ns, size = 1, prob = 0.2*(1 - icensored))
  Rcensored <- rbinom(ns, size = 1, prob = 0.3*(1 - icensored)*(1 - Lcensored))
  # icensored +  Lcensored + Rcensored
  d8$left <- d8$left * (1 - Lcensored) + (-1000) * Lcensored
  d8$right <- d8$right * (1 - Rcensored) + (1000) * Rcensored
  d8$right <- d8$right + delta * icensored
  d8$right[d8$right == 1000] <- +Inf
  d8$left[d8$left == -1000] <- -Inf
  
  
  d8.x2 <- fitdistrplus:::icendata(d8, w=1)
  d8.x2.x <- rbind(cbind(d8.x2$t, d8.x2$t), d8.x2$o)
  d8.x2.Delta <- fitdistrplus:::Deltamatrix(d8.x2.x)
  str(d8.x2.Delta)
  
  system.time(
    resd8 <- fitdistrplus:::npsurvminimal(d8, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
  )
  
  # if(FALSE)
  # {
  #   require(npsurv)
  #   d8bis <- d8
  #   d8bis$left[is.na(d8bis$left)] <- -Inf
  #   d8bis$right[is.na(d8bis$right)] <- Inf
  #   
  #   rescheck <- npsurv::npsurv(d8bis, w=1, maxit=100, tol=1e-6, verb=3)
  #   c(resd8$ll, rescheck$ll)
  #   
  #   cbind(resd8$f$left, rescheck$f$left)
  #   cbind(resd8$f$right, rescheck$f$right)
  #   head(cbind(resd8$f$p, rescheck$f$p))
  #   sum(abs(resd8$f$p- rescheck$f$p))
  # }
  
  
  #------------------------------------------------------
  # crash example with wrong interval censoring intervals
  
  set.seed(1232)
  ns <- 100
  r <- rnorm(ns)
  d8 <- data.frame(left = r, right = r)
  delta <- rlnorm(ns)
  icensored <- rbinom(ns, size = 1, prob = 0.2) 
  Lcensored <- rbinom(ns, size = 1, prob = 0.2*(1 - icensored))
  Rcensored <- rbinom(ns, size = 1, prob = 0.3*(1 - icensored)*(1 - Lcensored))
  # icensored +  Lcensored + Rcensored
  d8$left <- d8$left * (1 - Lcensored) + (-1000) * Lcensored
  d8$right <- d8$right * (1 - Rcensored) + (1000) * Rcensored
  d8$right <- d8$right + delta * icensored
  d8$right[d8$right == 1000] <- -Inf
  d8$left[d8$left == -1000] <- +Inf
  
  try(resd8 <- fitdistrplus:::npsurvminimal(d8, w=1, maxit=100, tol=1e-6, verb=2, pkg="stats"))
  
}