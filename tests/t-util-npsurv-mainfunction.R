library(fitdistrplus)

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

resap <- fitdistrplus:::npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 

#------------------------------------------------
# example with left censoring from package npsurv 
ap <- cbind(L=c(1:15, rep(-Inf, 15)), 
            R=c(1:15, 1:15),
            count=c(456, 226, 152, 171, 135, 125,  83,  74,  51,  42,  
                    43,  34,  18,   9,   6,  39,  22,  23,  24,  107, 
                    133, 102,  68,  64,  45,  53,  33,  27,  23,  30))
dim(ap)

resap <- fitdistrplus:::npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 


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

resap <- fitdistrplus:::npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 


#----------------------------------------------------------------------------------
# example with interval censoring leading to a single maximal intersection interval

LR <- matrix(1:100, ncol=2)
cnt <- round(1000*pgeom(1:NROW(LR)-1, prob=1/10, lower=FALSE))
fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)

fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x) # a single vector


resfk <- fitdistrplus:::npsurv.minimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 

#---------------------------------------------------------------------------------------
# geometric example with interval censoring leading to 50 maximal intersection intervals

LR <- matrix(1:100, ncol=2, byrow = TRUE)
theop <- dgeom(1:NROW(LR)-1, prob=1/10)
cnt <- round(1000*pgeom(1:NROW(LR)-1, prob=1/10, lower=FALSE))
fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)

fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x) #is diagonal


resfk <- fitdistrplus:::npsurv.minimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 

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


resfk <- fitdistrplus:::npsurv.minimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 

head(cbind(optimized=resfk$f$p, theoretical=theop[1:resfk$m]))


#----------------------------------------------------------------------------------------
# multimodal example with interval censoring leading to 43 maximal intersection intervals

n <- 100
set.seed(123)
LR <- sample.int(n)
LR <- cbind(LR, LR+sample.int(n)/10)
cnt <- round(450*dgeom(1:NROW(LR)-1, prob=1/10)
             +550*dbinom(1:NROW(LR)-1, size=NROW(LR), prob=4/5))
theop <- cnt/sum(cnt)
fakedata <- cbind(L=LR[,1], R=LR[,2], count=cnt)

fakedata.x2 <- fitdistrplus:::icendata(fakedata, w=1)
fakedata.x2.x <- rbind(cbind(fakedata.x2$t, fakedata.x2$t), fakedata.x2$o)
fakedata.x2.Delta <- fitdistrplus:::Deltamatrix(fakedata.x2.x)
str(fakedata.x2.Delta)

resfk <- fitdistrplus:::npsurv.minimal(fakedata, w=1, maxit=100, tol=1e-6, verb=mytrace, pkg="stats") 
