library(fitdistrplus)

ap <- cbind(L=c(0:14,1:15), 
            R=c(1:15, rep(Inf, 15)),
            count=c(456, 226, 152, 171, 135, 125,  83,  74,  51,  42,  
                    43,  34,  18,   9,   6,  39,  22,  23,  24,  107, 
                    133, 102,  68,  64,  45,  53,  33,  27,  23,  30))
npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=0, pkg="stats") 
npsurv.minimal(ap, w=1, maxit=100, tol=1e-6, verb=0, pkg="limSolve") 

