require("fitdistrplus")

nsample <- 10

# (1) Fit of a Weibull distribution to serving size data by maximum 
# goodness-of-fit estimation using all the distances available
# 

data(groundbeef)
serving <- groundbeef$serving
mgedist(serving, "weibull", gof="CvM")$estimate
mgedist(serving, "weibull", gof="CvM", silent = FALSE)$estimate
mgedist(serving, "weibull", gof="KS")$estimate
mgedist(serving, "weibull", gof="AD")$estimate
mgedist(serving, "weibull", gof="ADR")$estimate
mgedist(serving, "weibull", gof="ADL")$estimate
mgedist(serving, "weibull", gof="AD2R")$estimate
mgedist(serving, "weibull", gof="AD2L")$estimate
mgedist(serving, "weibull", gof="AD2")$estimate

#fitted values are around (2.123932 82.145177)

# (2) Fit of a uniform distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 

set.seed(1234)
u <- runif(nsample, min=5, max=10)
mgedist(u, "unif", gof="CvM")$estimate
mgedist(u, "unif", gof="KS")$estimate

#fitted values are around (5.495686 9.630573)

# (3) Fit of a triangular distribution using Cramer-von Mises or
# Kolmogorov-Smirnov distance
# 
require("mc2d")
set.seed(1234)
v <- rtriang(nsample, min=5, mode=6, max=10)
mgedist(v, "triang", start = list(min=4, mode=6,max=9), gof="CvM")$estimate
mgedist(v, "triang", start = list(min=4, mode=6,max=9), gof="KS")$estimate

#fitted values are around (4.896674 7.390675 8.494165 )

# (4) scaling problem
# the simulated dataset (below) has particularly small values, hence without scaling (10^0),
# the optimization raises an error. The for loop shows how scaling by 10^i
# for i=1,...,6 makes the fitting procedure work correctly.
if(FALSE)
{
  set.seed(1234)
  x2 <- rnorm(nsample, 1e-4, 2e-4)
  for(i in 6:0)
    cat(i, try(mgedist(x2*10^i,"cauchy")$estimate, silent=TRUE), "\n")
  
  # 6 19.86906 122.0622 
  # 5 1.986906 12.20622 
  # 4 0.1986906 1.220622 
  # 3 0.01986906 0.1220622 
  # 2 0.001986906 0.01220622 
  # 1 NA NA 
  # 0 NA NA 
}


# (5) scaling problem
#
if(FALSE)
{
  x <- c(-0.00707717, -0.000947418, -0.00189753, 
         -0.000474947, -0.00190205, -0.000476077, 0.00237812, 0.000949668, 
         0.000474496, 0.00284226, -0.000473149, -0.000473373, 0, 0, 0.00283688, 
         -0.0037843, -0.0047506, -0.00238379, -0.00286807, 0.000478583, 
         0.000478354, -0.00143575, 0.00143575, 0.00238835, 0.0042847, 
         0.00237248, -0.00142281, -0.00142484, 0, 0.00142484, 0.000948767, 
         0.00378609, -0.000472478, 0.000472478, -0.0014181, 0, -0.000946522, 
         -0.00284495, 0, 0.00331832, 0.00283554, 0.00141476, -0.00141476, 
         -0.00188947, 0.00141743, -0.00236351, 0.00236351, 0.00235794, 
         0.00235239, -0.000940292, -0.0014121, -0.00283019, 0.000472255, 
         0.000472032, 0.000471809, -0.0014161, 0.0014161, -0.000943842, 
         0.000472032, -0.000944287, -0.00094518, -0.00189304, -0.000473821, 
         -0.000474046, 0.00331361, -0.000472701, -0.000946074, 0.00141878, 
         -0.000945627, -0.00189394, -0.00189753, -0.0057143, -0.00143369, 
         -0.00383326, 0.00143919, 0.000479272, -0.00191847, -0.000480192, 
         0.000960154, 0.000479731, 0, 0.000479501, 0.000958313, -0.00383878, 
         -0.00240674, 0.000963391, 0.000962464, -0.00192586, 0.000481812, 
         -0.00241138, -0.00144963)
  
  
  #only i == 0, no scaling, should not converge.
  for(i in 6:0)
    cat(i, try(mgedist(x*10^i,"cauchy")$estimate, silent=TRUE), "\n")
  
  # 6 -283.1442 1185.903 
  # 5 -28.31442 118.5903 
  # 4 -2.831442 11.85903 
  # 3 -0.2831442 1.185903 
  # 2 -0.02831442 0.1185903 
  # 1 -0.002831442 0.01185903
  # 0 NA NA 
}


# (6) test error messages
#

dnorm2 <- pnorm2 <- function(x, a)
  "NA"
x <- rexp(10)

#should get a one-line error 
res <- mgedist(x, "norm2", start=list(a=1))
#as in 
attr(try(log("a"), silent=TRUE), "condition")


try(mgedist(c(serving, "a"), "gamma"))
try(mgedist(c(serving, NA), "gamma"))
try(mgedist(c(serving, Inf), "gamma"))
try(mgedist(c(serving, -Inf), "gamma"))
try(mgedist(c(serving, NaN), "gamma"))


# (7) test the component optim.message

x <- rnorm(1000)
#change parameter to obtain unsuccessful convergence
mgedist(x, "norm", control=list(maxit=2), start=list(mean=1e5, sd=1), optim.method="L-BFGS-B", lower=0)

# (8) test bounds

x <- rnorm(1000)
mgedist(x, "norm", optim.method="L-BFGS-B", lower=c(-Inf, 0)) #optim and L-BFGS-B
mgedist(x, "norm", optim.method="Nelder", lower=c(-Inf, 0))

# (9) numerical issue with log(1-p)

obs <- c(rep(1, 14), 10)
try(mgedist(obs, "weibull", gof="AD", control=list(trace=1, REPORT=1)))

# (10) large sample size issue

if(FALSE)
{
  set.seed(123)
  obs <- rlnorm(1e6, 3, 2)
  for(i in 2:6)
    cat(i, try(mgedist(obs[1:10^i], "lnorm", control=list(trace=0, REPORT=1))$estimate, silent=TRUE), "\n")
  
  # 2 3.143734 1.819549 
  # 3 3.023688 1.955416 
  # 4 2.995646 1.993392 
  # 5 3.002682 2.000394 
  # 6 2.999036 1.999887 
  # 7 3.000222 1.999981 
}

