require(fitdistrplus)

if(any(installed.packages()[, "Package"] == "actuar"))
{
  require(actuar)
  
data(endosulfan)
ATV <-endosulfan$ATV
library("actuar")
fBurr <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1, rate = 1))
llplot(fBurr)

fBurr2 <- fitdist(ATV, "burr", start = list(shape1 = 0.3, shape2 = 1),
                  fix.arg = list(rate = 1.5))
llplot(fBurr2)

fBurr3 <- fitdist(ATV, "burr", start = list(shape1 = 0.3, rate = 1),
                  fix.arg = list(shape2 = 1.5))
llplot(fBurr3)
}

### ADD tests