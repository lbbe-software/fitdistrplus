
library(CASdatasets)
data(frefictivetable)
frefictivetable$AgeOut <- as.numeric(frefictivetable$DateOut - frefictivetable$DateOfBirth)/365.25
frefictivetable$AgeIn <- as.numeric(frefictivetable$DateIn - frefictivetable$DateOfBirth)/365.25
frefictivetable$Death <- 1*(frefictivetable$Status == "deceased")

fre_male <- subset(frefictivetable, Gender == "Male")

rmcensor <- function(x, percent)
{
  if(percent <= 0)
    return(subset(x, Death == 1))
  nb_no <- sum(x$Death)
  nb_cens <- percent/(1-percent)*nb_no
  
  data_no <- subset(x, Death == 1)
  data_cens <- subset(x, Death == 0)
  rbind(data_no, data_cens[1:nb_cens,])
}
set.seed(123)
fre_male20 <- rmcensor(fre_male, .2)

mycol <- c("AgeIn", "AgeOut", "Death")
fre_male20final <- fre_male20[sample(1:NROW(fre_male20), 100), mycol]
rownames(fre_male20final) <- 1:100

fremale <- fre_male20final

save(fremale, file="./data/fremale.rda")
