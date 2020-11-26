#REPORT ON WINE DATASET, LIBRARY(CONTAMINATEDMIXT)
library(ContaminatedMixt)
data(wine)
dim(wine) 
attach(wine)
str(wine)

#UNIVARIATE ANALYSIS
#TYPE
table(Type)
frequencyType <- table(Type)
prop.table(frequencyType)*100 #*100 for percentage

#ALCOHOL
summary(Alcohol)
frequencyAlcohol <- table(Alcohol)
frequencyAlcohol
length(frequencyAlcohol)
names(frequencyAlcohol)[frequencyAlcohol == max(frequencyAlcohol)] #mode
library(labstatR)
labstatR::cv(Alcohol)
library(EnvStats)
kurtosis(Alcohol) #excess kurtosis
kurtosis(Alcohol, excess=FALSE)
skew(Alcohol)
boxplot(Alcohol, main = "Boxplot of Alcohol")
boxplot(Alcohol)$out
hist(Alcohol, breaks = 126, freq = FALSE)
lines(density(Alcohol), col = "red")

library(gamlss)
set.seed(123)
fit.GA.2 <- gamlssMXfits(n = 5, Alcohol~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.2[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.2[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.2[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.2[["models"]][[2]][["sigma.coefficients"]])
hist(Alcohol, breaks = 126, freq = FALSE)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.2[["prob"]][1]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.2[["prob"]][2]*dGA(seq(min(Alcohol),
      max(Alcohol),length=length(Alcohol)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.2[["prob"]][1]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      mu = mu.hat1, sigma = sigma.hat1) +
      fit.GA.2[["prob"]][2]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.2)

library(gamlss)
set.seed(123)
fit.GA.3 <- gamlssMXfits(n = 5, Alcohol~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.3[["models"]][[3]][["mu.coefficients"]])    
sigma.hat3 <- exp(fit.GA.3[["models"]][[3]][["sigma.coefficients"]])
hist(Alcohol, breaks = 126, freq = FALSE)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.3[["prob"]][1]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
                                mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.3[["prob"]][2]*dGA(seq(min(Alcohol),
                                    max(Alcohol),length=length(Alcohol)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.3[["prob"]][3]*dGA(seq(min(Alcohol),
                                    max(Alcohol),length=length(Alcohol)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
      fit.GA.3[["prob"]][1]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
                                mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.3[["prob"]][2]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
                                   mu = mu.hat2, sigma = sigma.hat2) +
      fit.GA.3[["prob"]][3]*dGA(seq(min(Alcohol),max(Alcohol),length=length(Alcohol)),
                                mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.3)

data.frame(row.names=c("Gamma mixture with K=2", "Gamma mixture with K=3"), AIC=c(fit.GA.2$aic,fit.GA.3$aic), SBC=c(fit.GA.2$sbc,fit.GA.3$sbc))

#MALIC
summary(Malic)
frequencyMalic <- table(Malic)
frequencyMalic
length(frequencyMalic)
names(frequencyMalic)[frequencyMalic==max(frequencyMalic)]
labstatR::cv(Malic)
library(EnvStats)
kurtosis(Malic) #excess kurtosis
kurtosis(Malic, excess=FALSE)
skew(Malic)
boxplot(Malic, main = "Boxplot of Malic")
boxplot(Malic)$out
hist(Malic, breaks = 133, freq = FALSE)
lines(density(Malic), col = "red")

library(gamlss)
fit.EXP <- histDist(Malic, family=EXP, nbins = 133, main="Exponential distribution")
fit.GA <- histDist(Malic, family=GA, nbins = 133, main="Gamma distribution")
fit.IG <- histDist(Malic, family=IG, nbins = 133, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Malic, family=LOGNO, nbins = 133, main="Log-Normal distribution")
fit.WEI <- histDist(Malic, family=WEI, nbins = 133, main="Weibull distribution")
fit.LO <- histDist(Malic, family=LO, nbins=133, main="Logistic distribution")
fit.sep2 <-histDist(Malic,family=SEP2, nbins=133, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Malic, family=ST5, nbins=133, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                   AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                   fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Mal <- gamlssMXfits(n = 5, Malic~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Mal[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Mal[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Mal[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Mal[["models"]][[2]][["sigma.coefficients"]])
hist(Malic, breaks = 133, freq = FALSE)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal[["prob"]][1]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                  mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal[["prob"]][2]*dGA(seq(min(Malic),
                                      max(Malic),length=length(Malic)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal[["prob"]][1]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                  mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Mal[["prob"]][2]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                     mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Mal)

library(gamlss)
set.seed(123)
fit.GA.Mal3 <- gamlssMXfits(n = 5, Malic~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Mal3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Mal3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Mal3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Mal3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Mal3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Mal3[["models"]][[3]][["sigma.coefficients"]])
hist(Malic, breaks = 133, freq = FALSE)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal3[["prob"]][1]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                   mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal3[["prob"]][2]*dGA(seq(min(Malic),
                                       max(Malic),length=length(Malic)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal3[["prob"]][2]*dGA(seq(min(Malic),
                                       max(Malic),length=length(Malic)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Malic),max(Malic),length=length(Malic)),
      fit.GA.Mal3[["prob"]][1]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                   mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Mal3[["prob"]][2]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                      mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Mal3[["prob"]][3]*dGA(seq(min(Malic),max(Malic),length=length(Malic)),
                                      mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Mal3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
            AIC=c(fit.GA.Mal$aic,fit.GA.Mal3$aic),
            BIC=c(fit.GA.Mal$sbc,fit.GA.Mal3$sbc))

#ASH
summary(Ash)
frequencyAsh <- table(Ash)
frequencyAsh
length(Ash)
names(frequencyAsh)[frequencyAsh==max(frequencyAsh)]
labstatR::cv(Ash)
library(EnvStats)
kurtosis(Ash) #excess kurtosis
kurtosis(Ash, excess=FALSE)
skew(Ash)
boxplot(Ash, main = "Boxplot of Ash")
boxplot(Ash)$out
hist(Ash, breaks = 178, freq = FALSE)
lines(density(Ash), col = "red")

library(gamlss)
fit.EXP <- histDist(Ash, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Ash, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Ash, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Ash, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Ash, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Ash, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Ash,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Ash, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Ash <- gamlssMXfits(n = 5, Ash~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Ash[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Ash[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Ash[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Ash[["models"]][[2]][["sigma.coefficients"]])
hist(Ash, breaks = 178, freq = FALSE)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash[["prob"]][1]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                  mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash[["prob"]][2]*dGA(seq(min(Ash),
                                      max(Ash),length=length(Ash)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash[["prob"]][1]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                  mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Ash[["prob"]][2]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                     mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Ash)

library(gamlss)
set.seed(123)
fit.GA.Ash3 <- gamlssMXfits(n = 5, Ash~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Ash3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Ash3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Ash3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Ash3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Ash3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Ash3[["models"]][[3]][["sigma.coefficients"]])
hist(Ash, breaks = 178, freq = FALSE)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash3[["prob"]][1]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                   mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash3[["prob"]][2]*dGA(seq(min(Ash),
                                       max(Ash),length=length(Ash)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash3[["prob"]][2]*dGA(seq(min(Ash),
                                       max(Ash),length=length(Ash)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Ash),max(Ash),length=length(Ash)),
      fit.GA.Ash3[["prob"]][1]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                   mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Ash3[["prob"]][2]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                      mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Ash3[["prob"]][3]*dGA(seq(min(Ash),max(Ash),length=length(Ash)),
                                      mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Ash3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
           AIC=c(fit.GA.Ash$aic,fit.GA.Ash3$aic),
           BIC=c(fit.GA.Ash$sbc,fit.GA.Ash3$sbc))

#ALCALINITY
summary(Alcalinity)
frequencyAlcalinity <- table(Alcalinity)
frequencyAlcalinity
length(Alcalinity)
names(frequencyAlcalinity)[frequencyAlcalinity==max(frequencyAlcalinity)]
labstatR::cv(Alcalinity)
library(EnvStats)
kurtosis(Alcalinity) #excess kurtosis
kurtosis(Alcalinity, excess=FALSE)
skew(Alcalinity)
boxplot(Alcalinity, main = "Boxplot of Alcalinity")
boxplot(Alcalinity)$out
hist(Alcalinity, breaks = 178, freq = FALSE)
lines(density(Alcalinity), col = "red")

library(gamlss)
fit.EXP <- histDist(Alcalinity, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Alcalinity, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Alcalinity, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Alcalinity, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Alcalinity, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Alcalinity, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Alcalinity,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Alcalinity, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Alcalinity <- gamlssMXfits(n = 5, Alcalinity~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Alcalinity[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Alcalinity[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Alcalinity[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Alcalinity[["models"]][[2]][["sigma.coefficients"]])
hist(Alcalinity, breaks = 178, freq = FALSE)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity[["prob"]][1]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                         mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity[["prob"]][2]*dGA(seq(min(Alcalinity),
                                             max(Alcalinity),length=length(Alcalinity)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity[["prob"]][1]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                         mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Alcalinity[["prob"]][2]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                            mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Alcalinity)

library(gamlss)
set.seed(123)
fit.GA.Alcalinity3 <- gamlssMXfits(n = 5, Alcalinity~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Alcalinity3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Alcalinity3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Alcalinity3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Alcalinity3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Alcalinity3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Alcalinity3[["models"]][[3]][["sigma.coefficients"]])
hist(Alcalinity, breaks = 178, freq = FALSE)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity3[["prob"]][1]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                          mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity3[["prob"]][2]*dGA(seq(min(Alcalinity),
                                              max(Alcalinity),length=length(Alcalinity)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity3[["prob"]][2]*dGA(seq(min(Alcalinity),
                                              max(Alcalinity),length=length(Alcalinity)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
      fit.GA.Alcalinity3[["prob"]][1]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                          mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Alcalinity3[["prob"]][2]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                             mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Alcalinity3[["prob"]][3]*dGA(seq(min(Alcalinity),max(Alcalinity),length=length(Alcalinity)),
                                             mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Alcalinity3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
           AIC=c(fit.GA.Alcalinity$aic,fit.GA.Alcalinity3$aic),
           BIC=c(fit.GA.Alcalinity$sbc,fit.GA.Alcalinity3$sbc))

#MAGNESIUM
summary(Magnesium)
frequencyMagnesium <- table(Magnesium)
frequencyMagnesium
length(Magnesium)
names(frequencyMagnesium)[frequencyMagnesium==max(frequencyMagnesium)]
labstatR::cv(Magnesium)
library(EnvStats)
kurtosis(Magnesium) #excess kurtosis
kurtosis(Magnesium, excess=FALSE)
skew(Magnesium)
boxplot(Magnesium, main = "Boxplot of Magnesium")
boxplot(Magnesium)$out
hist(Magnesium, breaks = 178, freq = FALSE)
lines(density(Magnesium), col = "red")

library(gamlss)
fit.EXP <- histDist(Alcalinity, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Alcalinity, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Alcalinity, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Alcalinity, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Alcalinity, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Alcalinity, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Alcalinity,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Alcalinity, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Magnesium <- gamlssMXfits(n = 5, Magnesium~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Magnesium[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Magnesium[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Magnesium[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Magnesium[["models"]][[2]][["sigma.coefficients"]])
hist(Magnesium, breaks = 178, freq = FALSE)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium[["prob"]][1]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                        mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium[["prob"]][2]*dGA(seq(min(Magnesium),
                                            max(Magnesium),length=length(Magnesium)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium[["prob"]][1]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                        mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Magnesium[["prob"]][2]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                           mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Magnesium)

library(gamlss)
set.seed(123)
fit.GA.Magnesium3 <- gamlssMXfits(n = 5, Magnesium~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Magnesium3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Magnesium3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Magnesium3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Magnesium3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Magnesium3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Magnesium3[["models"]][[3]][["sigma.coefficients"]])
hist(Magnesium, breaks = 178, freq = FALSE)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium3[["prob"]][1]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                         mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium3[["prob"]][2]*dGA(seq(min(Magnesium),
                                             max(Magnesium),length=length(Magnesium)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium3[["prob"]][2]*dGA(seq(min(Magnesium),
                                             max(Magnesium),length=length(Magnesium)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
      fit.GA.Magnesium3[["prob"]][1]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                         mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Magnesium3[["prob"]][2]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                            mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Magnesium3[["prob"]][3]*dGA(seq(min(Magnesium),max(Magnesium),length=length(Magnesium)),
                                            mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Magnesium3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
           AIC=c(fit.GA.Magnesium$aic,fit.GA.Magnesium3$aic),
           BIC=c(fit.GA.Magnesium$sbc,fit.GA.Magnesium3$sbc))

#PHENOLS
summary(Phenols)
frequencyPhenols <- table(Phenols)
frequencyPhenols
length(Phenols)
names(frequencyPhenols)[frequencyPhenols==max(frequencyPhenols)]
labstatR::cv(Phenols)
library(EnvStats)
kurtosis(Phenols) #excess kurtosis
kurtosis(Phenols, excess=FALSE)
skew(Phenols)
boxplot(Phenols, main = "Boxplot of Phenols")
boxplot(Phenols)$out
hist(Phenols, breaks = 178, freq = FALSE)
lines(density(Phenols), col = "red")

library(gamlss)
fit.EXP <- histDist(Phenols, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Phenols, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Phenols, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Phenols, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Phenols, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Phenols, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Phenols,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Phenols, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Phenols <- gamlssMXfits(n = 5, Phenols~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Phenols[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Phenols[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Phenols[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Phenols[["models"]][[2]][["sigma.coefficients"]])
hist(Phenols, breaks = 178, freq = FALSE)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols[["prob"]][1]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                      mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols[["prob"]][2]*dGA(seq(min(Phenols),
                                          max(Phenols),length=length(Phenols)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols[["prob"]][1]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                      mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Phenols[["prob"]][2]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                         mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Phenols)

library(gamlss)
set.seed(123)
fit.GA.Phenols3 <- gamlssMXfits(n = 5, Phenols~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Phenols3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Phenols3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Phenols3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Phenols3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Phenols3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Phenols3[["models"]][[3]][["sigma.coefficients"]])
hist(Phenols, breaks = 178, freq = FALSE)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols3[["prob"]][1]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                       mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols3[["prob"]][2]*dGA(seq(min(Phenols),
                                           max(Phenols),length=length(Phenols)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols3[["prob"]][2]*dGA(seq(min(Phenols),
                                           max(Phenols),length=length(Phenols)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Phenols),max(Phenols),length=length(Phenols)),
      fit.GA.Phenols3[["prob"]][1]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                       mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Phenols3[["prob"]][2]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                          mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Phenols3[["prob"]][3]*dGA(seq(min(Phenols),max(Phenols),length=length(Phenols)),
                                          mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Phenols3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
           AIC=c(fit.GA.Phenols$aic,fit.GA.Phenols3$aic),
           BIC=c(fit.GA.Phenols$sbc,fit.GA.Phenols3$sbc))

#FLAVANOIDS
summary(Flavanoids)
frequencyFlavanoids <- table(Flavanoids)
frequencyFlavanoids
length(Flavanoids)
names(frequencyFlavanoids)[frequencyFlavanoids==max(frequencyFlavanoids)]
labstatR::cv(Flavanoids)
library(EnvStats)
kurtosis(Flavanoids) #excess kurtosis
kurtosis(Flavanoids, excess=FALSE)
skew(Flavanoids)
boxplot(Flavanoids, main = "Boxplot of Flavanoids")
boxplot(Flavanoids)$out
hist(Flavanoids, breaks = 178, freq = FALSE)
lines(density(Flavanoids), col = "red")

library(gamlss)
fit.EXP <- histDist(Flavanoids, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Flavanoids, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Flavanoids, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Flavanoids, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Flavanoids, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Flavanoids, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Flavanoids,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Flavanoids, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Flavanoids <- gamlssMXfits(n = 5, Flavanoids~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Flavanoids[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Flavanoids[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Flavanoids[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Flavanoids[["models"]][[2]][["sigma.coefficients"]])
hist(Flavanoids, breaks = 178, freq = FALSE)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids[["prob"]][1]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                         mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids[["prob"]][2]*dGA(seq(min(Flavanoids),
                                             max(Flavanoids),length=length(Flavanoids)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids[["prob"]][1]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                         mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Flavanoids[["prob"]][2]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                            mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Flavanoids)

library(gamlss)
set.seed(123)
fit.GA.Flavanoids3 <- gamlssMXfits(n = 5, Flavanoids~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Flavanoids3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Flavanoids3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Flavanoids3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Flavanoids3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Flavanoids3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Flavanoids3[["models"]][[3]][["sigma.coefficients"]])
hist(Flavanoids, breaks = 178, freq = FALSE)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids3[["prob"]][1]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                          mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids3[["prob"]][2]*dGA(seq(min(Flavanoids),
                                              max(Flavanoids),length=length(Flavanoids)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids3[["prob"]][2]*dGA(seq(min(Flavanoids),
                                              max(Flavanoids),length=length(Flavanoids)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
      fit.GA.Flavanoids3[["prob"]][1]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                          mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Flavanoids3[["prob"]][2]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                             mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Flavanoids3[["prob"]][3]*dGA(seq(min(Flavanoids),max(Flavanoids),length=length(Flavanoids)),
                                             mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Flavanoids3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
           AIC=c(fit.GA.Flavanoids$aic,fit.GA.Flavanoids3$aic),
           BIC=c(fit.GA.Flavanoids$sbc,fit.GA.Flavanoids3$sbc))

#NONFLAVANOID
summary(Nonflavanoid)
frequencyNonflavanoid <- table(Nonflavanoid)
frequencyNonflavanoid
length(Nonflavanoid)
names(frequencyNonflavanoid)[frequencyNonflavanoid==max(frequencyNonflavanoid)]
labstatR::cv(Nonflavanoid)
library(EnvStats)
kurtosis(Nonflavanoid) #excess kurtosis
kurtosis(Nonflavanoid, excess=FALSE)
skew(Nonflavanoid)
boxplot(Nonflavanoid, main = "Boxplot of Nonflavanoid")
boxplot(Nonflavanoid)$out
hist(Nonflavanoid, breaks = 178, freq = FALSE)
lines(density(Nonflavanoid), col = "red")

library(gamlss)
fit.EXP <- histDist(Nonflavanoid, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Nonflavanoid, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Nonflavanoid, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Nonflavanoid, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Nonflavanoid, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Nonflavanoid, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Nonflavanoid,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Nonflavanoid, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

library(gamlss)
set.seed(123)
fit.GA.Nonflavanoid <- gamlssMXfits(n = 5, Nonflavanoid~1, family = GA, K = 2)
mu.hat1 <- exp(fit.GA.Nonflavanoid[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Nonflavanoid[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Nonflavanoid[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Nonflavanoid[["models"]][[2]][["sigma.coefficients"]])
hist(Nonflavanoid, breaks = 178, freq = FALSE)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid[["prob"]][1]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                           mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid[["prob"]][2]*dGA(seq(min(Nonflavanoid),
                                               max(Nonflavanoid),length=length(Nonflavanoid)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid[["prob"]][1]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                           mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Nonflavanoid[["prob"]][2]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                              mu = mu.hat2, sigma = sigma.hat2),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Nonflavanoid)

library(gamlss)
set.seed(123)
fit.GA.Nonflavanoid3 <- gamlssMXfits(n = 5, Nonflavanoid~1, family = GA, K = 3)
mu.hat1 <- exp(fit.GA.Nonflavanoid3[["models"]][[1]][["mu.coefficients"]])  
sigma.hat1 <- exp(fit.GA.Nonflavanoid3[["models"]][[1]][["sigma.coefficients"]])
mu.hat2 <- exp(fit.GA.Nonflavanoid3[["models"]][[2]][["mu.coefficients"]])    
sigma.hat2 <- exp(fit.GA.Nonflavanoid3[["models"]][[2]][["sigma.coefficients"]])
mu.hat3 <- exp(fit.GA.Nonflavanoid3[["models"]][[3]][["mu.coefficients"]])  
sigma.hat3 <- exp(fit.GA.Nonflavanoid3[["models"]][[3]][["sigma.coefficients"]])
hist(Nonflavanoid, breaks = 178, freq = FALSE)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid3[["prob"]][1]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                            mu = mu.hat1, sigma = sigma.hat1),lty=2,lwd=3,col=2)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid3[["prob"]][2]*dGA(seq(min(Nonflavanoid),
                                                max(Nonflavanoid),length=length(Nonflavanoid)), mu = mu.hat2, sigma = sigma.hat2),
      lty=2,lwd=3,col=3)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid3[["prob"]][2]*dGA(seq(min(Nonflavanoid),
                                                max(Nonflavanoid),length=length(Nonflavanoid)), mu = mu.hat3, sigma = sigma.hat3),
      lty=2,lwd=3,col=4)
lines(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
      fit.GA.Nonflavanoid3[["prob"]][1]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                            mu = mu.hat1, sigma = sigma.hat1) +
         fit.GA.Nonflavanoid3[["prob"]][2]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                               mu = mu.hat2, sigma = sigma.hat2) +
         fit.GA.Nonflavanoid3[["prob"]][3]*dGA(seq(min(Nonflavanoid),max(Nonflavanoid),length=length(Nonflavanoid)),
                                               mu = mu.hat3, sigma = sigma.hat3),
      lty = 1, lwd = 3, col = 1)
print(fit.GA.Nonflavanoid3)

data.frame(row.names=c("Gamma mixture with K=2","Gamma mixture with K=3"),
           AIC=c(fit.GA.Nonflavanoid$aic,fit.GA.Nonflavanoid3$aic),
           BIC=c(fit.GA.Nonflavanoid$sbc,fit.GA.Nonflavanoid3$sbc))

#PROANTHOCYANINS
summary(Proanthocyanins)
frequencyProanthocyanins <- table(Proanthocyanins)
frequencyProanthocyanins
length(Proanthocyanins)
names(frequencyProanthocyanins)[frequencyProanthocyanins==max(frequencyProanthocyanins)]
labstatR::cv(Proanthocyanins)
library(EnvStats)
kurtosis(Proanthocyanins) #excess kurtosis
kurtosis(Proanthocyanins, excess=FALSE)
skew(Proanthocyanins)
boxplot(Proanthocyanins, main = "Boxplot of Proanthocyanins")
boxplot(Proanthocyanins)$out
hist(Proanthocyanins, breaks = 178, freq = FALSE)
lines(density(Proanthocyanins), col = "red")

library(gamlss)
fit.EXP <- histDist(Proanthocyanins, family=EXP, nbins = 178, main="Exponential distribution")
fit.GA <- histDist(Proanthocyanins, family=GA, nbins = 178, main="Gamma distribution")
fit.IG <- histDist(Proanthocyanins, family=IG, nbins = 178, main="Inverse Gaussian distribution")
fit.LOGNO <- histDist(Proanthocyanins, family=LOGNO, nbins = 178, main="Log-Normal distribution")
fit.WEI <- histDist(Proanthocyanins, family=WEI, nbins = 178, main="Weibull distribution")
fit.LO <- histDist(Proanthocyanins, family=LO, nbins=178, main="Logistic distribution")
fit.sep2 <-histDist(Proanthocyanins,family=SEP2, nbins=178, main="Skew Power Exp. Type 2 distribution")
fit.st5 <- histDist(Proanthocyanins, family=ST5, nbins=178, main = "Skew t type 5")

data.frame(row.names=c("Exponential", "Gamma", "Inverse Gaussian", "Log-Normal", 
                       "Weibull", "Logistic", "Skew Power 2", "Skew t 5"),
           AIC=c(AIC(fit.EXP),AIC(fit.GA), AIC(fit.IG), AIC(fit.LOGNO),
                 AIC(fit.WEI), AIC(fit.LO), AIC(fit.sep2), AIC(fit.st5)),
           SBC=c(fit.EXP$sbc, fit.GA$sbc, fit.IG$sbc, fit.LOGNO$sbc,
                 fit.WEI$sbc, fit.LO$sbc, fit.sep2$sbc, fit.st5$sbc))

#COLOR, GUE, DILUTION, PROLINE
length(Color)
hist(Color, breaks = 178, freq = FALSE)
lines(density(Color), col = "red")
length(Hue)
hist(Hue, breaks = 178, freq = FALSE)
lines(density(Hue), col = "red")
length(Dilution)
hist(Dilution, breaks = 178, freq = FALSE)
lines(density(Dilution), col = "red")
length(Proline)
hist(Proline, breaks = 178, freq = FALSE)
lines(density(Proline), col = "red")

#UNIVARIATE FOR EACH TYPE
typeBarbera = wine[Type == "Barbera",]
typeBarolo = wine[Type == "Barolo",]
typeGrignolino = wine[Type == "Grignolino",]

dim(typeBarbera)

length(typeBarbera$Alcohol)
length(typeBarbera$Malic)
length(typeBarbera$Ash)
length(typeBarbera$Alcalinity)
length(typeBarbera$Magnesium)
length(typeBarbera$Phenols)
length(typeBarbera$Flavanoids)
length(typeBarbera$Nonflavanoid)
length(typeBarbera$Proanthocyanins)
length(typeBarbera$Color)
length(typeBarbera$Hue)
length(typeBarbera$Dilution)
length(typeBarbera$Proline)

par(mfrow = c(2,2)) #multiple graphs in a single plot by
hist(typeBarbera$Alcohol, breaks = 48, freq=F)
lines(density(typeBarbera$Alcohol), col="red")
hist(typeBarbera$Malic, breaks = 48, freq=F)
lines(density(typeBarbera$Malic), col="red")
hist(typeBarbera$Ash, breaks = 48, freq=F)
lines(density(typeBarbera$Ash), col="red")
hist(typeBarbera$Alcalinity, breaks = 48, freq=F)
lines(density(typeBarbera$Alcalinity), col="red")
hist(typeBarbera$Magnesium, breaks = 48, freq=F)
lines(density(typeBarbera$Magnesium), col="red")
hist(typeBarbera$Phenols, breaks = 48, freq=F)
lines(density(typeBarbera$Phenols), col="red")
hist(typeBarbera$Flavanoids, breaks = 48, freq=F)
lines(density(typeBarbera$Flavanoids), col="red")
hist(typeBarbera$Nonflavanoid, breaks = 48, freq=F)
lines(density(typeBarbera$Nonflavanoid), col="red")
hist(typeBarbera$Proanthocyanins, breaks = 48, freq=F)
lines(density(typeBarbera$Proanthocyanins), col="red")
hist(typeBarbera$Color, breaks = 48, freq=F)
lines(density(typeBarbera$Color), col="red")
hist(typeBarbera$Hue, breaks = 48, freq=F)
lines(density(typeBarbera$Hue), col="red")
hist(typeBarbera$Dilution, breaks = 48, freq=F)
lines(density(typeBarbera$Dilution), col="red")
hist(typeBarbera$Proline, breaks = 48, freq=F)
lines(density(typeBarbera$Proline), col="red")

dim(typeBarolo)

length(typeBarolo$Alcohol)
length(typeBarolo$Malic)
length(typeBarolo$Ash)
length(typeBarolo$Alcalinity)
length(typeBarolo$Magnesium)
length(typeBarolo$Phenols)
length(typeBarolo$Flavanoids)
length(typeBarolo$Nonflavanoid)
length(typeBarolo$Proanthocyanins)
length(typeBarolo$Color)
length(typeBarolo$Hue)
length(typeBarolo$Dilution)
length(typeBarolo$Proline)

par(mfrow = c(2,2)) #multiple graphs in a single plot by
hist(typeBarolo$Alcohol, breaks = 59, freq=F)
lines(density(typeBarolo$Alcohol), col="red")
hist(typeBarolo$Malic, breaks = 59, freq=F)
lines(density(typeBarolo$Malic), col="red")
hist(typeBarolo$Ash, breaks = 59, freq=F)
lines(density(typeBarolo$Ash), col="red")
hist(typeBarolo$Alcalinity, breaks = 59, freq=F)
lines(density(typeBarolo$Alcalinity), col="red")
hist(typeBarolo$Magnesium, breaks = 59, freq=F)
lines(density(typeBarolo$Magnesium), col="red")
hist(typeBarolo$Phenols, breaks = 59, freq=F)
lines(density(typeBarolo$Phenols), col="red")
hist(typeBarolo$Flavanoids, breaks = 59, freq=F)
lines(density(typeBarolo$Flavanoids), col="red")
hist(typeBarolo$Nonflavanoid, breaks = 59, freq=F)
lines(density(typeBarolo$Nonflavanoid), col="red")
hist(typeBarolo$Proanthocyanins, breaks = 59, freq=F)
lines(density(typeBarolo$Proanthocyanins), col="red")
hist(typeBarolo$Color, breaks = 59, freq=F)
lines(density(typeBarolo$Color), col="red")
hist(typeBarolo$Hue, breaks = 59, freq=F)
lines(density(typeBarolo$Hue), col="red")
hist(typeBarolo$Dilution, breaks = 59, freq=F)
lines(density(typeBarolo$Dilution), col="red")
hist(typeBarolo$Proline, breaks = 59, freq=F)
lines(density(typeBarolo$Proline), col="red")

dim(typeGrignolino)
length(typeGrignolino$Alcohol)
length(typeGrignolino$Malic)
length(typeGrignolino$Ash)
length(typeGrignolino$Alcalinity)
length(typeGrignolino$Magnesium)
length(typeGrignolino$Phenols)
length(typeGrignolino$Flavanoids)
length(typeGrignolino$Nonflavanoid)
length(typeGrignolino$Proanthocyanins)
length(typeGrignolino$Color)
length(typeGrignolino$Hue)
length(typeGrignolino$Dilution)
length(typeGrignolino$Proline)

par(mfrow = c(2,2)) #multiple graphs in a single plot by
hist(typeGrignolino$Alcohol, breaks = 71, freq=F)
lines(density(typeGrignolino$Alcohol), col="red")
hist(typeGrignolino$Malic, breaks = 71, freq=F)
lines(density(typeGrignolino$Malic), col="red")
hist(typeGrignolino$Ash, breaks = 71, freq=F)
lines(density(typeGrignolino$Ash), col="red")
hist(typeGrignolino$Alcalinity, breaks = 71, freq=F)
lines(density(typeGrignolino$Alcalinity), col="red")
hist(typeGrignolino$Magnesium, breaks = 71, freq=F)
lines(density(typeGrignolino$Magnesium), col="red")
hist(typeGrignolino$Phenols, breaks = 71, freq=F)
lines(density(typeGrignolino$Phenols), col="red")
hist(typeGrignolino$Flavanoids, breaks = 71, freq=F)
lines(density(typeGrignolino$Flavanoids), col="red")
hist(typeGrignolino$Nonflavanoid, breaks = 71, freq=F)
lines(density(typeGrignolino$Nonflavanoid), col="red")
hist(typeGrignolino$Proanthocyanins, breaks = 71, freq=F)
lines(density(typeGrignolino$Proanthocyanins), col="red")
hist(typeGrignolino$Color, breaks = 71, freq=F)
lines(density(typeGrignolino$Color), col="red")
hist(typeGrignolino$Hue, breaks = 71, freq=F)
lines(density(typeGrignolino$Hue), col="red")
hist(typeGrignolino$Dilution, breaks = 71, freq=F)
lines(density(typeGrignolino$Dilution), col="red")
hist(typeGrignolino$Proline, breaks = 71, freq=F)
lines(density(typeGrignolino$Proline), col="red")

#PRINCIPAL COMPONENT ANALYSIS
cor(wine[, -1])

pr.out = prcomp(wine[2:14], scale=TRUE)
pr.out
names(pr.out)
pr.out$sdev
pr.out$rotation
pr.out$center
pr.out$scale
pr.out$x
head(pr.out$x)

pr.out$rotation = -pr.out$rotation
pr.out$rotation
pr.out$x = -pr.out$x
head(pr.out$x)

fviz_pca_ind(pr.out, title = "PCA - Wine", legend = "top",
             habillage = wine$Type, palette = c("green","blue","red"), 
             geom = "point", ggtheme = theme_classic())

biplot(pr.out, scale=0, cex = c(0.5, 1))

pr.var= pr.out$sdev^2
pr.var
which(pr.var>1)

PVE = pr.var/sum(pr.var)
PVE
CPVE=cumsum(PVE)
CPVE 

par(mfrow=c(1,2))
plot(PVE, xlab="Principal Component", ylab="Proportion of Variance Explained", ylim=c(0,1), type='b')
plot(cumsum(PVE), xlab="Principal Component", ylab="Cumulative Proportion of Variance Explained", ylim=c(0,1),type='b')

library(scatterplot3d)
scatterplot3d(pr.out$x[,1:3],pch=16,type="h", lty.hplot=3)

#MULTIVARIATE ANALYSIS: CLUSTERING (CLUSTER ANALYSIS)
#ASSESSING CLUSTER TENDENCY

df <- wine.scaled <- scale(wine[2:14]) #scale the data
pairs(df, gap=0, pch=16) #graphical representation

# Random data generated from the iris data set

random_df <- apply(df, 2,
                   function(x){runif(length(x), min(x), (max(x)))})
random_df <- as.data.frame(random_df)
scaled.random_df = scale(random_df)
pairs(scaled.random_df, gap = 0, pch = 21)

library("factoextra")

# Plot the standardized df data 

fviz_pca_ind(prcomp(df), title = "PCA - Wine dataset",
             habillage = wine$Type, palette = "jco",
             geom = "point", ggtheme = theme_classic(),
             legend = "bottom")

# Plot the random df data

fviz_pca_ind(prcomp(random_df), title = "PCA - Random data",
             geom = "point", ggtheme = theme_classic())

## Why assessing clustering tendency?

library(factoextra)
set.seed(123)

# K-means on df data

km.res1 <- kmeans(df, 3)

# plot on the original space

cl1 <- km.res1$cluster
pairs(df, gap=0, pch=cl1, col=c("black", "red", "orange")[cl1])

# K-means on random df data

km.res2 <- kmeans(random_df, 3)

# plot on the original space

cl2 <- km.res2$cluster
pairs(random_df, gap=0, pch=cl2, col=c("black", "red", "orange")[cl2])

# +++++++++++++++++++ #
# Clustering tendency #
# +++++++++++++++++++ #

## Hopkins statistic

library(clustertend)

# Compute Hopkins statistic for the iris dataset

set.seed(123)
hopkins(df, n = nrow(df)-1)

# Compute Hopkins statistic for a random dataset

set.seed(123)
hopkins(random_df, n = nrow(random_df)-1)

## VAT algorithm 

fviz_dist(dist(df), show_labels = FALSE)+
   labs(title = "Wine data")

fviz_dist(dist(random_df), show_labels = FALSE)+
   labs(title = "Random data")


res.dist <- dist(df, method = "euclidean")
res.hc <- hclust(d = res.dist, method = "single")
# Dendrogram
library("factoextra")
fviz_dend(res.hc, cex = 0.5, main = "Single linkage method and euclidean distance") # cex: label size
# Compute cophenetic distance
cor(res.dist, cophenetic(res.hc))

res.dist = dist(df, method = "euclidean")
res.hc <- hclust(d = res.dist, method = "average")
fviz_dend(res.hc, cex = 0.5, main = "Average linkage method and euclidean distance")
cor(res.dist, cophenetic(res.hc))

res.dist = dist(df, method = "euclidean")
res.hc <- hclust(d = res.dist, method = "complete")
fviz_dend(res.hc, cex = 0.5, main = "Complete linkage method and euclidean distance")
cor(res.dist, cophenetic(res.hc))

res.dist2 = dist(df, method = "manhattan")
res.hc2 <- hclust(d = res.dist2, method = "single")
fviz_dend(res.hc2, cex = 0.5, main = "Single linkage method and manhattan distance")
cor(res.dist2, cophenetic(res.hc2))

res.dist2 = dist(df, method = "manhattan")
res.hc2 <- hclust(d = res.dist2, method = "average")
fviz_dend(res.hc2, cex = 0.5, main = "Average linkage method and manhattan distance")
cor(res.dist2, cophenetic(res.hc2))

res.dist2 = dist(df, method = "manhattan")
res.hc2 <- hclust(d = res.dist2, method = "complete")
fviz_dend(res.hc2, cex = 0.5, main = "Complete linkage method and manhattan distance")
cor(res.dist2, cophenetic(res.hc2))

res.dist = dist(df, method = "euclidean")
res.centroid <- hclust(d = res.dist, method = "centroid")
fviz_dend(res.centroid, cex = 0.5, main = "Centroid linkage method and euclidean distance")
cor(res.dist, cophenetic(res.centroid))

res.dist = dist(df, method = "euclidean")
res.wardd <- hclust(d = res.dist, method = "ward.D")
fviz_dend(res.wardd, cex = 0.5, main = "Ward D linkage method and euclidean distance")
cor(res.dist, cophenetic(res.wardd))

hc.dist = dist(scaled.olive, method = "euclidean")
hc.wardd2 <- hclust(d = hc.dist, method = "ward.D2")
fviz_dend(hc.wardd2, cex = 0.5, main = "Ward D2 linkage method and euclidean distance")

res.dist = dist(df, method = "euclidean")
res.wardd2 <- hclust(d = res.dist, method = "ward.D2")
fviz_dend(res.wardd2, cex = 0.5, main = "Ward D2 linkage method and euclidean distance")
cor(res.dist, cophenetic(res.wardd2))

fviz_dend(hclust(daisy(wine), method="single"),cex = 0.5, main="Single linkage method and Gower distance")
fviz_dend(hclust(daisy(wine), method="average"), cex= 0.5, main="Average linkage method and Gower distance")
fviz_dend(hclust(daisy(wine), method="complete"),cex=0.5, main="Complete linkage method and Gower distance")
cor(cophenetic(hclust(daisy(wine),method="single")), daisy(wine))
cor(cophenetic(hclust(daisy(wine),method="average")), daisy(wine))
cor(cophenetic(hclust(daisy(wine),method="complete")), daisy(wine))

fviz_nbclust(df, hcut, method="silhouette", hc_metric="euclidean", hc_method= "single") 
fviz_nbclust(df, hcut, method="silhouette", hc_metric="manhattan", hc_method= "single")
fviz_nbclust(df, hcut, method="silhouette", hc_metric="euclidean", hc_method= "average")
fviz_nbclust(df, hcut, method="silhouette", hc_metric="manhattan", hc_method= "average")
fviz_nbclust(df, hcut, method="silhouette", hc_metric="euclidean", hc_method= "complete")
fviz_nbclust(df, hcut, method="silhouette", hc_metric="manhattan", hc_method= "complete")
fviz_nbclust(df, hcut, method="silhouette", hc_metric="euclidean", hc_method= "centroid")

hclust_SG=c()
for (i in 2:14){hclust_SG=c(hclust_SG,mean(silhouette(dist=daisy(wine),x=cutree(hclust (daisy(wine),method="single"),k=i))[,3]))}
hclust_CG=c()
for (i in 2:14){hclust_CG=c(hclust_CG,mean(silhouette(dist=daisy(wine),x=cutree(hclust (daisy(wine),method="complete"),k=i))[,3]))}
hclust_AG=c()
for (i in 2:14){hclust_AG=c(hclust_AG,mean(silhouette(dist=daisy(wine),x=cutree(hclust (daisy(wine),method="average"),k=i))[,3]))}
par(mfrow=c(2,1))
plot(x=c(2:14),y=hclust_SG,xlab="Number of clusters k",ylab="Average silhouette width", main="Single linkage method with Gower distance",type="b")+ abline(v=2,lty="dotted") 
plot(x=c(2:14),y=hclust_CG,xlab="Number of clusters k",ylab="Average silhouette width", main="Complete linkage method with Gower distance", type="b") + abline(v=2,lty="dotted")
plot(x=c(2:14),y=hclust_AG,xlab="Number of clusters k",ylab="Average silhouette width", main="Average linkage method with Gower distance",type="b")+ abline(v=2,lty="dotted")


fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="euclidean", hc_method= "single")
fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="manhattan", hc_method= "single")
fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="euclidean", hc_method= "average")
fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="manhattan", hc_method= "average")
fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="euclidean", hc_method= "complete")
fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="manhattan", hc_method= "complete")
fviz_nbclust(df, hcut, method = "gap_stat", nboot = 500, hc_metric="euclidean", hc_method= "centroid")


library("NbClust")
library("factoextra")
nb <- NbClust(df, distance = "euclidean", method = "single")
fviz_nbclust(nb)
nb <- NbClust(df, distance = "euclidean", method = "average")
fviz_nbclust(nb)
nb <- NbClust(df, distance = "euclidean", method = "complete")
fviz_nbclust(nb)
nb <- NbClust(df, distance = "manhattan", method = "single")
fviz_nbclust(nb)
nb <- NbClust(df, distance = "manhattan", method = "average")
fviz_nbclust(nb)
nb <- NbClust(df, distance = "manhattan", method = "complete")
fviz_nbclust(nb)
nb <- NbClust(df, distance = "euclidean", method = "centroid")
fviz_nbclust(nb)

set.seed(123)
fviz_nbclust(df,kmeans,method="silhouette",nstart=50)
fviz_nbclust(df,kmeans,method="gap_stat", nboot=500, nstart=50)

set.seed(123)
fviz_nbclust(NbClust(df,method="kmeans"))

set.seed(123)
fviz_nbclust(df,cluster::pam,method = "silhouette",metric= "euclidean")
fviz_nbclust(df,cluster::pam,method = "silhouette",metric= "manhattan")
fviz_nbclust(df,cluster::pam,method = "gap_stat", nboot=500, metric= "euclidean")
fviz_nbclust(df,cluster::pam,method = "gap_stat", nboot=500, metric= "manhattan")

library(clValid)
rownames(df) <- rownames(wine) # add rownames
clmethods <- c("hierarchical","kmeans","pam")

single_eucl <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "euclidean", method = "single")
summary(single_eucl)

single_man <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "manhattan", method = "single")
summary(single_man)

average_eucl <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "euclidean", method = "average")
summary(average_eucl)

average_man <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "manhattan", method = "average")
summary(average_man)

complete_eucl <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "euclidean", method = "complete")
summary(complete_eucl)

complete_man <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "manhattan", method = "complete")
summary(complete_man)

ward_eucl <- clValid(df, nClust = c(2,3), clMethods = clmethods, validation = c("internal", "stability"), metric = "euclidean", method = "ward")
summary(ward_eucl)

km.res <- eclust(df, "kmeans", k = 3, nstart = 25, graph = FALSE)
table(wine$Type, km.res$cluster)
library("fpc")
# Compute cluster stats
type <- as.numeric(wine$Type)
clust_stats <- cluster.stats(d = dist(df), type, km.res$cluster)
# Corrected Rand index
clust_stats$corrected.rand  
# Meila's VI index
clust_stats$vi

# Agreement between species and pam clusters
pam.res <- eclust(df, "pam", k = 3, graph = FALSE)
table(wine$Type, pam.res$cluster)
cluster.stats(d = dist(df), type, pam.res$cluster)$vi
# Agreement between species and HC clusters
res.hc <- eclust(df, "hclust", k = 3, graph = FALSE)
table(wine$Type, res.hc$cluster)
cluster.stats(d = dist(df), type, res.hc$cluster)$vi

library(gridExtra)
grid.arrange(fviz_dend(hclust(dist(df,method="euclidean"),method="single"), cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Single linkage method and euclidean distance"), 
             fviz_dend(hclust(dist(df,method="manhattan"),method="single"), cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Single linkage method and manhattan distance"),
             fviz_dend(hclust(dist(df,method="euclidean"),method="complete"),cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Complete linkage method and euclidean distance"), 
             fviz_dend(hclust(dist(df,method="manhattan"),method="complete"),cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Complete linkage method and manhattan distance"), 
             fviz_dend(hclust(dist(df,method="euclidean"),method="average"),cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T, main="Average linkage method and euclidean distance"), 
             fviz_dend(hclust(dist(df,method="manhattan"),method="average"),cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Average linkage method and manhattan distance"),
             ncol=2,nrow=3)

fviz_dend(hclust(dist(df,method="euclidean"),method="ward.D2"),cex=0.3,k=2, k_colors="jco",color_labels_by_k = T,rect=T,main="Ward D.2 linkage method")
fviz_dend(hclust(dist(df,method="euclidean"),method="ward.D"),cex=0.3,k=2, k_colors="jco",color_labels_by_k = T,rect=T,main="Ward D linkage method")

fviz_dend(hclust(daisy(wine),method="single"),cex=0.3,k=2,k_colors= "jco",color_labels_by_k = T,rect=T,main="Single linkage method and Gower distance", lower_rect=0)
fviz_dend(hclust(daisy(wine),method="complete"),cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Complete linkage method and Gower distance",lower_rect=0)
fviz_dend(hclust(daisy(wine),method="average"),cex=0.3,k=2,k_colors="jco",color_labels_by_k = T,rect=T,main="Average linkage method and Gower distance",lower_rect=0)


grid.arrange(fviz_dend(hclust(dist(df,method="euclidean"),method="single"), cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Single linkage method and euclidean distance"), 
             fviz_dend(hclust(dist(df,method="manhattan"),method="single"), cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Single linkage method and manhattan distance"),
             fviz_dend(hclust(dist(df,method="euclidean"),method="complete"),cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Complete linkage method and euclidean distance"), 
             fviz_dend(hclust(dist(df,method="manhattan"),method="complete"),cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Complete linkage method and manhattan distance"), 
             fviz_dend(hclust(dist(df,method="euclidean"),method="average"),cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T, main="Average linkage method and euclidean distance"), 
             fviz_dend(hclust(dist(df,method="manhattan"),method="average"),cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Average linkage method and manhattan distance"),
             ncol=2,nrow=3)

fviz_dend(hclust(dist(df,method="euclidean"),method="ward.D2"),cex=0.3,k=3, k_colors="jco",color_labels_by_k = T,rect=T,main="Ward D.2 linkage method")
fviz_dend(hclust(dist(df,method="euclidean"),method="ward.D"),cex=0.3,k=3, k_colors="jco",color_labels_by_k = T,rect=T,main="Ward D linkage method")

fviz_dend(hclust(daisy(wine),method="single"),cex=0.3,k=3,k_colors= "jco",color_labels_by_k = T,rect=T,main="Single linkage method and Gower distance", lower_rect=0)
fviz_dend(hclust(daisy(wine),method="complete"),cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Complete linkage method and Gower distance",lower_rect=0)
fviz_dend(hclust(daisy(wine),method="average"),cex=0.3,k=3,k_colors="jco",color_labels_by_k = T,rect=T,main="Average linkage method and Gower distance",lower_rect=0)

library(mclust)
library(gclus) # for the wine data

data(wine, package = "gclus")
head(wine)
df <- scale(wine[2:14])
mod <- Mclust(df)
fviz_mclust(mod, "BIC",palette = "jco")
summary(mod)
summary(mod$BIC)
plot(mod, what = "BIC", ylim = range(mod$BIC, na.rm = TRUE), legendArgs = list(x = "bottomleft"))
mod$modelName                # Selected parsimonious configuration
mod$G                        # Optimal number of clusters
head(round(mod$z, 6), 30)    # Probality to belong to a given cluster
head(mod$classification, 30) # Cluster assignement of each observation

# pairs plot with classification

pairs(df, gap=0, pch = 16, col = mod$classification)

table(Type, mod$classification)

# Adjusted (or Corrected) Rand Index

adjustedRandIndex(Type, mod$classification)

# Visulaing the obtained results in the PC space

library(factoextra)

# BIC values used for choosing the number of clusters

fviz_mclust(mod, "BIC", palette = "jco")

# Classification: plot showing the clustering

fviz_mclust(mod, "classification", geom = "point", pointsize = 1.5, palette = "jco")

# Classification uncertainty
# larger symbols indicate the more uncertain observations

fviz_mclust(mod, "uncertainty", palette = "jco")

# Plot the data using only two variables of interest, let say here "Alcohol" and "Magnesium"

fviz_mclust(mod, "uncertainty", palette = "jco", choose.vars = c("Alcohol", "Magnesium"))

