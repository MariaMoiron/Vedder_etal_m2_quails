# Code for "Rapid decline of prenatal maternal effects with age is independent of postnatal environment in a precocial bird"
# Published in Evolution, 2023, doi: tba
# Vedder O, Tschirren B, Postma E, Moiron M

# The code provided here is sufficient to replicate the analyses presented in the above paper

######################################################
# DATA ANALYSIS OF BODY MASS GROWTH DATA (Table S2)
######################################################

# Loading packages
library(nadiv)
library(MCMCglmm)
library(tidyr)
library(dplyr)
library(ggplot2)

# Loading phenotypic data
Data <- read.table("data.txt", header=T)

#Response variable
Data$mass = as.numeric(Data$mass)
hist(Data$mass)

#Fixed effects
Data$diet=as.factor(Data$chick.diet)
Data$sex= as.factor(Data$chick.sex)
Data$year = as.factor(Data$year)
Data$f = as.numeric(Data$f)
Data$motherL=as.factor(Data$mother.type)
Data$fatherL=as.factor(Data$father.type)
Data$motherR=as.factor(Data$mother.replicate)
Data$fatherR=as.factor(Data$father.replicate)
Data$ageC=as.factor(Data$age)

#Random effects
Data$ID = as.factor(Data$Offspring.ID)
Data$animal = as.factor(Data$Offspring.ID)
Data$pair=as.factor(paste(Data$father.ID, Data$mother.ID))
Data$mother.ID = as.factor(Data$mother.ID)
Data$father.ID=as.factor(Data$father.ID)

#Load pedigree info
ped<- read.table("quail.ped.txt",header=TRUE)
colnames(ped)[1] <- "animal"
ped3=prunePed(ped, Data$animal, make.base=TRUE)
str(ped3)
my_inverse <- inverseA(ped3)$Ainv

# Setting number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 10000
(NITT <- BURN + THIN*nsamp)

# Setting prior
prior <- list(R = list(V = diag(1)*100, nu = 1.002),
                      G = list(G1 = list(V = diag(1)*100, nu = 1.002),
                               G2 = list(V = diag(1)*100, nu = 1.002),
                               G3 =  list(V = diag(1)*5, nu = 1.002),
                               G4 =  list(V = diag(1)*5, nu = 1.002),
                               G5 =  list(V = diag(1)*5, nu = 1.002)))

# Running model
mod<- MCMCglmm(mass~ageC*sex+f+ageC*diet+year+motherL*fatherL*diet+motherR+fatherR,
                   random = ~ ID +animal+pair+mother.ID+father.ID,
                   rcov = ~ units,
                   data = Data,
                   prior =prior,
                   ginverse=list(animal = my_inverse),
                   family = "gaussian",
                   nitt = NITT, thin = THIN, burnin = BURN,
                   #pr=TRUE
)

summary(mod)

# Assesing convergence 
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

#Fixed effects
posterior.mode(mod$Sol)
HPDinterval(mod$Sol)
#plot(mod$Sol)

#Random effects
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV), 3)
plot(mod$VCV)

h2=mod$VCV[,"animal"]/rowSums(mod$VCV)  #different way of estimating Rpt
mean(h2)###mean
HPDinterval(h2) ###Calculate 95%CI
plot(h2)

rpt=(mod$VCV[,"animal"]+mod$VCV[,"ID"])/rowSums(mod$VCV)  #different way of estimating Rpt
mean(rpt)###mean
HPDinterval(rpt) ###Calculate 95%CI
plot(rpt)
