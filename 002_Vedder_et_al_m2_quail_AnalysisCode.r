# Code for "Rapid decline of prenatal maternal effects with age is independent of postnatal environment in a precocial bird"
# Unpublished manuscript, doi: tba
# Vedder O, Tschirren B, Postma E, Moiron M

# The code provided here is sufficient to replicate the simulations presented in the above paper


######################################################
# DATA ANALYSIS OF BODY MASS DATA (Table S3)
######################################################

# Loading packages
library(nadiv)
library(MCMCglmm)
library(tidyr)
library(dplyr)
library(ggplot2)

# Loading phenotypic data
Dataset <- read.table("data.txt", header=T)

#Response variable
Data = as.data.frame(Dataset %>% group_by(diet) %>% mutate(mass = as.numeric(scale(mass, scale = FALSE))))
Data$mass = as.numeric(Data$mass)
hist(Data$mass)

#Fixed effects
Data$diet=as.factor(Data$chick.diet)
Data$sex= as.factor(Data$chick.sex)
Data$year = as.factor(Data$year)
Data$f = as.numeric(Data$f)
Data$motherL=as.factor(Data$mother.type)
Data$fatherL=as.factor(Data$father.type)
Data$egg.mass=as.numeric(Data$egg.mass)
Data$motherR=as.factor(Data$mother.replicate)
Data$fatherR=as.factor(Data$father.replicate)

#Random effects
Data$ID = as.factor(Data$Offspring.ID)
Data$animal = as.factor(Data$Offspring.ID)
Data$mother = as.factor(Data$mother.ID)
Data$father=as.factor(Data$fatherID)

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
prior1 <- list(R = list(V = diag(2)*5, nu = 2),
                      G = list(G1 = list(V = diag(2)*5, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                               G2 = list(V = diag(2)*5, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                               G3 = list(V = diag(2)*5, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                               G4 = list(V = diag(2)*5, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))

prior2 <- list(R = list(V = diag(2)*100, nu = 2),
               G = list(G1 = list(V = diag(2)*100, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                        G2 = list(V = diag(2)*100, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                        G3 = list(V = diag(2)*100, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                        G4 = list(V = diag(2)*100, nu = 2, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))

mod <- MCMCglmm(mass~sex+year+f+fatherL+motherR+fatherR,
                   random = ~ us(diet):animal + us(diet):dam+us(diet):mother+ us(diet):father,
                   rcov = ~ idh(diet):units,
                   ginverse=list(animal=my_inverse, dam=my_inverse),
                   data = Data,
                   prior =prior1, #change prior according to data
                   family = "gaussian",
                   nitt = NITT, thin = THIN, burnin = BURN#,
                   #pr=TRUE
)

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
