# Code for "Rapid decline of prenatal maternal effects with age is independent of postnatal environment in a precocial bird"
# Published in Evolution, 2023, doi: tba
# Vedder O, Tschirren B, Postma E, Moiron M

# The code provided here is sufficient to replicate the analyses presented in the above paper

################################################################################################
# DATA ANALYSES OF HYBRID APPROACH TO TEST FOR MATERNAL BY ENVIRONMENT INTERACTIONS (Table S4)
################################################################################################

# Load packages
library(nadiv)
library(MCMCglmm)
library(tidyr)
library(dplyr)
library(ggplot2)

# Load phenotypic data
Dataset <- read.table("maternal_by_environemnt_interactions.txt", header=T)

# Fixed effects
Data$diet=as.factor(Data$chick.diet)
Data$sex= as.factor(Data$chick.sex)
Data$year = as.factor(Data$year)
Data$f = as.numeric(Data$f)
Data$motherL=as.factor(Data$mother.type)
Data$fatherL=as.factor(Data$father.type)
Data$motherR=as.factor(Data$mother.replicate)
Data$fatherR=as.factor(Data$father.replicate)

# Random effects
Data$ID = as.factor(Data$offspring.ID)
Data$animal = as.factor(Data$offspring.ID)
Data$dam = as.factor(Data$mother.ID)
Data$mother = as.factor(Data$mother.ID)
Data$father=as.factor(Data$father.ID)

# Response variable
Data$body.mass=as.numeric(Data$mass.d0)
#Data$body.mass=as.numeric(Data$mass.d7)
hist(Data$body.mass)

# within-subject centering approach
# to estimate average line mass
df_scaled <- as.data.frame(Data %>%group_by(Data$motherL) %>% mutate(mean.mass.Line = as.numeric(mean(egg.mass))))
Data=df_scaled

# to estimate average replicate mass
df_scaled <- as.data.frame(Data %>% group_by(motherR, motherL) %>% mutate(mean.mass.R = as.numeric(mean(egg.mass))))
Data=df_scaled

# to estimate average mother mass
df_scaled <- as.data.frame(Data %>%group_by(motherR, motherL,mother) %>% mutate(mean.mass.Mother = as.numeric(mean(egg.mass))))
Data=df_scaled

# to estimate delta (average replicate - average line)
Data$sd.mass.R = as.numeric(Data$mean.mass.R-Data$mean.mass.Line)

# to estimate delta (average mother - average replicate) 
Data$sd.mass.Mother = as.numeric(Data$mean.mass.Mother-Data$mean.mass.R)

# to estimate delta (observation - average mother)
Data$sd.mass.within.mother = as.numeric(Data$egg.mass-Data$mean.mass.Mother)

# Load pedigree info
ped<- read.table("pedigree.txt",header=TRUE)
colnames(ped)[1] <- "animal"
ped3=prunePed(ped, Data$animal, make.base=TRUE)
my_inverse <- inverseA(ped3)$Ainv

# Set number of samples and iterations
nsamp <- 1000 
BURN <- 10000; THIN <- 20000
(NITT <- BURN + THIN*nsamp)

# Set prior
prior <- list(R = list(V = diag(2)*0.05, nu = 3),
              G = list(G1 = list(V = diag(2)*0.05, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                       G2 = list(V = diag(2)*0.05, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                       G3 = list(V = diag(2)*0.05, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000),
                       G4 = list(V = diag(2)*0.05, nu = 3, alpha.mu = c(0,0), alpha.V = diag(2)*1000)))

# Run model
mod <- MCMCglmm(mass~sex+year+f+motherL*fatherL+motherR+fatherR,
                random = ~ us(diet):animal + us(diet):dam+us(diet):mother+ us(diet):father,
                rcov = ~ idh(diet):units,
                ginverse=list(animal=my_inverse, dam=my_inverse),
                data = Data,
                prior =prior,
                family = "gaussian",
                nitt = NITT, thin = THIN, burnin = BURN#,
                #pr=TRUE
)

summary(mod)

# Assess model convergence 
effectiveSize(mod$VCV)
heidel.diag(mod$VCV)
autocorr.diag(mod$VCV)

# Fixed effects
posterior.mode(mod$Sol)
HPDinterval(mod$Sol)
#plot(mod$Sol)

# Random effects
round(posterior.mode(mod$VCV),3)
round(HPDinterval(mod$VCV),3)
plot(mod$VCV)
