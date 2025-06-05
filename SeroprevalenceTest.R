# script to estimate power to compare seroprevalence and force of infection between predicted high and low risk areas

# assume a time constant and age constant FOI


start.time <- Sys.time()

# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)
library(binom)
library(plyr)
library(ggplot2)

#*******Example simulate serology data*******************

#************seroprevalence as a function of FOI and age*****************
propImmFunc <- function(lambda, age){
return(1 - exp(-lambda*age))
}
propImmFuncV <- Vectorize(propImmFunc)
#********************************************************************

#******************'true' FOI and Study design parameters*************
lowLambda <- 0.02
highLambda <- 0.04
sdLogFOI <- 0.2
# sampling / design choices
n.village <- 18    # number of villages
n.hh <- 10         # number of households per village
people.in.household <- 2

dat <- expand.grid(hh = 1:n.hh, village = 1:n.village,people=1:people.in.household)
# allocate villages to high and low prevalence in 1:1 ratio 
dat$risk.level <- dat$village %% 2 - 0.5
dat$FOI <- NA
dat$FOI[dat$risk.level==-0.5] <- lowLambda
dat$FOI[dat$risk.level==0.5] <- highLambda

# add village-level FOI variation
dat$villageFOI <- exp(rnorm(length(dat$village),mean=log(dat$FOI),sd=sdLogFOI))

hist(dat$villageFOI)

# randomly assign an age  - simplest option uniform
dat$age <- NA
dat$age<-sapply(1:length(dat$age),function(x){round(runif(1,1,70),0)})
# calculate true prob immune
dat$probImm <- propImmFuncV(lambda=dat$villageFOI,age=dat$age) #

# from 'true' proportion immune to random variation                                    
dat$Infected <- sapply(1:length(dat$probImm),function(x){
rbinom(n=1,size=1,prob=dat$probImm[x])
})

# summarise
summary <- ddply(dat,.(risk.level),summarise,num=sum(Infected),denom=length(dat$hh))
meanPrev <- binom.confint(summary$num,summary$denom,methods="exact")
meanPrev$risk.level <- c("Low","High")

ggplot(meanPrev) +
  geom_point(aes(x=risk.level,y=mean*100)) +
  geom_errorbar(aes(x=risk.level,ymin=lower*100,ymax=upper*100))

# ****put into function - SeroprevalenceFOISimFunc.R
