
# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)
library(binom)
library(plyr)
library(ggplot2)
library(lmtest)

source("SeroprevalenceFOISimFunc.R")

#*************Simulate, model fit and calculate p-value************
# function to simulate data and estimate p-value for null hypothesis 
# that high and low risk areas have the same seroprevalence
res.tab.fn <- function(...) {
  # create template data set
  dat <- simulateSeroprevalence(lowLambda=0.01
                         ,highLambda=0.02
                         ,n.village=18
                         ,n.hh=15
                         ,people.in.household=2
                         ,ageMin=1
                         ,ageMax=70
                         ,ageMean=10
                         ,sdLogFOI=0.1
                         ,ageSD=2
                         ,distAge="uniform")
  
  fit <- glm(Infected~risk.level
      ,family=binomial(link="cloglog")
      ,offset=log(age)
      ,data=dat)
  
  fit0 <- update(fit, ~ . - risk.level)
  pval <- lrtest(fit, fit0)[2,5]
  coefs <- coef(fit)
  return(c(pval,coefs))
}


# no of data sets to simulate (1000 takes ~ 4 min)
nsim <- 1000


# repeat simulations many times and calculate p-value
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterExport(cl, varlist=c("res.tab.fn","simulateSeroprevalence","propImmFuncV","lrtest"),
              envir=environment())               


start <- Sys.time()
sim.res <- parLapply(cl,1:nsim,res.tab.fn) 
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


# estimate power
pvals <- c(sapply(sim.res,"[[",1))
length(pvals[pvals<=0.05])

risk.level <- c(sapply(sim.res,"[[",3))
plot(risk.level)
var(risk.level)
range(risk.level)
# 360  - 0.15 - 1.3
