
# load packages
library(GLMMmisc) # available via devtools::install_github("pcdjohnson/GLMMmisc")
library(lme4)
library(parallel)
library(binom)
library(plyr)
library(ggplot2)
library(lmtest)
library(grid)
library(gridExtra)

source("SeroprevalenceFOISimFunc.R")


# accuracy and precision
accPrecFunc <- function(risk.level=risk,llambda=0.01,hlambda=0.02){
  bias <- mean(risk.level-(llambda/hlambda))
  variance <- mean((risk.level-mean(risk.level))^2)
  mse <- bias^2+variance
  return(mse)
}



#*************Simulate, model fit ************

numHH <- c(5:30)

# repeat simulations many times and calculate p-value
cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterExport(cl, varlist=c("simulateSeroprevalence","propImmFuncV","accPrecFunc"),
              envir=environment())               


  start <- Sys.time()
  sim.res <- parLapply(cl,numHH,function(x){
    nsim <- 1
    risk <- numeric(nsim)
    pval <- numeric(nsim)
    while(nsim<1000){
      dat <- simulateSeroprevalence(lowLambda=0.01
                                    ,highLambda=0.02
                                    ,n.village=18
                                    ,n.hh=x
                                    ,people.in.household=1
                                    ,ageMin=5
                                    ,ageMax=60
                                    ,ageMean=20
                                    ,sdLogFOI=0.5
                                    ,ageSD=2
                                    ,distAge="uniform")
    
      fit <- glm(Infected~risk.level
                 ,family=binomial(link="cloglog")
                 ,offset=log(age)
                 ,data=dat)
     pvals <- coef(summary(fit))[,4]
     pval[nsim] <- pvals[2]
     risk[nsim] <- coef(fit)[2]
     nsim <- nsim+1
    }
    mse <- accPrecFunc(risk.level=risk)
    ss <- x*18
    pwr <- length(pval[pval<=0.05])/nsim
  return(c(ss,mse,pwr))
  }) 
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start


sims <- do.call(rbind.data.frame,sim.res)
names(sims) <- c("samplesize","MSE","pwr")

msePlot <- ggplot(sims) +
  geom_point(aes(x=samplesize,y=MSE)) +
  xlab("Sample size") +
  ylab("Mean squared error") +
  theme(text = element_text(size = 20))

pwerPlot <- ggplot(sims) +
  geom_point(aes(x=samplesize,y=pwr*100)) +
  xlab("Sample size") +
  ylab("Power (%)") +
  theme(text = element_text(size = 20))

grid.arrange(msePlot,pwerPlot,ncol=2)
