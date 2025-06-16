library(stringr)
library(serofoi)
# set up data frame with varying FOIs

myTryCatch <- function(expr) {
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error=function(e) {
      err <<- e
      NULL
    }), warning=function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value=value
       #, warning=warn, error=err
  )
}



simulateTVFOI <- function(foiByYears=yf
                            ,ageMin=1
                            ,ageMax=55
                            ,sampleSizeByAge=10 # or a vector of sample length as ageMin:ageMax
                            
){

  
    ages <- c(ageMin:ageMax)
    birthYears <- max(yf$years) - ages
    yf$probSus <- exp(-yf$fois)
    
    probExp <- numeric(length(birthYears))
    for(i in 1:length(birthYears)){probExp[i] <- 1-prod(yf$probSus[yf$years>birthYears[i]])}

   ageInf <- cbind.data.frame(ages,probExp)
   ageInf$seropos <- sapply(1:length(ageInf$probExp),function(x){
     sum((rbinom(n=sampleSizeByAge,size=1,prob=ageInf$probExp[x])))
   })
    ageInf$n <- sampleSizeByAge


  labs <- c(paste(seq(0, ageMax, by = 3), seq(0 + 3 - 1, ageMax+3 - 1, by = 3),
                sep = "-"), paste(ageMax+3, "+", sep = ""))


  ageInf$AgeGroup <- cut(ageInf$ages, breaks = c(seq(0, ageMax+3, by = 3),Inf), labels = labs, right = FALSE)

  simsero <- ddply(ageInf,.(AgeGroup),summarise
                  ,n_seropositive=sum(seropos)
                  ,n_sample=sum(n))
  ages <- do.call(rbind,(str_split(simsero$AgeGroup,pattern="-")))
  simsero$age_min <- as.numeric(ages[,1])
  simsero$age_max <- as.numeric(ages[,2])
  simsero <- simsero[,-1]
  simsero$survey_year <- max(yf$years)
  
  foii <- get_foi_index(simsero,group_size=5,model_type="time")
  test <- myTryCatch(fit_seromodel(
    serosurvey = simsero,
    model_type = "time",
    foi_index = data.frame(
      year = min(foii$year):max(foii$year),
      foi_index = foii[,2]
    ),
    iter = 3000 )
  )
  modOutput <- extract_central_estimates(test$value)
  
  
  labs <- c(paste(seq(min(yf$years), max(yf$years), by = 5), seq(min(yf$years) + 5 - 1, max(yf$years), by = 5),
                  sep = "-"))
  
  
  yf$yrGrps<- cut(yf$years,breaks = c(seq(min(yf$years), max(yf$years)+1, by = 5)), labels = labs, right = FALSE)
  trueFOISummary <- ddply(yf,.(yrGrps),summarise,meanFOI=mean(fois))
  
  modOutput$trueVal <- trueFOISummary[,2]
  modOutput$yearGrp <- trueFOISummary[,1]
  
  
  acc <- sum(modOutput$trueVal-modOutput$median)^2
  ciRange <- mean(modOutput$upper - modOutput$lower)
  nsamples <- sum(simsero$n_sample)
  return(c(acc,ciRange,nsamples))

}



#*********************test out function********************************
years <- seq(1969,2023,1)
fois <- rep(0.005,length(years))
fois[33:48] <- 0.01
fois[49:length(fois)] <- 0.04

yf <- cbind.data.frame(years,fois)
ami <- 1
ama <- 55

simulateTVFOI(foiByYears=yf
              ,ageMin=ami
              ,ageMax=ama
              ,sampleSizeByAge=10)



#****************************************************************
library(parallel)

cl <- makeCluster(detectCores()-1)               # get cores from your computer - 1
clusterExport(cl, varlist=c("simulateTVFOI"
                            ,"yf"
                            ,"myTryCatch"
                            ,"get_foi_index"
                            ,"fit_seromodel"
                            ),
              envir=environment())     

clusterEvalQ(cl, {
  library(stringr)
  library(serofoi)
  library(plyr)   
})


start <- Sys.time()
sim.res <- parLapply(cl,1:100, 
function(i){
  simulateTVFOI(foiByYears=yf,
                ageMin=1,
                ageMax=55,
                sampleSizeByAge=35)
})
stopCluster(cl)   # stop the clusters
end <- Sys.time()
end - start

ss40 <- do.call(rbind.data.frame,sim.res)
ss40Means <- colMeans(ss40)
names(ss40Means) <- c("acc","uncert","sampleSize")

ss35 <- do.call(rbind.data.frame,sim.res)
ss35Means <- colMeans(ss35)
names(ss35Means) <- c("acc","uncert","sampleSize")

ss30 <- do.call(rbind.data.frame,sim.res)
ss30Means <- colMeans(ss30)
names(ss30Means) <- c("acc","uncert","sampleSize")

ss25 <- do.call(rbind.data.frame,sim.res)
ss25Means <- colMeans(ss25)
names(ss25Means) <- c("acc","uncert","sampleSize")

ss20 <- do.call(rbind.data.frame,sim.res)
ss20Means <- colMeans(ss20)
names(ss20Means) <- c("acc","uncert","sampleSize")

ss15 <- do.call(rbind.data.frame,sim.res)
ss15Means <- colMeans(ss15)
names(ss15Means) <- c("acc","uncert","sampleSize")

ss10 <- do.call(rbind.data.frame,sim.res)
ss10Means <- colMeans(ss10)
names(ss10Means) <- c("acc","uncert","sampleSize")

ss5 <- do.call(rbind.data.frame,sim.res)
ss5Means <- colMeans(ss5)
names(ss5Means) <- c("acc","uncert","sampleSize")

compareRes <- rbind.data.frame(ss40Means,ss35Means,ss30Means,ss25Means,ss20Means,ss15Means,ss10Means,ss5Means)
names(compareRes) <- c("acc","uncert","sampleSize")

saveRDS(compareRes,"compareRes.rds")

ggplot(compareRes) +
  geom_point(aes(x=sampleSize,y=acc))


ggplot(compareRes) +
  geom_point(aes(x=sampleSize,y=uncert))

