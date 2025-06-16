library(stringr)
# set up data frame with varying FOIs

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

  return(simsero)
}

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
  list(value=value, warning=warn, error=err)
}

#*********************test out function********************************
years <- seq(1969,2023,1)
fois <- rep(0.005,length(years))
fois[33:48] <- 0.01
fois[49:length(fois)] <- 0.04

yf <- cbind.data.frame(years,fois)
ami <- 1
ama <- 55

sims <- simulateTVFOI(foiByYears=yf
                      ,ageMin=ami
                      ,ageMax=ama
                      ,sampleSizeByAge=10 # or a vector of sample length as ageMin:ageMax
                      
)

foii <- get_foi_index(sims,group_size=5,model_type="time")

seromodel <- fit_seromodel(
  serosurvey = sims,
  model_type = "time",
  foi_index = data.frame(
    year = min(foii$year):max(foii$year),
    foi_index = foii[,2]
  ),
  iter = 3000
)
plot_seromodel(seromodel=seromodel,serosurvey=sims)
plot_foi_estimates(seromodel=seromodel,serosurvey=sims)

#**********************place within try catch**********************


test <- myTryCatch(fit_seromodel(
  serosurvey = sims,
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

sum(modOutput$trueVal-modOutput$median)^2

