
library(serofoi)
library(plyr)

#***********************Example from serofoi package***************
tc <- fit_seromodel(serosurvey=chagas2012)

seromodel <- fit_seromodel(
  serosurvey = chagas2012
)

plot_seromodel(seromodel=seromodel,serosurvey=chagas2012)


seromodel <- fit_seromodel(
  serosurvey = chagas2012,
  model_type = "time",
  foi_index = data.frame(
    year = 1935:2011,
    foi_index = c(rep(1, 46), rep(2, 31))
  ),
  iter = 1000
)

plot_seromodel(seromodel=seromodel,serosurvey=chagas2012)
plot_foi_estimates(seromodel=seromodel,serosurvey=chagas2012)



#*******************CHIKV data from Lim et al 2023*********

n_sample <- c(61,143,144,124,114,116,73,64,53,48,59)
survey_year <- c(2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014)
n_seropositive <- c(1,4,7,16,36,49,34,38,33,30,43)
age_min <- c(1,5,10,15,20,25,30,35,40,45,50)
age_max <- c(4,9,14,19,24,29,34,39,44,49,55)

chikv<-cbind.data.frame(n_sample,survey_year,n_seropositive,age_min,age_max)

chikvfoii <- get_foi_index(chikv,group_size=10,model_type="time")

seromodel <- fit_seromodel(
  serosurvey = chikv,
  model_type = "time",
  foi_index = data.frame(
    year = 1959:2013,
    foi_index = chikvfoii[,2]
  ),
  iter = 3000
)
plot_seromodel(seromodel=seromodel,serosurvey=chikv)
plot_foi_estimates(seromodel=seromodel,serosurvey=chikv)

#***********************************************************************


n_sample <- round(c(61,143,144,124,114,116,73,64,53,48,59)/2,0)
survey_year <- c(2014,2014,2014,2014,2014,2014,2014,2014,2014,2014,2014)
n_seropositive <- round(c(1,4,7,16,36,49,34,38,33,30,43)/2,0)
age_min <- c(1,5,10,15,20,25,30,35,40,45,50)
age_max <- c(4,9,14,19,24,29,34,39,44,49,55)

chikv<-cbind.data.frame(n_sample,survey_year,n_seropositive,age_min,age_max)

chikvfoii <- get_foi_index(chikv,group_size=5,model_type="time")

seromodel <- fit_seromodel(
  serosurvey = chikv,
  model_type = "time",
  foi_index = data.frame(
    year = 1959:2013,
    foi_index = chikvfoii[,2]
  ),
  iter = 3000
)
plot_seromodel(seromodel=seromodel,serosurvey=chikv)
plot_foi_estimates(seromodel=seromodel,serosurvey=chikv)


#*********************************************************
# if we wanted to know how exposure may have changed
# over the last 30 years
# what would that look like?
# 2025 to 1995

years <- seq(1970,2025,1)
fois <- rep(0.00,length(years))
fois[20:30] <- 0.04
fois[50:56] <- 0.06

yf <- cbind.data.frame(years,fois)

# if a person is 25 in 2025
# work out year of birth
pers1YrBirth <- 2025 - 25
# born in 2000 - so probability immune
# sum of the probabilities over the years
yf$probSus <- exp(-yf$fois)


ages <- c(1:55)
birthYears <- 2025 - ages
probExp <- numeric(length(birthYears))
for(i in 1:length(birthYears)){probExp[i] <- 1-prod(yf$probSus[yf$years>birthYears[i]])}
plot(ages,probExp)

denom <- 10

ageInf <- cbind.data.frame(ages,probExp)
ageInf$seropos <- sapply(1:length(ageInf$probExp),function(x){
  sum((rbinom(n=denom,size=1,prob=ageInf$probExp[x])))
})
ageInf$n <- denom


labs <- c(paste(seq(0, 55, by = 5), seq(0 + 5 - 1, 60- 1, by = 5),
                sep = "-"), paste(60, "+", sep = ""))
labs

ageInf$AgeGroup <- cut(ageInf$ages, breaks = c(seq(0, 60, by = 5),Inf), labels = labs, right = FALSE)

simsero <- ddply(ageInf,.(AgeGroup),summarise
                 ,n_seropositive=sum(seropos)
                 ,n_sample=sum(n))
ages <- do.call(rbind,(str_split(pattern="-",simsero$AgeGroup)))
simsero$age_min <- as.numeric(ages[,1])
simsero$age_max <- as.numeric(ages[,2])
simsero <- simsero[,-1]
simsero$survey_year <- 2025

sum(simsero$n_sample)

seromodel <- fit_seromodel(
  serosurvey = simsero
)

plot_seromodel(seromodel=seromodel,serosurvey=simsero)
plot_foi_estimates(seromodel=seromodel,serosurvey=simsero)


foii <- get_foi_index(simsero,group_size=5,model_type="time")

seromodel <- fit_seromodel(
  serosurvey = simsero,
  model_type = "time",
  foi_index = data.frame(
    year = 1966:2024,
    foi_index = foii[,2]
  ),
  iter = 3000
)
plot_seromodel(seromodel=seromodel,serosurvey=simsero)
plot_foi_estimates(seromodel=seromodel,serosurvey=simsero)

