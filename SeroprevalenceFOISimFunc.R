


#************seroprevalence as a function of FOI and age*****************
propImmFunc <- function(lambda, age){
  return(1 - exp(-lambda*age))
}
propImmFuncV <- Vectorize(propImmFunc)
#********************************************************************

#************function to simulate seroprevalence data******************
simulateSeroprevalence <- function(lowLambda=0.02      # lower force of infection
                                   ,highLambda=0.03    # high force of infection
                                   ,sdLogFOI=0.2       # to add in variation between sites
                                   ,n.village=50       
                                   ,n.hh=10
                                   ,people.in.household=2
                                   ,ageMin=1
                                   ,ageMax=70
                                   ,ageMean=25
                                   ,ageSD=3
                                   ,distAge="uniform"  
                                   ){


  dat <- expand.grid(hh = 1:n.hh, village = 1:n.village,people=1:people.in.household)
  # allocate villages to high and low prevalence in 1:1 ratio 
  dat$risk.level <- dat$village %% 2 - 0.5
  dat$FOI <- NA
  dat$FOI[dat$risk.level==-0.5] <- lowLambda
  dat$FOI[dat$risk.level==0.5] <- highLambda

  # add village-level FOI variation
  dat$villageFOI <- rlnorm(length(dat$village),mean=log(dat$FOI),sd=sdLogFOI)
  # add household-level FOI variation
  dat$hhFOI <- rlnorm(length(dat$village),mean=log(dat$villageFOI),sd=sdLogFOI)
  
  # randomly assign an age  - simplest option uniform
  dat$age <- NA
  if(distAge=="uniform"){
  dat$age<-sapply(1:length(dat$age),function(x){round(runif(1,ageMin,ageMax),0)})}else{
    dat$age<-sapply(1:length(dat$age),function(x){round(rnorm(1,25,0))})
  }
  # calculate true prob immune
  dat$probImm <- propImmFuncV(lambda=dat$hhFOI,age=dat$age) #

  # from 'true' proportion immune to random variation                                    
  dat$Infected <- sapply(1:length(dat$probImm),function(x){
   rbinom(n=1,size=1,prob=dat$probImm[x])
  })

return(dat)
}
#************************************************************************



dat <- simulateSeroprevalence(lowLambda=0.01
                       ,highLambda=0.02
                       ,n.village=5
                       ,n.hh=5
                       ,people.in.household=2
                       ,ageMin=1
                       ,ageMax=60
                       ,ageMean=20
                       ,sdLogFOI=0.5
                       ,ageSD=2
                       ,distAge="uniform")


ggplot(dat) +
  geom_point(aes(x=age,y=probImm,col=risk.level))



fit <- glm(Infected~risk.level
           ,family=binomial(link="cloglog")
           ,offset=log(age)
           ,data=dat)

fit0 <- update(fit, ~ . - risk.level)
pval <- lrtest(fit, fit0)[2,5]
coefs <- coef(fit)  
coefs
