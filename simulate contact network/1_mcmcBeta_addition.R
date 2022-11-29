set.seed(9876)
room_table <- read.csv('output/pre-calculation for simulate contact matrix/S_matrix_room_count_RandT_new.csv')
room_table <- room_table[complete.cases(room_table),]
logl = function(dt,pars)
{
  
  LL = sum(dbeta(dt,pars$alpha,pars$beta,log=TRUE))
  pars$LL = LL
  return(pars)
}

mh = function(oldp,newp,dt)
{
  reject = FALSE
  if(newp$alpha<0)reject=TRUE
  if(newp$beta<0)reject=TRUE
  if(!reject)
  {
    newp = logl(dt,newp)
    la = newp$LL - oldp$LL
    lu = -rexp(1)
    if(lu>la)reject=TRUE
  }
  if(reject)return(oldp)
  return(newp)
}

mcmc <- function(data,pars,MCMCITS)
{
  
  storage=list(
    alpha=rep(0,MCMCITS),
    beta=rep(0,MCMCITS),
    LL=rep(0,MCMCITS)
  )
  SUBITS=10
  for(iteration in 1:MCMCITS)
  {
    if(iteration%%100==0)cat(iteration,'in',MCMCITS,'\n')
    for(subit in 1:SUBITS)
    {
      oldpars = pars
      pars$alpha = rnorm(1,pars$alpha,1)
      pars$beta = rnorm(1,pars$beta,1) 
      pars = mh(oldpars,pars,data)
    }
    storage$alpha[iteration] = pars$alpha
    storage$beta[iteration] = pars$beta
    storage$LL[iteration] = pars$LL
    
  }
  
  storage <- as.data.frame(storage)
  storage <- storage[-c(1:1000),]
  storage
}



MCMCITS=100000
dat <- room_table%>%filter(dorm=="A")
dat$prop_Rc[which(dat$prop_Rc==1)] <- 0.9999
dat$prop_Rc[which(dat$prop_Rc==0)] <- 0.0001
# dat$prop_Tc[which(dat$prop_Tc==1)] <- 0.9999
# dat$prop_Tc[which(dat$prop_Tc==0)] <- 0.0001
pars <- list(alpha=0.1,beta=0.1)
pars <- logl(dat$prop_Rc,pars)
storage_ARc<- mcmc(data=dat$prop_Rc,pars=pars,MCMCITS)
write.csv(storage_ARc,'working/simulate contact matrix/Acacia_regular_Room_BetaMCMC.csv',row.names = F)
mean(storage_ARc$alpha)
mean(storage_ARc$beta)

# pars <- list(alpha=1,beta=2)
# pars <- logl(dat$prop_Tc,pars)
# storage_ATc<- mcmc(data=dat$prop_Tc,pars=pars,MCMCITS)
# mean(storage_ATc$alpha)
# mean(storage_ATc$beta)
# write.csv(storage_ATc,'working/simulate contact matrix/Acacia_transient_Room_BetaMCMC.csv',row.names = F)

dat <- room_table%>%filter(dorm=="C")
dat$prop_Rc[which(dat$prop_Rc==1)] <- 0.9999
dat$prop_Rc[which(dat$prop_Rc==0)] <- 0.0001
# dat$prop_Tc[which(dat$prop_Tc==1)] <- 0.9999
# dat$prop_Tc[which(dat$prop_Tc==0)] <- 0.0001
pars <- list(alpha=0.1,beta=0.1)
pars <- logl(dat$prop_Rc,pars)
storage_CRc<- mcmc(data=dat$prop_Rc,pars=pars,MCMCITS)
mean(storage_CRc$alpha)
mean(storage_CRc$beta)
write.csv(storage_CRc,'working/simulate contact matrix/Cassia_regular_Room_BetaMCMC.csv',row.names = F)


# pars <- list(alpha=1,beta=2)
# pars <- logl(dat$prop_Tc,pars)
# storage_CTc<- mcmc(data=dat$prop_Tc,pars=pars,MCMCITS)
# mean(storage_CTc$alpha)
# mean(storage_CTc$beta)
# write.csv(storage_CTc,'working/simulate contact matrix/Cassia_transient_Room_BetaMCMC.csv',row.names = F)




level_table <- read.csv('output/S_matrix_level_count_RandT_new.csv')
level_table <- level_table[complete.cases(level_table),]
dat <- level_table%>%filter(dorm=="A")
dat <- dat[-which(dat$prop_Rc==0),]
dat$prop_Rc[which(dat$prop_Rc==1)] <- 0.9999
dat$prop_Rc[which(dat$prop_Rc==0)] <- 0.0001
dat$prop_Tc[which(dat$prop_Tc==1)] <- 0.9999
dat$prop_Tc[which(dat$prop_Tc==0)] <- 0.0001
pars <- list(alpha=1,beta=2)
pars <- logl(dat$prop_Tc,pars)

storage_AFRc<- mcmc(data=dat$prop_Rc,pars=pars,MCMCITS)
plot(storage_AFRc$alpha,type='l')
plot(storage_AFRc$beta,type='l')
rbeta(6,mean(storage_AFRc$alpha),mean(storage_AFRc$beta))
write.csv(storage_AFRc,'working/simulate contact matrix/Acacia_regular_Level_BetaMCMC.csv',row.names = F)

storage_AFTc<- mcmc(data=dat$prop_Tc,pars=pars,MCMCITS)
plot(storage_AFTc$alpha,type='l')
plot(storage_AFTc$beta,type='l')
rbeta(6,mean(storage_AFTc$alpha),mean(storage_AFTc$beta))
write.csv(storage_AFTc,'working/simulate contact matrix/Acacia_transient_Level_BetaMCMC.csv',row.names = F)



dat <- level_table%>%filter(dorm=="C")
dat <- dat[-which(dat$prop_Rc==0),]
storage_CFRc<- mcmc(data=dat$prop_Rc,pars=pars,MCMCITS)
# plot(storage_CFRc$alpha,type='l')
# plot(storage_CFRc$beta,type='l')
# rbeta(6,mean(storage_CFRc$alpha),mean(storage_CFRc$beta))
write.csv(storage_CFRc,'working/simulate contact matrix/Cassia_regular_Level_BetaMCMC.csv',row.names = F)


storage_CFTc<- mcmc(data=dat$prop_Tc,pars=pars,MCMCITS)
# plot(storage_AFTc$alpha,type='l')
# plot(storage_AFTc$beta,type='l')
# rbeta(6,mean(storage_CFTc$alpha),mean(storage_CFTc$beta))
write.csv(storage_CFTc,'working/simulate contact matrix/Cassia_transient_Level_BetaMCMC.csv',row.names = F)


blk_table <- read.csv('output/S_matrix_blk_count_RandT_new.csv')
blk_table <- blk_table[complete.cases(blk_table),]
dat <- blk_table%>%filter(dorm=="A")
storage_ABRc<- mcmc(data=dat$prop_Rc,pars=pars,MCMCITS)
plot(storage_ABRc$alpha,type='l')
plot(storage_ABRc$beta,type='l')
write.csv(storage_ABRc,'working/simulate contact matrix/Acacia_regular_Blk_BetaMCMC.csv',row.names = F)
storage_ABTc<- mcmc(data=dat$prop_Tc,pars=pars,MCMCITS)
write.csv(storage_ABTc,'working/simulate contact matrix/Acacia_transient_Blk_BetaMCMC.csv',row.names = F)


dat <- blk_table%>%filter(dorm=="C")
storage_CBRc<- mcmc(data=dat$prop_Rc,pars=pars,MCMCITS)
plot(storage_CBRc$alpha,type='l')
plot(storage_CBRc$beta,type='l')
write.csv(storage_CBRc,'working/simulate contact matrix/Cassia_regular_Blk_BetaMCMC.csv',row.names = F)
storage_CBTc<- mcmc(data=dat$prop_Tc,pars=pars,MCMCITS)
write.csv(storage_CBTc,'working/simulate contact matrix/Cassia_transient_Blk_BetaMCMC.csv',row.names = F)



# reply review, mcmc convergence ------------------------------------------

# Geweke diagnostic
install.packages('coda')
library(coda)
# temp <- geweke.diag(storage_ARc$alpha,frac1 = 0.1,frac2 = 0.2)
# pnorm(abs(temp$z),lower.tail=FALSE)*2
# dat1 <- read.csv("working/simulate contact matrix/before Nov 2022/Acacia_regular_Level_BetaMCMC.csv")
temp11 <- geweke.diag(storage_ARc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp11$z),lower.tail=FALSE)*2

temp12 <- geweke.diag(storage_ARc$beta[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp12$z),lower.tail=FALSE)*2

temp21 <- geweke.diag(storage_AFRc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp21$z),lower.tail=FALSE)*2

temp22 <- geweke.diag(storage_AFRc$beta[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp22$z),lower.tail=FALSE)*2

temp31 <- geweke.diag(storage_ABRc$alpha[c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp31$z),lower.tail=FALSE)*2

temp32 <- geweke.diag(storage_ABRc$beta[-c(1:60000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp32$z),lower.tail=FALSE)*2


tempc11 <- geweke.diag(storage_CRc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc11$z),lower.tail=FALSE)*2

tempc12 <- geweke.diag(storage_CRc$beta[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc12$z),lower.tail=FALSE)*2

tempc21 <- geweke.diag(storage_CFRc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc21$z),lower.tail=FALSE)*2

tempc22 <- geweke.diag(storage_CFRc$beta[c(1:51000)],frac1 = 0.1,frac2 = 0.5)
pnorm(abs(tempc22$z),lower.tail=FALSE)*2

tempc31 <- geweke.diag(storage_CBRc$alpha[c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc31$z),lower.tail=FALSE)*2

tempc32 <- geweke.diag(storage_CBRc$beta[c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc32$z),lower.tail=FALSE)*2






temp11 <- geweke.diag(storage_ATc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp11$z),lower.tail=FALSE)*2

temp12 <- geweke.diag(storage_ATc$beta[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp12$z),lower.tail=FALSE)*2

temp21 <- geweke.diag(storage_AFTc$alpha[-c(1:61000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp21$z),lower.tail=FALSE)*2

temp22 <- geweke.diag(storage_AFTc$beta[-c(1:75000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp22$z),lower.tail=FALSE)*2

temp31 <- geweke.diag(storage_ABTc$alpha[c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp31$z),lower.tail=FALSE)*2

temp32 <- geweke.diag(storage_ABTc$beta[c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(temp32$z),lower.tail=FALSE)*2


tempc11 <- geweke.diag(storage_CTc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc11$z),lower.tail=FALSE)*2

tempc12 <- geweke.diag(storage_CTc$beta[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc12$z),lower.tail=FALSE)*2

tempc21 <- geweke.diag(storage_CFTc$alpha[-c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc21$z),lower.tail=FALSE)*2

tempc22 <- geweke.diag(storage_CFTc$beta[c(1:51000)],frac1 = 0.1,frac2 = 0.5)
pnorm(abs(tempc22$z),lower.tail=FALSE)*2

tempc31 <- geweke.diag(storage_CBTc$alpha[c(10000:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc31$z),lower.tail=FALSE)*2

tempc32 <- geweke.diag(storage_CBTc$beta[c(1:51000)],frac1 = 0.1,frac2 = 0.6)
pnorm(abs(tempc32$z),lower.tail=FALSE)*2
