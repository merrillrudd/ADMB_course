setwd("C:\\merrill\\ADMB_course\\Day3\\3a_Simulation\\agestructure_sim1")

library(R2admb)

source("read.admb.R")
source("C:\\merrill\\admb_course\\tools.R")

run.Simulation=function(N=10)
{
    out <- NULL 
    theta <<- NULL
    seeds <- sample(1:1000,N)
    estF <- NULL
    for(i in 1:length(seeds))
    {
        arg = paste("-sim", seeds[i])
        run_admb("SCA", extra.args=arg)
        print(arg)
        P=read.fit("SCA")
        theta<<-rbind(theta, P$est[1:6])
        estF[[i]] <- readVec("mfexp(log_fbar+log_ft_dev)", file="SCA.rep")
    }
    out$theta <- theta
    out$estF <- estF
    return(out)
}
compile_admb("SCA")
sim <- run.Simulation(10)

## ro= initial recruitment
## cr= recruitment compensation ratio
## rbar= average recruitment over time series
## fbar= average fishing mortality over time seires
## ahat= age at 50% selectivity
## ghat= shape of logistic selectivity curve
names=c("ro","cr","rbar","fbar","ahat","ghat")
trueF <- readVec(string="#simF input fishing mortality", file="Simdata.ctl")
tvalues=c(12000,10,12000,mean(trueF),2,0.75)

trtheta=cbind(exp(theta[,1:4]),theta[,5:6])
bias=trtheta*0
for(i in 1:6)
{
    bias[,i]=(trtheta[,i]-tvalues[i])/tvalues[i]*100
}
    
boxplot(bias,names=names)
abline(h=0, col="red")

plot(trueF, ylim=c(0,0.6))
for(i in 1:length(sim$estF)){
	lines(sim$estF[[i]])
}