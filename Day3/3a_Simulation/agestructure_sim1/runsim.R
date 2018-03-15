setwd("C:\\merrill\\ADMB_course\\Day3\\3a_Simulation\\agestructure_sim1")

library(R2admb)

source("read.admb.R")

run.Simulation=function(N=10)
{
    theta <<- NULL
    for(i in sample(1:1000,N))
    {
        arg = paste("-sim", i)
        run_admb("SCA", extra.args=arg)
        print(arg)
        P=read.fit("SCA")
        theta<<-rbind(theta, P$est[1:6])
    }

}
compile_admb("SCA")
run.Simulation(50)

## ro= initial recruitment
## cr= recruitment compensation ratio
## rbar= average recruitment over time series
## fbar= average fishing mortality over time seires
## ahat= age at 50% selectivity
## ghat= shape of logistic selectivity curve
names=c("ro","cr","rbar","fbar","ahat","ghat")
tvalues=c(12000,10,12000,0.5,2,0.75)
trtheta=cbind(exp(theta[,1:4]),theta[,5:6])
bias=trtheta*0
for(i in 1:6)
{
    bias[,i]=(trtheta[,i]-tvalues[i])/tvalues[i]*100
}
    
boxplot(bias,names=names)
abline(h=0, col="red")
