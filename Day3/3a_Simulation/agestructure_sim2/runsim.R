source("read.admb.R")
run.Simulation=function(N=10)
{
    theta <<- NULL
    for(i in sample(1:1000,N))
    {
        arg = paste("./sca -sim", i)
        system(arg)
        print(arg)
        P=read.fit("sca")
        theta<<-rbind(theta, P$est[1:6])
    }

}
run.Simulation(10)
tvalues=c(12000,10,12000,0.5,2,0.75)
names=c("ro","cr","rbar","fbar","ahat","ghat")
trtheta=cbind(exp(theta[,1:4]),theta[,5:6])
bias=trtheta*0
for(i in 1:6)
{
    bias[,i]=(trtheta[,i]-tvalues[i])/tvalues[i]*100
}
    
boxplot(bias,names=names)
