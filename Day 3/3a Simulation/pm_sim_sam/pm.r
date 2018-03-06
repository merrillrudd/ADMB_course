# R-code for the production model.
require(MCMCpack)
require(hacks)
source("read.admb.r")
A = read.admb("pm")
par(mfcol=c(2, 2), las=1)

with(A, {
	plot(iyr, bt, type="l", xlab="Year", ylab="Biomass")
	gletter(1)
	
	plot(iyr,ct/bt,type="o", xlab="Year", ylab="Exploitation rate")
	gletter(2)
	
	plot(iyr,it,type="l",xlab="Year",ylab="CPUE", 
		ylim=c(0, max(it, yt)))
	points(iyr, yt, pch=20)
	gletter(3)
	
	plot(iyr,nu,type="h",xlab="Year", ylab="CPUE residual")
	gletter(4)


	if(exists("post.samp")){
		ps=as.data.frame(post.samp)
		colnames(ps)=fit$names[1:fit$nopar]
		panel.hist <- function(x, ...)
		{
		    usr <- par("usr"); on.exit(par(usr))
		    par(usr = c(usr[1:2], 0, 1.5) )
		    h <- hist(x, plot = FALSE)
		    breaks <- h$breaks; nB <- length(breaks)
		    y <- h$counts; y <- y/max(y)
		    rect(breaks[-nB], 0, breaks[-1], y, col=colr("blue", 0.5))
		}
		
		n=dim(ps)[1]
		ix = (n/2):n		#discard 1st half for burnin
		pairs(ps[ix,1:5],pch=".",gap=0,col=c("orange"), 
			upper.panel=panel.smooth, 
			lower.panel=fried.egg, 
			diag.panel=panel.hist)
		
		
		#now plot priors and martinal posteriors.
		plot.marg <- function(ps,p1=0,p2=1,prior="dunif", mle=NULL, ...)
		{
			xl=range(ps)
			hist(ps,breaks=30, prob=T, col=colr("tan",0.5), 
			ylab="Probability density", main="", ...)
			
			fn=match.fun(prior)
			curve(unlist(lapply(x,fn,p1,p2)),
				xl[1],xl[2],add=T, col=4, lty=1, lwd=2)
			
			abline(v=mle, col=colr("red",0.5), lwd=5)
			
		}
		par(mfcol=c(3, 2),mar=c(4, 4, 2, 2), oma=c(2, 2, 2, 2))
		plot.marg(exp(ps[ix,1]),8.0, 0.5,"dlnorm", xlab="K", mle=exp(fit$est[1]))
		plot.marg(ps[ix,2],-1.38,0.51,"dlnorm", xlab="r", mle=fit$est[2])
		plot.marg(exp(ps[ix,3])*1e4,0,1,"dunif", xlab="q (1e4)", mle=exp(fit$est[3])*1e4)
		plot.marg(1/exp(ps[ix,4]),1.71, 0.0086,"dinvgamma",xlab="tau",mle=1/exp(fit$est[4]))
		plot.marg(1/exp(ps[ix,5]),3.79,0.0102,"dinvgamma", xlab="sig",mle=1/exp(fit$est[5]))
		plot.marg(ps[ix,2]*exp(ps[ix,1])/4, -15, 15, xlab="MSY", mle=fit$est[2]*exp(fit$est[1])/4)
	}

})


run.Simulation<-
function(N=10)
{
	theta <<- NULL
	for(i in sample(1:1000,N))
	{
		arg = paste("./pm -sim", i)
		system(arg)
		print(arg)
		P=read.fit("pm")
		theta<<-rbind(theta, P$est[1:5])
	}
	boxplot(theta)
}


