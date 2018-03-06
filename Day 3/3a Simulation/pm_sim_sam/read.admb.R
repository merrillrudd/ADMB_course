read.admb <-
function(ifile)
{	
	ret=read.fit(ifile)
	
	fn=paste(ifile,'.rep', sep='')
	A=read.rep(fn)
	A$fit=ret
	
	pfn=paste(ifile,'.psv',sep='')
	if(file.exists(pfn))
		A$post.samp=read.psv(pfn)
	
	return(A)
}

read.fit <-
function(ifile)
{
	# __Example:             
	#	file <-("~/admb/simple")
	#	A <- reptoRlist(file)
	#	Note there is no extension on the file name.
	
	## The following is a contribution from:
	## Anders Nielsen that reads the par & cor files.
	ret<-list() 
	parfile<-as.numeric(scan(paste(ifile,'.par', sep=''),   
	 what='', n=16, quiet=TRUE)[c(6,11,16)]) 
	ret$nopar<-as.integer(parfile[1]) 
	ret$nlogl<-parfile[2] 
	ret$maxgrad<-parfile[3] 
	file<-paste(ifile,'.cor', sep='') 
	lin<-readLines(file) 
	ret$npar<-length(lin)-2 
	ret$logDetHess<-as.numeric(strsplit(lin[1], '=')[[1]][2]) 
	sublin<-lapply(strsplit(lin[1:ret$npar+2], ' '),function(x)x[x!='']) 
	ret$names<-unlist(lapply(sublin,function(x)x[2])) 
	ret$est<-as.numeric(unlist(lapply(sublin,function(x)x[3]))) 
	ret$std<-as.numeric(unlist(lapply(sublin,function(x)x[4]))) 
	ret$cor<-matrix(NA, ret$npar, ret$npar) 
	corvec<-unlist(sapply(1:length(sublin), function(i)sublin[[i]][5:(4+i)])) 
	ret$cor[upper.tri(ret$cor, diag=TRUE)]<-as.numeric(corvec) 
	ret$cor[lower.tri(ret$cor)] <- t(ret$cor)[lower.tri(ret$cor)] 
	ret$cov<-ret$cor*(ret$std%o%ret$std)
	return(ret)
}

read.rep <- 
function(fn)
{
	# The following reads a report file
	# Then the 'A' object contains a list structure
	# with all the elemements in the report file.
	# In the REPORT_SECTION of the AMDB template use 
	# the following format to output objects:
	#  	report<<"object \n"<<object<<endl;
	#
	# The part in quotations becomes the list name.
	# Created By Steven Martell
	options(warn=-1)  #Suppress the NA message in the coercion to double
	
	
	ifile=scan(fn,what="character",flush=TRUE,blank.lines.skip=FALSE,quiet=TRUE)
	idx=sapply(as.double(ifile),is.na)
	vnam=ifile[idx] #list names
	nv=length(vnam) #number of objects
	A=list()
	ir=0
	for(i in 1:nv)
	{
		ir=match(vnam[i],ifile)
		if(i!=nv) irr=match(vnam[i+1],ifile) else irr=length(ifile)+1 #next row
		dum=NA
		if(irr-ir==2) dum=as.double(scan(fn,skip=ir,nlines=1,quiet=TRUE,what=""))
		if(irr-ir>2) dum=as.matrix(read.table(fn,skip=ir,nrow=irr-ir-1,fill=TRUE))

		if(is.numeric(dum))#Logical test to ensure dealing with numbers
		{
			A[[vnam[i]]]=dum
		}
	}
	options(warn=0)
	
	return(A)
}

read.psv <-
function(fn, nsamples=10000)
{
	#This function reads the binary output from ADMB
	#-mcsave command line option.
	#fn = paste(ifile,'.psv',sep='')
	filen <- file(fn, "rb")
	nopar <- readBin(filen, what = integer(), n = 1)
	mcmc <- readBin(filen, what = numeric(), n = nopar * nsamples)
	mcmc <- matrix(mcmc, byrow = TRUE, ncol = nopar)
	close(filen)
	return(mcmc)
}

gletter <-
function(i=1){
	usr=par("usr"); inset.x=0.95*(usr[2]-usr[1]); inset.y=0.05*(usr[4]-usr[3])
		text(usr[1]+inset.x,usr[4]-inset.y,paste("(",letters[i],")",sep=""),cex=1.,font=1)
	}
	
	
require(MASS)
require( KernSmooth)
fried.egg<-
function(xx,yy,...)
{
	bw=25
	bwx=diff(extendrange(xx))/bw; bwy=diff(extendrange(yy))/bw
	#bwx=(max(xx)-min(xx))/bw
	#bwy=(max(yy)-min(yy))/bw
	est <- bkde2D(cbind(xx,yy),bandwidth=c(bwx,bwy),gridsize=c(81, 81))
	est$fhat=est$fhat/max(est$fhat)
	#plot(xx,yy,pch=".",col="dark grey",xlab=NA,ylab=NA,type="n")
	#text(max(xx),max(yy),labels="D",adj=c(1,1))
	lvs=c(0.05,0.25,0.75,0.95)
	maxct=max(lvs)
	nlvs=length(lvs)
	thelines=contourLines(est$x1,est$x2,est$fhat,levels=lvs)
	iclr=colr("black", 0.3)
	polygon(thelines[[nlvs-3]]$x,thelines[[nlvs-3]]$y,col=iclr,border=iclr,lwd=1)
	iclr=colr("snow", 0.9)
	polygon(thelines[[nlvs-2]]$x,thelines[[nlvs-2]]$y,col=iclr,border=iclr,lwd=2)
	iclr=colr("yellow", 0.9)
	polygon(thelines[[nlvs-1]]$x,thelines[[nlvs-1]]$y,col=iclr,border=iclr,lwd=3)
	polygon(thelines[[nlvs]]$x,thelines[[nlvs]]$y,col="lightyellow",border="yellow",lwd=1)
	#contour(est$x1,est$x2,est$fhat,drawlabels=T,add=T,levels=lvs,lty=1,lwd=1,labcex= 0.7)
	#Add salt and pepper
	#xi=sample(1:length(xx),300)
	#points(xx[xi],yy[xi],pch=".",col=grey(0:10/10))
}
