LMC=read.psv("pm.psv")

recanom=LMC[,6:length(LMC[1,])]
dfile=recanom
#create vectors for plotting and histogram function
tms=1:length(dfile[1,])
rng=range(dfile)*c(1.05,1.05)
brks=seq(rng[1],rng[2],length=20)
#define functions that are nedded
normalize=function(x) x/sum(x) #ensures a vector sums to 1
getprobs=function(x) hist(x,breaks=brks,plot=FALSE)$density  # uses the hist function of get the densities 
ncol=21 #number of colors
#Define colour spectrum
rmp=c("red","orange","yellow","lightgreen","green","steelblue","slateblue3")
clr = colorRampPalette(rmp, space = "Lab")
cols=rev(clr(ncol))
#bin
binned=apply(dfile,2,FUN=getprobs) # apply the getprobs function to each column of the dfile
#Normalize data for plotting
normdfile=apply(binned,2,FUN=normalize) # although each column should sum to 1 make sure by applying the normalize function to each column. 
maxp=max(normdfile) # needed to scale plotting
if (maxp==1) maxp=1.1
#contour(x=1:100,y=0:300,z=t(normdfile),levels=seq(0.02,0.98,length=40),
#   drawlabels=FALSE,xaxs="i",yaxs="i",col=cols)  # this one make pretty colored lines
filled.contour(x=tms,y=seq(rng[1],rng[2],length=length(binned[,1])),z=t(normdfile),levels=seq(0.001,maxp,length=ncol),xaxs="i",yaxs="i",col=cols)
