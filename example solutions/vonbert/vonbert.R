wd <- "C:\\merrill\\ADMB_course"
project_dir <- file.path(wd, "Day2", "2a_Nonlinear", "vonbert")
setwd(project_dir)

### read in results from report file using code in tools.R
source(file.path(wd, "tools.R"))

## read in data
nrow <- readVec("# n", file="vonbert2.dat")
vb_data <- readMat("# Age Length", file="vonbert2.dat", nrow=nrow)
ages <- vb_data[,1]
Lobs <- vb_data[,2]

## predicted values
vb_fit <- readFit(file="vonbert2")
	
	## find length estimates
	Lpred_est <- vb_fit$est[which(vb_fit$names == "Lpred")]
	Lpred_sd <- vb_fit$std[which(vb_fit$names == "Lpred")]
	Lpred_lcl <- Lpred_est - 1.96 * Lpred_sd
	Lpred_ucl <- Lpred_est + 1.96 * Lpred_sd

## plot uncertainty in predicted length
plot(x=1, y=1, type="n", xlim=c(min(ages), max(ages)), ylim=c(0, max(Lpred_ucl)*1.1), xaxs="i", yaxs="i", xaxt="n", yaxt="n", xlab="Age", ylab="Length")
polygon(x=c(ages, rev(ages)), y=c(Lpred_lcl, rev(Lpred_ucl)), col="#AA000050", border=NA)

## plot predicted length estimates
lines(x=ages, y=Lpred_est, col="#AA0000")

## plot observed data
points(x=ages, y=Lobs, pch=19, xpd=NA)

## add axes
axis(1, pretty(ages))
axis(2, pretty(Lobs), las=2)