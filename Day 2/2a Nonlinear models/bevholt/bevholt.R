wd <- "C:\\merrill\\ADMB_course"
project_dir <- file.path(wd, "Day 2", "2a Nonlinear models", "bevholt")
setwd(project_dir)

### read in results from report file using code in tools.R
source(file.path(wd, "tools.R"))

## read in data
nrow <- readVec("# n", file="bevholt1.dat")
sr_data <- readMat("#     S       R", file="bevholt1.dat", nrow=nrow)
Sobs <- sr_data[,1][order(sr_data[,1])]
Robs <- sr_data[,2][order(sr_data[,1])]

## plot observed data
plot(x=Sobs, y=Robs, pch=19, xpd=NA, xlab="Spawners", ylab="Recruits", xaxs="i", yaxs="i")

## predicted values
bh_fit <- readFit(file="bevholt1")
	
	## find length estimates
	Rpred_est <- bh_fit$est[which(bh_fit$names == "Rpred")][order(sr_data[,1])]
	Rpred_sd <- bh_fit$std[which(bh_fit$names == "Rpred")][order(sr_data[,1])]
	Rpred_lcl <- Rpred_est - 1.96 * Rpred_sd
	Rpred_ucl <- Rpred_est + 1.96 * Rpred_sd

## plot uncertainty in predicted recruits
polygon(x=c(Sobs, rev(Sobs)), y=c(Rpred_lcl, rev(Rpred_ucl)), col="#AA000050", border=NA)

## plot predicted recruitment estimates
lines(x=Sobs, y=Rpred_est, col="#AA0000")
