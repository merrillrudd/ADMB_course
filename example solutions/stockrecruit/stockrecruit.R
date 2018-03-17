wd <- "C:\\merrill\\ADMB_course"
project_dir <- file.path(wd, "Day2", "2a_Nonlinear")
d1 <- file.path(project_dir, "stockrecruit")
setwd(d1)

### read in results from report file using code in tools.R
source(file.path(wd, "tools.R"))

######################################
## read in data for one Bev-Holt model
######################################
nrow <- readVec("# n", file="bevholt1.dat")
sr_data <- readMat("#     S       R", file="bevholt1.dat", nrow=nrow)
Sobs <- sr_data[,1][order(sr_data[,1])]
Robs <- sr_data[,2][order(sr_data[,1])]

## predicted values
fit <- readFit(file="bevholt1")
	
	## find length estimates
	Rpred_est <- fit$est[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_sd <- fit$std[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_lcl <- Rpred_est - 1.96 * Rpred_sd
	Rpred_ucl <- Rpred_est + 1.96 * Rpred_sd

## plot observed data
plot(x=Sobs, y=Robs, ylim=c(0, max(c(Robs, Rpred_ucl))*1.1), pch=19, xpd=NA, xlab="Spawners", ylab="Recruits", xaxs="i", yaxs="i")

## plot uncertainty in predicted recruits
polygon(x=c(Sobs, rev(Sobs)), y=c(Rpred_lcl, rev(Rpred_ucl)), col=paste0(gray(0.3),"50"), border=NA)

## plot predicted recruitment estimates
lines(x=Sobs, y=Rpred_est, col=gray(0.3), lwd=2)

###############################################
## plot estimates for multiple Bev-Holt models
###############################################
library(RColorBrewer)
## add models to compare

modnames <- c("bevholt1", "bevholt2", "bevholt3")
Rpred_est <- Rpred_sd <- Rpred_lcl <- Rpred_ucl <- nll <- npar <- aic <- list()
for(i in 1:length(modnames)){
	## predicted values
	fit <- readFit(file=modnames[i])
	nobs <- readVec("# n", paste0(modnames[i],".dat"))
	nll <- fit$nlogl
	npar <- fit$npar
	aic[[i]] <- 2 * npar - 2 * nll + ((2 * npar^2 + 2 * npar)/(nobs - npar - 1))
	
	## find length estimates
	Rpred_est[[i]] <- fit$est[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_sd[[i]] <- fit$std[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_lcl[[i]] <- Rpred_est[[i]] - 1.96 * Rpred_sd[[i]]
	Rpred_ucl[[i]] <- Rpred_est[[i]] + 1.96 * Rpred_sd[[i]]
}


## plot observed data
plot(x=Sobs, y=Robs, ylim=c(0, max(c(Robs))*1.1), pch=19, xpd=NA, xlab="Spawners", ylab="Recruits", xaxs="i", yaxs="i")

cols <- brewer.pal(length(modnames), "Set1")
for(i in 1:length(modnames)){
	polygon(x=c(Sobs, rev(Sobs)), y=c(Rpred_lcl[[i]], rev(Rpred_ucl[[i]])), col=paste0(cols[i],"50"), border=NA)
	lines(x=Sobs, y=Rpred_est[[i]], col=cols[i], lwd=2)
}
legend("topleft", legend=modnames, col=cols[1:length(modnames)], pch=15)

#################################################
## plot estimates for Bev-Holt and Ricker models
#################################################
## Rmax parameterization
modnames <- c("bevholt1", "ricker1")
Rpred_est <- Rpred_sd <- Rpred_lcl <- Rpred_ucl <- nll <- npar <- aic <- list()
for(i in 1:length(modnames)){
	## predicted values
	fit <- readFit(file=modnames[i])
	nobs <- readVec("# n", paste0(modnames[i],".dat"))
	nll <- fit$nlogl
	npar <- fit$npar
	aic[[i]] <- 2 * npar - 2 * nll + ((2 * npar^2 + 2 * npar)/(nobs - npar - 1))
	
	## find length estimates
	Rpred_est[[i]] <- fit$est[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_sd[[i]] <- fit$std[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_lcl[[i]] <- Rpred_est[[i]] - 1.96 * Rpred_sd[[i]]
	Rpred_ucl[[i]] <- Rpred_est[[i]] + 1.96 * Rpred_sd[[i]]
}


## plot observed data
plot(x=Sobs, y=Robs, ylim=c(0, max(c(Robs))*1.1), pch=19, xpd=NA, xlab="Spawners", ylab="Recruits", xaxs="i", yaxs="i")

cols <- brewer.pal(length(modnames), "Set1")
for(i in 1:length(modnames)){
	polygon(x=c(Sobs, rev(Sobs)), y=c(Rpred_lcl[[i]], rev(Rpred_ucl[[i]])), col=paste0(cols[i],"50"), border=NA)
	lines(x=Sobs, y=Rpred_est[[i]], col=cols[i], lwd=2)
}
legend("topleft", legend=modnames, col=cols[1:length(modnames)], pch=15)


## alpha-beta parameterization
modnames <- c("bevholt2", "ricker2")
Rpred_est <- Rpred_sd <- Rpred_lcl <- Rpred_ucl <- nll <- npar <- aic <- list()
for(i in 1:length(modnames)){
	## predicted values
	fit <- readFit(file=modnames[i])
	nobs <- readVec("# n", paste0(modnames[i],".dat"))
	nll <- fit$nlogl
	npar <- fit$npar
	aic[[i]] <- 2 * npar - 2 * nll + ((2 * npar^2 + 2 * npar)/(nobs - npar - 1))
	
	## find length estimates
	Rpred_est[[i]] <- fit$est[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_sd[[i]] <- fit$std[which(fit$names == "Rec_pred")][order(sr_data[,1])]
	Rpred_lcl[[i]] <- Rpred_est[[i]] - 1.96 * Rpred_sd[[i]]
	Rpred_ucl[[i]] <- Rpred_est[[i]] + 1.96 * Rpred_sd[[i]]
}

## plot observed data
plot(x=Sobs, y=Robs, ylim=c(0, max(c(Robs))*1.1), pch=19, xpd=NA, xlab="Spawners", ylab="Recruits", xaxs="i", yaxs="i")

cols <- brewer.pal(length(modnames), "Set1")
for(i in 1:length(modnames)){
	polygon(x=c(Sobs, rev(Sobs)), y=c(Rpred_lcl[[i]], rev(Rpred_ucl[[i]])), col=paste0(cols[i],"50"), border=NA)
	lines(x=Sobs, y=Rpred_est[[i]], col=cols[i], lwd=2)
}
legend("topleft", legend=modnames, col=cols[1:length(modnames)], pch=15)
