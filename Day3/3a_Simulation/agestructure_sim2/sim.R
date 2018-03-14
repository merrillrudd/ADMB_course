rm(list=ls())

# install.packages("R2admb")
library(R2admb)

main_dir <- "C:\\merrill\\ADMB_course"
proj_dir <- file.path(main_dir, "Day3", "3a_Simulation", "agestructure_sim")
setwd(proj_dir)

source(file.path(main_dir, "tools.R"))

######################################
## Setup R2admb
######################################
files <- list.files(file.path(proj_dir))

## find tpl
tpl <- file.path(proj_dir, files[which(grepl(".TPL", files))])
name <- strsplit(files[grepl("TPL",files)],".TPL")[[1]][1]

# ## compile
compile_admb(name, verbose=TRUE)

# ## check executable
exe <- paste0(name,".exe")
file.exists(exe)

####################################
## simulate deterministic population
#####################################
sim_det <- simPop(linf=60,
				vbk=0.2,
				t0=-0.01,
				lwa=3e-3,
				lwb=3,
				M=0.2,
				S50=4,
				S95=6,
				M50=3,
				M95=4,
				R0=1000,
				h=0.7,
				SigmaR=0,
				SigmaF=0,
				Nyear=20,
				InitialDepl=0.8,
				PropFcrash=0.8,
				EffortDyn="Two-way",
				seed=123)

par(mfrow=c(3,3))
plot(sim_det$F_t, type="l", lwd=2, ylim=c(0,1), xlab="Year", ylab="Fishing mortality")
plot(sim_det$C_t, type="l", lwd=2, ylim=c(0,max(sim_det$C_t)*1.1), xlab="Year", ylab="Catch")
plot(sim_det$CPUE_t, type="l", lwd=2, ylim=c(0, max(sim_det$CPUE_t)*1.1), xlab="Year", ylab="CPUE")
plot(sim_det$R_t, type="l", lwd=2, ylim=c(0,3000), xlab="Year", ylab="Recruitment")
plot(sim_det$SB_t/sim_det$SB0, type="l", lwd=2, ylim=c(0, max(sim_det$SB_t/sim_det$SB0)*1.1), xlab="Year", ylab="Relative spawning biomass")
plot(sim_det$N_t, type="l", lwd=2, ylim=c(0, max(sim_det$N_t)*1.1), xlab="Year", ylab="Abundance")
plot(sim_det$page[,1], type="h", lwd=4, ylim=c(0, max(sim_det$page)), col="gray", xlab="Age", ylab="Proportion in catch")
	lines(sim_det$page[,ncol(sim_det$page)], type="h", lwd=4)
	legend('topright', legend=c("First year", "Last year"), col=c("gray","black"), lwd=4)
plot(sim_det$L_a, type="l", lwd=2, ylim=c(0, max(sim_det$L_a)*1.1), xlab="Age", ylab="Length")
plot(sim_det$Mat_a, type="l", lwd=2, ylim=c(0,1), xlab="Age", ylab="Proportion")
	lines(sim_det$S_a, lwd=2, lty=2)
	legend("bottomright", legend=c("Maturity", "Selectivity"), lty=c(1,2), lwd=2)

#########################################
## create deterministic ADMB data file
#########################################
cat("# Age-structured simulation data file", file="AS_DET.dat", sep="\n")	
## numbers and integers - use sep="\n"
cat("# Number of years", file="AS_DET.dat", sep="\n", append=TRUE)
cat(sim_det$Nyear, file="AS_DET.dat", sep="\n", append=TRUE)
cat("# Number of ages", file="AS_DET.dat", sep="\n", append=TRUE)
cat(sim_det$Nage, file="AS_DET.dat", sep="\n", append=TRUE)
cat("# Natural mortality", file="AS_DET.dat", sep="\n", append=TRUE)
cat(sim_det$M, file="AS_DET.dat", sep="\n", append=TRUE)	

## vectors - use sep=" "
cat("# Weight-at-age", file="AS_DET.dat", sep="\n", append=TRUE)
cat(sim_det$W_a, file="AS_DET.dat", sep="\t", append=TRUE)
cat("\n", file="AS_DET.dat", append=TRUE)	

## back to numbers and integers
cat("# logR0", file="AS_DET.dat", sep="\n", append=TRUE)
cat(log(sim_det$R0), file="AS_DET.dat", sep="\n", append=TRUE)

## use errors to simulate stochastic data within ADMB
cat("# SigmaR", file="AS_DET.dat", sep="\n", append=TRUE)
cat(0.7, file="AS_DET.dat", sep="\n", append=TRUE)
cat("# SigmaF", file="AS_DET.dat", sep="\n", append=TRUE)
cat(0.3, file="AS_DET.dat", sep="\n", append=TRUE)
cat("# CV (Catch)", file="AS_DET.dat", sep="\n", append=TRUE)
cat(0.05, file="AS_DET.dat", sep="\n", append=TRUE)
cat("# Sigma (log-CPUE)", file="AS_DET.dat", sep="\n", append=TRUE)
cat(0.2, file="AS_DET.dat", sep="\n", append=TRUE)
cat("# Omega (proportion data)", file="AS_DET.dat", sep="\n", append=TRUE)
cat(50, file="AS_DET.dat", sep="\n", append=TRUE)	

## matrix - save by row
cat("# Year Catch CPUE", file="AS_DET.dat", sep="\n", append=TRUE)
cmat <- matrix(c(1:sim_det$Nyear, sim_det$C_t, sim_det$CPUE_t), nrow=sim_det$Nyear, ncol=3)
for(i in 1:nrow(cmat)){
	cat(cmat[i,], file="AS_DET.dat", sep="\t", append=TRUE)
	cat("\n", file="AS_DET.dat", append=TRUE)
}
cat("# catch-at-age (year=rows, age=column)", file="AS_DET.dat", sep="\n", append=TRUE)
pmat <- t(sim_det$page)
for(i in 1:nrow(pmat)){
	cat(pmat[i,], file="AS_DET.dat", sep="\t", append=TRUE)
	cat("\n", file="AS_DET.dat", append=TRUE)
}	
cat("# end of file", file="AS_DET.dat", sep="\n", append=TRUE)
cat(999, file="AS_DET.dat", sep="\n", append=TRUE)

################################
## create ADMB pin file
################################
## parameter starting values
cat("# Age-structured simulation pin file", file="AS_DET.pin", sep="\n")
cat("# logR", file="AS_DET.PIN", sep="\n", append=TRUE)
cat(log(sim_det$R_t), file="AS_DET.PIN", sep="\t", append=TRUE)
cat("\n", file="AS_DET.PIN", append=TRUE)
cat("# Sel50", file="AS_DET.PIN", sep="\n", append=TRUE)
cat(sim_det$S50, file="AS_DET.PIN", sep="\n", append=TRUE)
cat("# Sel95", file="AS_DET.PIN", sep="\n", append=TRUE)
cat(sim_det$S95, file="AS_DET.PIN", sep="\n", append=TRUE)
cat("# logF", file="AS_DET.PIN", sep="\n", append=TRUE)
cat(log(sim_det$F_t), file="AS_DET.PIN", sep="\t", append=TRUE)
cat("\n", file="AS_DET.PIN", append=TRUE)
cat("# logq", file="AS_DET.PIN", sep="\n", append=TRUE)
cat(log(sim_det$q), file="AS_DET.PIN", sep="\n", append=TRUE)

## parameter starting values
cat("# Age-structured simulation pin file", file="AS_default.pin", sep="\n")
cat("# logR", file="AS_default.pin", sep="\n", append=TRUE)
cat(rep(log(1000),sim_det$Nyear), file="AS_default.pin", sep="\t", append=TRUE)
cat("\n", file="AS_default.pin", append=TRUE)
cat("# Sel50", file="AS_default.pin", sep="\n", append=TRUE)
cat(3, file="AS_default.pin", sep="\n", append=TRUE)
cat("# Sel95", file="AS_default.pin", sep="\n", append=TRUE)
cat(4, file="AS_default.pin", sep="\n", append=TRUE)
cat("# logF", file="AS_default.pin", sep="\n", append=TRUE)
cat(rep(log(1),sim_det$Nyear), file="AS_default.pin", sep="\t", append=TRUE)
cat("\n", file="AS_default.pin", append=TRUE)
cat("# logq", file="AS_default.pin", sep="\n", append=TRUE)
cat(log(1e-2), file="AS_default.pin", sep="\n", append=TRUE)

################################
## create ADMB ctl file
################################
## control options so we don't have to re-compile each time we want to run a different setup
## estimate all parameters
cat("# Age-structured simulation ctl file", file="AS_EST.ctl", sep="\n")
cat("# R phase", file="AS_EST.ctl", sep="\n", append=TRUE)
cat(1, file="AS_EST.ctl", sep="\n", append=TRUE)
cat("# Sel50 phase", file="AS_EST.ctl", sep="\n", append=TRUE)
cat(1, file="AS_EST.ctl", sep="\n", append=TRUE)
cat("# Sel95 phase", file="AS_EST.ctl", sep="\n", append=TRUE)
cat(1, file="AS_EST.ctl", sep="\n", append=TRUE)
cat("# F phase", file="AS_EST.ctl", sep="\n", append=TRUE)
cat(1, file="AS_EST.ctl", sep="\n", append=TRUE)
cat("# q phase", file="AS_EST.ctl", sep="\n", append=TRUE)
cat(1, file="AS_EST.ctl", sep="\n", append=TRUE)
cat("# end of file", file="AS_EST.ctl", sep="\n", append=TRUE)
cat(999, file="AS_EST.ctl", sep="\n", append=TRUE)

#########################################################################
## Model run 1: deterministic data file, estimation mode, no simulation
#########################################################################
setwd(proj_dir)
dir1 <- file.path(proj_dir, "run1")
dir.create("run1", showWarnings=FALSE)
setwd(dir1)

## copy TPL to run directory
file.copy(from=tpl, to=dir1, overwrite=TRUE)

## copy exe
file.copy(from=file.path(proj_dir,exe), to=dir1, overwrite=TRUE)

## choose data file
det_data <- file.path(proj_dir, "AS_DET.dat")
file.copy(from=det_data, to=dir1, overwrite=TRUE)

## choose pin file -- default values - away from true values to make sure ADMB is estimating properly
det_pin <- file.path(proj_dir, "AS_DET.pin")
def_pin <- file.path(proj_dir, "AS_default.pin")
file.copy(from=def_pin, to=dir1, overwrite=TRUE)
file.rename(from=file.path(dir1, "AS_default.pin"), to="AS_sim.pin")

## choose control file
est_ctl <- file.path(proj_dir, "AS_EST.ctl")
file.copy(from=est_ctl, to=dir1, overwrite=TRUE)


	#######################################################
	## create ADMB .dat file to ID data and control files
	#######################################################
	cat("# Age-structured simulation files", file="AS_sim.dat", sep="\n")
	cat("# Deterministic data file", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("AS_DET.dat", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("# control file calling to estimate all parameters", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("AS_EST.ctl", file="AS_sim.dat", sep="\n", append=TRUE)

	#######################################################
	## run
	#######################################################
	run_admb(name)

	rep <- readFit("AS_sim")

	## input files
	F_inp <- exp(readVec("# logF", "AS_sim.pin"))
	R_inp <- exp(readVec("# logR", "AS_sim.pin"))
	data_mat  <- readMat("# Year Catch CPUE", "AS_DET.dat", nrow=sim_det$Nyear)
	Catch_inp <- data_mat[,2]
	CPUE_inp <- data_mat[,3]
	S50_inp <- readVec("# Sel50", "AS_sim.pin")
	S95_inp <- readVec("# Sel95", "AS_sim.pin")
	S_inp <- 1/ (1 + exp(-log(19) * ((sim_det$ages - 1) - S50_inp)/(S95_inp - S50_inp)))
	page_inp <- readMat("# catch-at-age (year=rows, age=column)", "AS_DET.dat", nrow=sim_det$Nyear)


	## read true values from simulation file
	F_sim <- readVec("F", "AS_sim.sim")
	R_sim <- readVec("R", "AS_sim.sim")	
	VB_sim <- readVec("VulBio", "AS_sim.sim")
	S_sim <- readVec("Selex", "AS_sim.sim")
	N_sim <- readVec("Nt", "AS_sim.sim")
	logF_devs <- readVec("log_F_devs", "AS_sim.sim")
	logR_devs <- readVec("log_R_devs", "AS_sim.sim")
	Catch_sim <- readVec("Catch_obs", "AS_sim.sim")
	CPUE_sim <- readVec("CPUE_obs", "AS_sim.sim")
	page_sim <- readMat("Propn_obs", "AS_sim.sim", nrow=sim_det$Nyear)

	## read output values from estimation model
	F_est <- readVec("F", "AS_sim.rep")
	F_sd <- rep$std[which(rep$names=="logF")]
	R_est <- readVec("Recruits", "AS_sim.rep")	
	R_sd <- rep$std[which(rep$names=="logR")]
	VB_est <- readVec("VulBio", "AS_sim.rep")
	N_est <- readVec("Nt", "AS_sim.rep")
	S_est <- readVec("Selex", "AS_sim.rep")
	Weight <- readVec("Weight", "AS_sim.rep")
	Catch_obs <- readVec("Catch_obs", "AS_sim.rep")
	Catch_pred <- readVec("Catch_pred", "AS_sim.rep")
	CPUE_obs <- readVec("CPUE_obs", "AS_sim.rep")
	page_obs <- readMat("Propn_obs", "AS_sim.rep", nrow=sim_det$Nyear)
	CPUE_pred <- readVec("CPUE_pred", "AS_sim.rep")
	page_pred <- readMat("Propn_pred", "AS_sim.rep", nrow=sim_det$Nyear)

	par(mfrow=c(2,3))
	plot(sim_det$F_t, ylim=c(0,1), col="red", cex=1.5)
	points(F_inp)
	lines(F_est, col="red")
	lines(F_sim, lty=2)
	polygon(x=c(1:length(F_est), length(F_est):1), y=c(F_est - 1.96*F_sd, rev(F_est + 1.96*F_sd)), col="#AA000040", border=NA)

	plot(sim_det$R_t, ylim=c(0, max(c(sim_det$R_t, R_sim, R_est))), col="red", cex=1.5)
	points(R_inp)
	lines(R_est, col="red")
	lines(R_sim, lty=2)
	polygon(x=c(1:length(R_est), length(R_est):1), y=c(R_est - 1.96*R_sd, rev(R_est + 1.96*R_sd)), col="#AA000040", border=NA)

	# plot(sim_det$VB_t, ylim=c(0, max(c(sim_det$VB_t, VB_sim, VB_est))), col="red", cex=1.5)
	# lines(VB_est, col="red")
	# lines(VB_sim, lty=2)

	plot(sim_det$N_t, ylim=c(0, max(c(sim_det$N_t, N_sim, N_est))), col="red", cex=1.5)
	lines(N_est, col="red")
	lines(N_sim, lty=2)

	plot(sim_det$C_t, ylim=c(0, max(c(sim_det$C_t, Catch_obs, Catch_pred))), col="red", cex=1.5)
	points(Catch_inp)
	lines(Catch_pred, col="red")
	lines(Catch_sim, lty=2)
	points(Catch_obs, pch=19)

	plot(sim_det$CPUE_t, ylim=c(0, max(c(sim_det$CPUE_t, CPUE_obs, CPUE_pred))), col="red", cex=1.5)
	points(CPUE_inp)
	lines(CPUE_pred, col="red")
	lines(CPUE_sim, lty=2)
	points(CPUE_obs, pch=19)

	plot(sim_det$page[,20], col="red", cex=1.5, ylim=c(0, max(c(sim_det$page[,20], page_inp[20,], page_pred[20,], page_sim[20,], page_obs[20,]))))
	points(page_inp[20,])
	lines(page_pred[20,], col="red")
	lines(page_sim[20,], lty=2)
	points(page_obs[20,], pch=19)

	# plot(sim_det$S_a, ylim=c(0,1), col="red", cex=1.5)
	# points(S_inp)
	# lines(S_est, col="red")
	# lines(S_sim, lty=2)
	legend("topright", legend=c("Simulated in R", "Input to ADMB", "Estimated", "ADMB simulation", "Observed data"), pch=c(1,1,NA,NA,19), lty=c(NA,NA,1,2,NA), col=c("red", "black", "red", "black","black"))

### can estimate simulated values from R
### ADMB only sees PIN file with default values (black line)

##########################################################################################################
## Model run 2: deterministic data file, estimation mode, simulation mode (ADD STOCHASTICITY WITHIN ADMB)
##########################################################################################################
setwd(proj_dir)
dir2 <- file.path(proj_dir, "run2")
dir.create("run2", showWarnings=FALSE)
setwd(dir2)

## copy TPL to run directory
file.copy(from=tpl, to=dir2, overwrite=TRUE)

## copy exe
file.copy(from=file.path(proj_dir,exe), to=dir2, overwrite=TRUE)

## choose data file
det_data <- file.path(proj_dir, "AS_DET.dat")
file.copy(from=det_data, to=dir2, overwrite=TRUE)

## choose pin file -- deterministic model - make sure starting values are at truth for adding error
det_pin <- file.path(proj_dir, "AS_DET.pin")
def_pin <- file.path(proj_dir, "AS_default.pin")
# file.copy(from=def_pin, to=dir2, overwrite=TRUE)
# file.rename(from=file.path(dir2, "AS_default.pin"), to="AS_sim.pin")
file.copy(from=det_pin, to=dir2, overwrite=TRUE)
file.rename(from=file.path(dir2, "AS_DET.pin"), to="AS_sim.pin")

## choose control file
est_ctl <- file.path(proj_dir, "AS_EST.ctl")
file.copy(from=est_ctl, to=dir2, overwrite=TRUE)


	#######################################################
	## create ADMB .dat file to ID data and control files
	#######################################################
	cat("# Age-structured simulation files", file="AS_sim.dat", sep="\n")
	cat("# Deterministic data file", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("AS_DET.dat", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("# control file calling to estimate all parameters", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("AS_EST.ctl", file="AS_sim.dat", sep="\n", append=TRUE)

	#######################################################
	## run
	#######################################################
	run_admb(name, extra.args="-sim 123")

	rep <- readFit("AS_sim")

	## input files
	F_inp <- exp(readVec("# logF", "AS_sim.pin"))
	R_inp <- exp(readVec("# logR", "AS_sim.pin"))
	data_mat  <- readMat("# Year Catch CPUE", "AS_DET.dat", nrow=sim_det$Nyear)
	Catch_inp <- data_mat[,2]
	CPUE_inp <- data_mat[,3]
	S50_inp <- readVec("# Sel50", "AS_sim.pin")
	S95_inp <- readVec("# Sel95", "AS_sim.pin")
	S_inp <- 1/ (1 + exp(-log(19) * ((sim_det$ages - 1) - S50_inp)/(S95_inp - S50_inp)))
	page_inp <- readMat("# catch-at-age (year=rows, age=column)", "AS_DET.dat", nrow=sim_det$Nyear)


	## read true values from simulation file
	F_sim <- readVec("F", "AS_sim.sim")
	R_sim <- readVec("R", "AS_sim.sim")	
	VB_sim <- readVec("VulBio", "AS_sim.sim")
	S_sim <- readVec("Selex", "AS_sim.sim")
	N_sim <- readVec("Nt", "AS_sim.sim")
	logF_devs <- readVec("log_F_devs", "AS_sim.sim")
	logR_devs <- readVec("log_R_devs", "AS_sim.sim")
	Catch_sim <- readVec("Catch_obs", "AS_sim.sim")
	CPUE_sim <- readVec("CPUE_obs", "AS_sim.sim")
	page_sim <- readMat("Propn_obs", "AS_sim.sim", nrow=sim_det$Nyear)

	## read output values from estimation model
	F_est <- readVec("F", "AS_sim.rep")
	F_sd <- rep$std[which(rep$names=="logF")]
	R_est <- readVec("Recruits", "AS_sim.rep")	
	R_sd <- rep$std[which(rep$names=="logR")]
	VB_est <- readVec("VulBio", "AS_sim.rep")
	N_est <- readVec("Nt", "AS_sim.rep")
	S_est <- readVec("Selex", "AS_sim.rep")
	Weight <- readVec("Weight", "AS_sim.rep")
	Catch_obs <- readVec("Catch_obs", "AS_sim.rep")
	Catch_pred <- readVec("Catch_pred", "AS_sim.rep")
	CPUE_obs <- readVec("CPUE_obs", "AS_sim.rep")
	page_obs <- readMat("Propn_obs", "AS_sim.rep", nrow=sim_det$Nyear)
	CPUE_pred <- readVec("CPUE_pred", "AS_sim.rep")
	page_pred <- readMat("Propn_pred", "AS_sim.rep", nrow=sim_det$Nyear)

	par(mfrow=c(2,3))
	plot(sim_det$F_t, ylim=c(0,1), col="red", cex=1.5)
	points(F_inp)
	lines(F_est, col="red")
	lines(F_sim, lty=2)
	polygon(x=c(1:length(F_est), length(F_est):1), y=c(F_est - 1.96*F_sd, rev(F_est + 1.96*F_sd)), col="#AA000040", border=NA)

	plot(sim_det$R_t, ylim=c(0, max(c(sim_det$R_t, R_sim, R_est))), col="red", cex=1.5)
	points(R_inp)
	lines(R_est, col="red")
	lines(R_sim, lty=2)
	polygon(x=c(1:length(R_est), length(R_est):1), y=c(R_est - 1.96*R_sd, rev(R_est + 1.96*R_sd)), col="#AA000040", border=NA)

	# plot(sim_det$VB_t, ylim=c(0, max(c(sim_det$VB_t, VB_sim, VB_est))), col="red", cex=1.5)
	# lines(VB_est, col="red")
	# lines(VB_sim, lty=2)

	plot(sim_det$N_t, ylim=c(0, max(c(sim_det$N_t, N_sim, N_est))), col="red", cex=1.5)
	lines(N_est, col="red")
	lines(N_sim, lty=2)

	plot(sim_det$C_t, ylim=c(0, max(c(sim_det$C_t, Catch_obs, Catch_pred))), col="red", cex=1.5)
	points(Catch_inp)
	lines(Catch_pred, col="red")
	lines(Catch_sim, lty=2)
	points(Catch_obs, pch=19)

	plot(sim_det$CPUE_t, ylim=c(0, max(c(sim_det$CPUE_t, CPUE_obs, CPUE_pred))), col="red", cex=1.5)
	points(CPUE_inp)
	lines(CPUE_pred, col="red")
	lines(CPUE_sim, lty=2)
	points(CPUE_obs, pch=19)

	plot(sim_det$page[,20], col="red", cex=1.5, ylim=c(0, max(c(sim_det$page[,20], page_inp[20,], page_pred[20,], page_sim[20,], page_obs[20,]))))
	points(page_inp[20,])
	lines(page_pred[20,], col="red")
	lines(page_sim[20,], lty=2)
	points(page_obs[20,], pch=19)

	# plot(sim_det$S_a, ylim=c(0,1), col="red", cex=1.5)
	# points(S_inp)
	# lines(S_est, col="red")
	# lines(S_sim, lty=2)
	legend("topright", legend=c("Simulated in R", "Input to ADMB", "Estimated", "ADMB simulation", "Observed data"), pch=c(1,1,NA,NA,19), lty=c(NA,NA,1,2,NA), col=c("red", "black", "red", "black","black"))


##########################################################################################################
## Model run 3: loop over stochastic runs within ADMB
##########################################################################################################
setwd(proj_dir)

itervec <- 1:10
dir3 <- file.path(proj_dir, "run_iters")
dir.create("run_iters", showWarnings=FALSE)

for(i in 1:length(itervec)){
	iterdir <- file.path(dir3, itervec[i])
	dir.create(iterdir, showWarnings=FALSE)
	setwd(iterdir)

	## copy TPL to run directory
	file.copy(from=tpl, to=iterdir, overwrite=TRUE)	

	## copy exe
	file.copy(from=file.path(proj_dir,exe), to=iterdir, overwrite=TRUE)	

	## choose data file
	det_data <- file.path(proj_dir, "AS_DET.dat")
	file.copy(from=det_data, to=iterdir, overwrite=TRUE)	

	## choose pin file -- deterministic model - make sure starting values are at truth for adding error
	det_pin <- file.path(proj_dir, "AS_DET.pin")
	def_pin <- file.path(proj_dir, "AS_default.pin")
	# file.copy(from=def_pin, to=iterdir, overwrite=TRUE)
	# file.rename(from=file.path(iterdir, "AS_default.pin"), to="AS_sim.pin")
	file.copy(from=det_pin, to=iterdir, overwrite=TRUE)
	file.rename(from=file.path(iterdir, "AS_DET.pin"), to="AS_sim.pin")	

	## choose control file
	est_ctl <- file.path(proj_dir, "AS_EST.ctl")
	file.copy(from=est_ctl, to=iterdir, overwrite=TRUE)

	#######################################################
	## create ADMB .dat file to ID data and control files
	#######################################################
	cat("# Age-structured simulation files", file="AS_sim.dat", sep="\n")
	cat("# Deterministic data file", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("AS_DET.dat", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("# control file calling to estimate all parameters", file="AS_sim.dat", sep="\n", append=TRUE)
	cat("AS_EST.ctl", file="AS_sim.dat", sep="\n", append=TRUE)

	#######################################################
	## run
	#######################################################
	## set seed at itervec value
	run_admb(name, extra.args=paste("-sim", itervec[i]))

}

	## input files
	F_inp <- exp(readVec("# logF", "AS_sim.pin"))
	R_inp <- exp(readVec("# logR", "AS_sim.pin"))
	data_mat  <- readMat("# Year Catch CPUE", "AS_DET.dat", nrow=sim_det$Nyear)
	Catch_inp <- data_mat[,2]
	CPUE_inp <- data_mat[,3]
	S50_inp <- readVec("# Sel50", "AS_sim.pin")
	S95_inp <- readVec("# Sel95", "AS_sim.pin")
	S_inp <- 1/ (1 + exp(-log(19) * ((sim_det$ages - 1) - S50_inp)/(S95_inp - S50_inp)))
	page_inp <- readMat("# catch-at-age (year=rows, age=column)", "AS_DET.dat", nrow=sim_det$Nyear)


	## read true values from simulation file
	F_sim <- readVec("F", "AS_sim.sim")
	R_sim <- readVec("R", "AS_sim.sim")	
	VB_sim <- readVec("VulBio", "AS_sim.sim")
	S_sim <- readVec("Selex", "AS_sim.sim")
	N_sim <- readVec("Nt", "AS_sim.sim")
	logF_devs <- readVec("log_F_devs", "AS_sim.sim")
	logR_devs <- readVec("log_R_devs", "AS_sim.sim")
	Catch_sim <- readVec("Catch_obs", "AS_sim.sim")
	CPUE_sim <- readVec("CPUE_obs", "AS_sim.sim")
	page_sim <- readMat("Propn_obs", "AS_sim.sim", nrow=sim_det$Nyear)

	## read output values from estimation model
	F_est <- readVec("F", "AS_sim.rep")
	F_sd <- rep$std[which(rep$names=="logF")]
	R_est <- readVec("Recruits", "AS_sim.rep")	
	R_sd <- rep$std[which(rep$names=="logR")]
	VB_est <- readVec("VulBio", "AS_sim.rep")
	N_est <- readVec("Nt", "AS_sim.rep")
	S_est <- readVec("Selex", "AS_sim.rep")
	Weight <- readVec("Weight", "AS_sim.rep")
	Catch_obs <- readVec("Catch_obs", "AS_sim.rep")
	Catch_pred <- readVec("Catch_pred", "AS_sim.rep")
	CPUE_obs <- readVec("CPUE_obs", "AS_sim.rep")
	page_obs <- readMat("Propn_obs", "AS_sim.rep", nrow=sim_det$Nyear)
	CPUE_pred <- readVec("CPUE_pred", "AS_sim.rep")
	page_pred <- readMat("Propn_pred", "AS_sim.rep", nrow=sim_det$Nyear)

	par(mfrow=c(2,3))
	plot(sim_det$F_t, ylim=c(0,1), col="red", cex=1.5)
	points(F_inp)
	lines(F_est, col="red")
	lines(F_sim, lty=2)
	polygon(x=c(1:length(F_est), length(F_est):1), y=c(F_est - 1.96*F_sd, rev(F_est + 1.96*F_sd)), col="#AA000040", border=NA)

	plot(sim_det$R_t, ylim=c(0, max(c(sim_det$R_t, R_sim, R_est))), col="red", cex=1.5)
	points(R_inp)
	lines(R_est, col="red")
	lines(R_sim, lty=2)
	polygon(x=c(1:length(R_est), length(R_est):1), y=c(R_est - 1.96*R_sd, rev(R_est + 1.96*R_sd)), col="#AA000040", border=NA)

	# plot(sim_det$VB_t, ylim=c(0, max(c(sim_det$VB_t, VB_sim, VB_est))), col="red", cex=1.5)
	# lines(VB_est, col="red")
	# lines(VB_sim, lty=2)

	plot(sim_det$N_t, ylim=c(0, max(c(sim_det$N_t, N_sim, N_est))), col="red", cex=1.5)
	lines(N_est, col="red")
	lines(N_sim, lty=2)

	plot(sim_det$C_t, ylim=c(0, max(c(sim_det$C_t, Catch_obs, Catch_pred))), col="red", cex=1.5)
	points(Catch_inp)
	lines(Catch_pred, col="red")
	lines(Catch_sim, lty=2)
	points(Catch_obs, pch=19)

	plot(sim_det$CPUE_t, ylim=c(0, max(c(sim_det$CPUE_t, CPUE_obs, CPUE_pred))), col="red", cex=1.5)
	points(CPUE_inp)
	lines(CPUE_pred, col="red")
	lines(CPUE_sim, lty=2)
	points(CPUE_obs, pch=19)

	plot(sim_det$page[,20], col="red", cex=1.5, ylim=c(0, max(c(sim_det$page[,20], page_inp[20,], page_pred[20,], page_sim[20,], page_obs[20,]))))
	points(page_inp[20,])
	lines(page_pred[20,], col="red")
	lines(page_sim[20,], lty=2)
	points(page_obs[20,], pch=19)

	# plot(sim_det$S_a, ylim=c(0,1), col="red", cex=1.5)
	# points(S_inp)
	# lines(S_est, col="red")
	# lines(S_sim, lty=2)
	legend("topright", legend=c("Simulated in R", "Input to ADMB", "Estimated", "ADMB simulation", "Observed data"), pch=c(1,1,NA,NA,19), lty=c(NA,NA,1,2,NA), col=c("red", "black", "red", "black","black"))
