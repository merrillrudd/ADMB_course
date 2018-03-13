rm(list=ls())

# install.packages("R2admb")
library(R2admb)

main_dir <- "C:\\merrill\\ADMB_course"
proj_dir <- file.path(main_dir, "Day3", "3b_Retrospective", "agestructure_retro")
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
sim_stoch <- simPop(linf=60,
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
				SigmaR=0.6,
				SigmaF=0.3,
				Nyear=20,
				InitialDepl=0.8,
				PropFcrash=0.8,
				EffortDyn="Two-way",
				seed=123)

par(mfrow=c(3,3))
plot(sim_stoch$F_t, type="l", lwd=2, ylim=c(0,1), xlab="Year", ylab="Fishing mortality")
plot(sim_stoch$C_t, type="l", lwd=2, ylim=c(0,max(sim_stoch$C_t)*1.1), xlab="Year", ylab="Catch")
plot(sim_stoch$CPUE_t, type="l", lwd=2, ylim=c(0, max(sim_stoch$CPUE_t)*1.1), xlab="Year", ylab="CPUE")
plot(sim_stoch$R_t, type="l", lwd=2, ylim=c(0,3000), xlab="Year", ylab="Recruitment")
plot(sim_stoch$SB_t/sim_stoch$SB0, type="l", lwd=2, ylim=c(0, max(sim_stoch$SB_t/sim_stoch$SB0)*1.1), xlab="Year", ylab="Relative spawning biomass")
plot(sim_stoch$N_t, type="l", lwd=2, ylim=c(0, max(sim_stoch$N_t)*1.1), xlab="Year", ylab="Abundance")
plot(sim_stoch$page[,1], type="h", lwd=4, ylim=c(0, max(sim_stoch$page)), col="gray", xlab="Age", ylab="Proportion in catch")
	lines(sim_stoch$page[,ncol(sim_stoch$page)], type="h", lwd=4)
	legend('topright', legend=c("First year", "Last year"), col=c("gray","black"), lwd=4)
plot(sim_stoch$L_a, type="l", lwd=2, ylim=c(0, max(sim_stoch$L_a)*1.1), xlab="Age", ylab="Length")
plot(sim_stoch$Mat_a, type="l", lwd=2, ylim=c(0,1), xlab="Age", ylab="Proportion")
	lines(sim_stoch$S_a, lwd=2, lty=2)
	legend("bottomright", legend=c("Maturity", "Selectivity"), lty=c(1,2), lwd=2)

#########################################
## create deterministic ADMB data file
#########################################
cat("# Age-structured simulation data file", file="AS_STOCH.dat", sep="\n")	
## numbers and integers - use sep="\n"
cat("# Number of years of data", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(sim_stoch$Nyear, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# Number of years in model", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(sim_stoch$Nyear, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# Number of ages", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(sim_stoch$Nage, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# Natural mortality", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(sim_stoch$M, file="AS_STOCH.dat", sep="\n", append=TRUE)	

## vectors - use sep=" "
cat("# Weight-at-age", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(sim_stoch$W_a, file="AS_STOCH.dat", sep="\t", append=TRUE)
cat("\n", file="AS_STOCH.dat", append=TRUE)	

## back to numbers and integers
cat("# logR0", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(log(sim_stoch$R0), file="AS_STOCH.dat", sep="\n", append=TRUE)

## use errors to simulate stochastic data within ADMB
cat("# SigmaR", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(0.7, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# SigmaF", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(0.3, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# CV (Catch)", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(0.05, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# Sigma (log-CPUE)", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(0.2, file="AS_STOCH.dat", sep="\n", append=TRUE)
cat("# Omega (proportion data)", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(50, file="AS_STOCH.dat", sep="\n", append=TRUE)	

## matrix - save by row
cat("# Year Catch CPUE", file="AS_STOCH.dat", sep="\n", append=TRUE)
cmat <- matrix(c(1:sim_stoch$Nyear, sim_stoch$C_t, sim_stoch$CPUE_t), nrow=sim_stoch$Nyear, ncol=3)
for(i in 1:nrow(cmat)){
	cat(cmat[i,], file="AS_STOCH.dat", sep="\t", append=TRUE)
	cat("\n", file="AS_STOCH.dat", append=TRUE)
}
cat("# catch-at-age (year=rows, age=column)", file="AS_STOCH.dat", sep="\n", append=TRUE)
pmat <- t(sim_stoch$page)
for(i in 1:nrow(pmat)){
	cat(pmat[i,], file="AS_STOCH.dat", sep="\t", append=TRUE)
	cat("\n", file="AS_STOCH.dat", append=TRUE)
}	
cat("# end of file", file="AS_STOCH.dat", sep="\n", append=TRUE)
cat(999, file="AS_STOCH.dat", sep="\n", append=TRUE)

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
## Run retrospectives
#########################################################################

setwd(proj_dir)

retro <- 0:10
rundir <- file.path(proj_dir, "run_retro")
dir.create("run_retro", showWarnings=FALSE)

for(i in 1:length(retro)){
	retrodir <- file.path(rundir, paste0("retro_", retro[i]))
	dir.create(retrodir, showWarnings=FALSE)
	setwd(retrodir)

	## copy TPL to run directory
	file.copy(from=tpl, to=retrodir, overwrite=TRUE)	

	## copy exe
	file.copy(from=file.path(proj_dir,exe), to=retrodir, overwrite=TRUE)	

	## choose data file
	det_data <- file.path(proj_dir, "AS_STOCH.dat")
	file.copy(from=det_data, to=retrodir, overwrite=TRUE)	

	## choose pin file -- deterministic model - make sure starting values are at truth for adding error
	cat("# Age-structured simulation pin file", file="AS_retro.pin", sep="\n")
	cat("# logR", file="AS_retro.pin", sep="\n", append=TRUE)
	cat(rep(log(1000),(sim_stoch$Nyear-retro[i])), file="AS_retro.pin", sep="\t", append=TRUE)
	cat("\n", file="AS_retro.pin", append=TRUE)
	cat("# Sel50", file="AS_retro.pin", sep="\n", append=TRUE)
	cat(3, file="AS_retro.pin", sep="\n", append=TRUE)
	cat("# Sel95", file="AS_retro.pin", sep="\n", append=TRUE)
	cat(4, file="AS_retro.pin", sep="\n", append=TRUE)
	cat("# logF", file="AS_retro.pin", sep="\n", append=TRUE)
	cat(rep(log(1),(sim_stoch$Nyear-retro[i])), file="AS_retro.pin", sep="\t", append=TRUE)
	cat("\n", file="AS_retro.pin", append=TRUE)
	cat("# logq", file="AS_retro.pin", sep="\n", append=TRUE)
	cat(log(1e-2), file="AS_retro.pin", sep="\n", append=TRUE)

	## choose control file
	est_ctl <- file.path(proj_dir, "AS_EST.ctl")
	file.copy(from=est_ctl, to=retrodir, overwrite=TRUE)

	#######################################################
	## create ADMB .dat file to ID data and control files
	#######################################################
	cat("# Age-structured simulation files", file="AS_retro.dat", sep="\n")
	cat("# Deterministic data file", file="AS_retro.dat", sep="\n", append=TRUE)
	cat("AS_STOCH.dat", file="AS_retro.dat", sep="\n", append=TRUE)
	cat("# control file calling to estimate all parameters", file="AS_retro.dat", sep="\n", append=TRUE)
	cat("AS_EST.ctl", file="AS_retro.dat", sep="\n", append=TRUE)

	#######################################################
	## run
	#######################################################
	## set seed at retro value
	run_admb(name, extra.args=paste("-retro", retro[i]))

}

	library(RColorBrewer)
	Fish <- Rec <- Abund <- VB <- Catch <- CPUE <- list()
	for(i in 1:length(retro)){
		retrodir <- file.path(rundir, paste0("retro_", retro[i]))
		setwd(retrodir)

		Fish[[i]] <- readVec("F", "AS_retro.rep")
		Rec[[i]] <- readVec("Recruits", "AS_retro.rep")
		Abund[[i]] <- readVec("Nt", "AS_retro.rep")
		VB[[i]] <- readVec("VulBio", "AS_retro.rep")
	}
	
	col_vec <- brewer.pal(length(retro)-1, "Reds")
	par(mfrow=c(2,2))
	plot(sim_stoch$F_t, ylim=c(0,1), cex=1.5, pch=19)
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(Fish[[i]], col=col)
	}
	plot(sim_stoch$R_t, ylim=c(0,max(c(sim_stoch$R_t))*1.5), cex=1.5, pch=19)
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(Rec[[i]], col=col)
	}
	plot(sim_stoch$N_t, ylim=c(0, max(sim_stoch$N_t)*1.5), cex=1.5, pch=19)
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(Abund[[i]], col=col)
	}
	plot(sim_stoch$VB_t, ylim=c(0, max(sim_stoch$VB_t)*1.5), cex=1.5, pch=19)
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(VB[[i]], col=col)
	}

	