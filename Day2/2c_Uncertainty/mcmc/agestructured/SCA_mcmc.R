wd <- "C:\\merrill\\ADMB_course"

project_dir <- file.path(wd, "Day2", "2c_Uncertainty", "mcmc", "agestructured")
setwd(project_dir)

library(R2admb)

compile_admb("SCA")

nsamp <- 10000
nsave <- 10
run_admb("SCA", extra.args=paste("-mcmc", nsamp, "-mcsave", nsave))

run_admb("SCA", extra.args="-mceval")

chains <- read.table("refpar.mcmc", header=TRUE)

par(mfrow=c(3,2))
for(i in 1:ncol(chains)){
	acf(chains[,i])
}

par(mfrow=c(3,2))
for(i in 1:ncol(chains)){
	plot(chains[,1], xlab="Sample", ylab=colnames(chains)[i])
}

nburn <- 0.5
chain_new <- chains[-c(1:(nburn*nrow(chains))),]
par(mfrow=c(3,2))
for(i in 1:ncol(chain_new)){
	acf(chain_new[,i])
}

par(mfrow=c(3,2))
for(i in 1:ncol(chain_new)){
	plot(chain_new[,1], xlab="Sample", ylab=colnames(chains)[i])
}