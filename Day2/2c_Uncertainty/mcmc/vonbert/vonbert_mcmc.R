wd <- "C:\\merrill\\ADMB_course"
project_dir <- file.path(wd, "Day2", "2c_Uncertainty", "mcmc", "vonbert")
setwd(project_dir)

library(R2admb)

compile_admb("vonbert1")

nsamp <- 10000
nsave <- 10
run_admb("vonbert1", extra.args=paste("-mcmc", nsamp, "-mcsave", nsave))
run_admb("vonbert1", extra.args="-mceval")

chains <- read.table("refpar.mcmc", header=TRUE)

par(mfrow=c(1,2))
acf(chains[,1])
acf(chains[,2])

par(mfrow=c(1,2))
plot(chains[,1], xlab="Sample", ylab=colnames(chains)[1])
plot(chains[,2], xlab="Sample", ylab=colnames(chains)[2])
