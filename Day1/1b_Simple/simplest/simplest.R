main_dir <- "C:\\merrill\\ADMB_course"
proj_dir <- file.path(main_dir, "Day1", "1b_Simple", "simplest")
setwd(proj_dir)

source(file.path(main_dir, "tools.R"))

## use readVec function in tools.R file
## read from report file
obs <- readVec("Observed", "simplest.rep")
pred <- readVec("Predicted", "simplest.rep")

plot(obs, pch=19, cex=1.5, ylim=c(0, max(obs)*1.1), xlab="Observation", ylab="Value")
lines(pred, col="red", lwd=2)
