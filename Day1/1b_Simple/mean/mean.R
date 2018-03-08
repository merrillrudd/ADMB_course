wd <- "C:\\merrill\\ADMB_course"
ex_wd <- file.path(wd, "example solutions", "mean")
setwd(ex_wd)
library(RColorBrewer)

source(file.path(wd, "tools.R"))

## read observed and estimated average
read_vars <- c("Observed", "Average")

## read results by model
models <- c("rss", "like_full", "like_conc")

## find results from rep file
results <- lapply(1:length(models), function(x){
	estimates <- lapply(1:length(read_vars), function(y){
		readVec(read_vars[y], paste0("mean_", models[x], ".rep"))
	})
	names(estimates) <- read_vars
	return(estimates)
})
names(results) <- models

colors <- brewer.pal(length(models), "Set1")
boxplot(results[[1]]$Observed)
for(m in 1:length(models)){
	abline(h=results[[m]]$Average, col=colors[m], lty=m, lwd=5)
}