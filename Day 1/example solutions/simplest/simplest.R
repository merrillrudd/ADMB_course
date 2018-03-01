wd <- "C:\\merrill\\ADMB_course\\Day1\\examples\\simplest"
setwd(wd)

library(ggplot2)
library(dplyr)

readVec <- function(string, file)
  ### Find 'string' in 'file' and read vector from next line
{
  txt <- readLines(file)
  skip <- match(string, txt)
  vec <- scan(file, quiet=TRUE, skip=skip, nlines=1)
  return(vec)
}

obs <- readVec("Observed", "simplest.rep")
pred <- readVec("Predicted", "simplest.rep")

plot(obs, pch=19, cex=1.5, ylim=c(0, max(obs)*1.1), xlab="Observation", ylab="Value")
lines(pred, col="red", lwd=2)
