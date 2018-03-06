rec <- read.table("ricker2.dat", skip=4, col.names=c("S","R"))

fm <- lm(log(R/S)~S, data=rec)

a <- exp(coef(fm)[[1]])
b <- -coef(fm)[[2]]
