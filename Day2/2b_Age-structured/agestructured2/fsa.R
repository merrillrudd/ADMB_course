################################################################################
###                                                                            #
### Script:    fsa.R                                                           #
###                                                                            #
### Purpose:   Plot and discuss simple statistical catch-at-age model          #
###                                                                            #
### Functions: readVec, readMat, importFSA                                     #
###                                                                            #
### Requires:  package:scape, fsa.dat, fsa.par                                 #
###                                                                            #
### Notes:     fsa.dat and fsa.par are ADMB input/output files                 #
###                                                                            #
### History:   2013-20-02 Arni Magnusson created for ICES ADMB workshop        #
###                                                                            #
################################################################################

readVec <- function(string, file)
### Find 'string' in 'file' and read vector from next line
{
  txt <- readLines(file)
  skip <- match(string, txt)
  vec <- scan(file, quiet=TRUE, skip=skip, nlines=1)
  return(vec)
}

readMat <- function(string, file, nrow)
### Find 'string' in 'file' and read matrix with 'nrow' rows from next line
{
  txt <- readLines(file)
  skip <- match(string, txt)
  mat <- as.matrix(read.table(file, skip=skip, nrows=nrow))
  dimnames(mat) <- NULL
  return(mat)
}

require(scape)
importFSA <- function(run="C:\\merrill\\ADMB_course\\Day2\\2b_Age-structured\\agestructured2\\fsa", info="")
### Import model data and fit
{
  ## 1a  Import model dimensions
  datfile <- paste0(run, ".dat")
  minAge <- readVec("# minAge", datfile)
  maxAge <- readVec("# maxAge", datfile)
  minYear <- readVec("# minYear", datfile)
  maxYear <- readVec("# maxYear", datfile)
  ages <- minAge:maxAge
  nages <- length(ages)
  years <- minYear:maxYear
  nyears <- length(years)

  ## 1b  Import survey dimensions
  minAgeS <- readVec("# minAgeS", datfile)
  maxAgeS <- readVec("# maxAgeS", datfile)
  minYearS <- readVec("# minYearS", datfile)
  maxYearS <- readVec("# maxYearS", datfile)
  ages.s <- minAgeS:maxAgeS
  nages.s <- length(ages.s)
  years.s <- minYearS:maxYearS
  nyears.s <- length(years.s)

  ## 1c  Import data
  catch <- readMat("# catch in numbers", datfile, nages)
  weight <- readMat("# stock mean weight", datfile, nages)
  maturity <- readMat("# prop mature", datfile, nages)
  natmort <- readMat("# Assumed M", datfile, nages)
  surveyTime <- readVec("# survey time (fraction into year)", datfile)
  survey <- readMat("# Q1 survey", datfile, nages.s)

  ## 1d  Import parameter estimates
  parfile <- paste0(run, ".par")
  logN1Y <- readVec("# logN1Y:", parfile)
  logN1A <- readVec("# logN1A:", parfile)
  logFY <- readVec("# logFY:", parfile)
  logFA <- readVec("# logFA:", parfile)
  logVarC <- readVec("# logVarLogCatch:", parfile)
  logQ <- readVec("# logQ:", parfile)
  logVarS <- readVec("# logVarLogSurvey:", parfile)

  ## 2a  Calculate F, M, Z
  Fa <- c(exp(logFA), 1, 1, 1)
  names(Fa) <- ages
  Ft <- exp(logFY)
  names(Ft) <- years
  Fmat <- outer(Ft, Fa)
  M <- t(natmort)
  dimnames(M) <- dimnames(Fmat)
  Z <- Fmat+M

  ## 2b  Calculate Fbar24, Sel
  Fbar24 <- rowMeans(Fmat[,2:4])
  SelC <- c(exp(logFA), 1, 1, 1)
  SelC <- SelC / max(SelC)
  names(SelC) <- ages
  Q <- exp(logQ)
  SelS <- Q / max(Q)

  ## 2c  Calculate N
  N <- matrix(0, nrow=nyears, ncol=nages, dimnames=dimnames(Fmat))
  N[1,] <- exp(logN1Y)
  N[-1,1] <- exp(logN1A)
  for(t in 1:(nyears-1))
  {
    N[t+1,-1] <- N[t,-maxAge] * exp(-Z[t,-maxAge])
  }

  ## 2d  Calculcate SB, VB, R
  mat <- t(maturity)
  dimnames(mat) <- dimnames(N)
  w <- t(weight)
  dimnames(w) <- dimnames(N)
  SB <- rowSums(N * w * mat)               # spawning biomass (mature)
  VB <- rowSums(sweep(N*w, 2, SelC, "*"))  # vulnerable biomass (fleet selected)
  R <- c(N[-1,1], NA)                      # recruitment (by birth year)

  ## 2e  Calculate CA, q, sigma
  CAc.obs <- t(catch)
  dimnames(CAc.obs) <- list(years, ages)
  CAc.fit <- Fmat/Z * N*(1-exp(-Z))
  sigmaC <- sqrt(exp(logVarC))
  CAs.obs <- t(survey)
  dimnames(CAs.obs) <- list(years.s, ages.s)
  sigmaS <- sqrt(exp(logVarS))
  Nsurvey <- N[as.character(years.s),ages.s]
  Zsurvey <- Z[as.character(years.s),ages.s]
  CAs.fit <- sweep(Nsurvey * exp(-Zsurvey*surveyTime), 2, Q, "*")

  ## 3  Construct scape
  N.out <- data.frame(Sex="Unisex", Year=rep(years,each=nages),
                      Age=rep(ages,nyears), N=c(t(N)))
  B.out <- data.frame(Year=years, VB=VB, SB=SB, Y=NA_real_, R=R)
  Sel.out <- data.frame(Series=rep(c("Fleet","Survey","Maturity"),
                          c(nages,nages.s,nages)), Sex="Unisex",
                        Age=c(ages,ages.s,ages), P=c(SelC,SelS,colMeans(mat)))
  div <- 1e3
  CAc.out <- data.frame(Series="Fleet", Year=rep(years,each=nages), SS=sigmaC,
                        Sex="Unisex", Age=rep(ages,nyears),
                        Obs=c(t(CAc.obs))/div, Fit=c(t(CAc.fit))/div)
  CAs.out <- data.frame(Series="Survey", Year=rep(years.s,each=nages.s),
                        SS=sigmaS, Sex="Unisex", Age=rep(ages.s,nyears.s),
                        Obs=c(t(CAs.obs)), Fit=c(t(CAs.fit)))
  model <- list(N=N.out, B=B.out, Sel=Sel.out, CAc=CAc.out, CAs=CAs.out)
  attr(model,"call") <- match.call()
  attr(model,"scape.version") <- installed.packages()["scape","Version"]
  attr(model,"info") <- info
  class(model) <- "scape"

  return(model)
}

## Import and plot

model <- importFSA()

plotSel(model, together=T, strip=F, main="\nSelectivity and maturity")

plotCA(model, fit=F)
plotCA(model, main="\nCommercial CA")
plotCA(model, same=F, main="\nCommercial CA")
plotCA(model, log=T, main="\nCommercial CA")
plotCA(model, log=T, swap=T, main="\nCommercial CA")

plotCA(model, "s", fit=F)
plotCA(model, "s", log=T, main="\nSurvey CA")
plotCA(model, "s", log=T, swap=T, main="\nSurvey CA")

plotN(model, div=1e3)
plotN(model, "b", xlab="Age", cex.points=.6, main="\nPopulation numbers at age")
plotN(x.cod, "b", main="\nIcelandic cod")

plotB(model, div=1e3, main="\nBiomass", ylab="kt")

round(model$N[model$N$Year==2007,]$N/1000, 1)
