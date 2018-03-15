wd <- "C:\\merrill\\ADMB_course"
project_dir <- file.path(wd, "Day2", "2c_Uncertainty", "likeprof", "vonbert")
setwd(project_dir)

### read in results from report file using code in tools.R
source(file.path(wd, "tools.R"))

## run ADMB in R
library(R2admb)
compile_admb("vonbert1", verbose=TRUE)
run_admb("vonbert1", extra.args="-lprof")

## read in likelihood profile
prof <- readMat(string="Profile likelihood", file="Linf_pro.plt", nrow=87)
profk <- readMat(string="Profile likelihood", file="k_prof.plt", nrow=85)

par(mfrow=c(1,2))
Linf_vec <- prof[,1]
like_vec <- prof[,2]
nll_vec <- -log(like_vec)
plot(Linf_vec, nll_vec)
# plot(Linf_vec, -nll_vec)
# plot(Linf_vec, -exp(nll_vec))

K_vec <- profk[,1]
klike_vec <- profk[,2]
knll_vec <- -log(klike_vec)
plot(K_vec, knll_vec)