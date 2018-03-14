wd <- "C:\\merrill\\ADMB_course"
project_dir <- file.path(wd, "Day2", "2c_Uncertainty", "likeprof", "vonbert")
setwd(project_dir)

### read in results from report file using code in tools.R
source(file.path(wd, "tools.R"))

## read in likelihood profile
prof <- readMat(string="Profile likelihood", file="Linf_pro.plt", nrow=87)

Linf_vec <- prof[,1]
like_vec <- prof[,2]
nll_vec <- -log(like_vec)

plot(Linf_vec, nll_vec)

