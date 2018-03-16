wd <- "C:\\merrill\\ADMB_course\\HER"
setwd(wd)

library(R2admb)

compile_admb("her", verbose=TRUE)
run_admb("her")