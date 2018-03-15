main_dir <- "C:\\merrill\\ADMB_course\\Day2\\2b_Age-structured"
proj_dir <- file.path(main_dir, "agestructured3")

library(R2admb)

setwd(proj_dir)
compile_admb("SCA", verbose=TRUE)

run_admb("SCA")
