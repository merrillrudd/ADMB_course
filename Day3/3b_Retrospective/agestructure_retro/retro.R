rm(list=ls())

# install.packages("R2admb")
library(R2admb)

main_dir <- "C:\\merrill\\ADMB_course"
proj_dir <- file.path(main_dir, "Day3", "3b_Retrospective", "agestructure_retro")
setwd(proj_dir)

source(file.path(main_dir, "tools.R"))

######################################
## Setup R2admb
######################################
files <- list.files(file.path(proj_dir))

## find tpl
tpl <- file.path(proj_dir, files[which(grepl(".TPL", files))])
name <- strsplit(files[grepl("TPL",files)],".TPL")[[1]][1]

# ## compile
compile_admb(name, verbose=TRUE)

# ## check executable
exe <- paste0(name,".exe")
file.exists(exe)

#########################################################################
## Run retrospectives
#########################################################################

setwd(proj_dir)

## loop over retrospective years
retro <- 0:10
rundir <- file.path(proj_dir, "run_retro")
dir.create("run_retro", showWarnings=FALSE)

for(i in 1:length(retro)){
	retrodir <- file.path(rundir, paste0("retro_", retro[i]))
	dir.create(retrodir, showWarnings=FALSE)
	setwd(retrodir)

	## copy TPL to run directory
	file.copy(from=tpl, to=retrodir, overwrite=TRUE)	

	## copy exe
	file.copy(from=file.path(proj_dir,exe), to=retrodir, overwrite=TRUE)	

	## choose data file
	data <- file.path(proj_dir, paste0(name, ".dat"))
	file.copy(from=data, to=retrodir, overwrite=TRUE)	

	#######################################################
	## run
	#######################################################
	## set seed at retro value
	run_admb(name, extra.args=paste("-retro", retro[i]))

}

	library(RColorBrewer)
	Fish <- Rec <- VB <- list()
	for(i in 1:length(retro)){
		retrodir <- file.path(rundir, paste0("retro_", retro[i]))
		setwd(retrodir)

		Fish[[i]] <- readVec("F", paste0(name, ".rep"))
		Rec[[i]] <- readVec("Recruits", paste0(name, ".rep"))
		VB[[i]] <- readVec("VulBio", paste0(name, ".rep"))
	}
	
	## colors
	col_vec <- rev(brewer.pal(length(retro)-1, "Reds"))
	
	## figure margins
	par(mfrow=c(3,1), mar=c(0,0,0,0), omi=c(1,1,0.5,0.2))

	## set up plot
	plot(x=1,y=1,type="n", xlim=c(0,length(Fish[[1]])), ylim=c(0,max(Fish[[1]])*1.5), cex.axis=1.5, xaxt="n", las=2)
	mtext(side=2, "Fishing mortality", line=3)

	## add lines for each retrospective year 
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(Fish[[i]], col=col)
	}

	## set up plot
	plot(x=1,y=1,type="n", xlim=c(0,length(Rec[[1]])), ylim=c(0,max(Rec[[1]])*1.5), cex.axis=1.5, xaxt="n", las=2)
	mtext(side=2, "Recruitment", line=3)

	## add lines for each retropsective year
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(Rec[[i]], col=col)
	}

	## set up plots
	plot(x=1,y=1,type="n", xlim=c(0,length(VB[[1]])), ylim=c(0, max(VB[[1]])*1.5), cex.axis=1.5, xaxt="n", las=2)
	mtext(side=2, "Vulnerable biomass", line=3)
	mtext(side=1, "Year", line=3)
	axis(1, cex.axis=1.5)

	## add lines for each retrospective year
	for(i in 1:length(retro)){
		col <- ifelse(i==1, "black", col_vec[i])
		lines(VB[[i]], col=col)
	}

	