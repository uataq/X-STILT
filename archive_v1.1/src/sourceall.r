#***************************************************************************************************
# source all required R functions
#***************************************************************************************************
# $Id: sourceall.r,v 1.22 2009-11-25 08:51:03 gerbig Exp $
#---------------------------------------------------------------------------------------------------
# added by Dien Wu

go <- graphics.off

# Default sourcepath is current directory
#if(!is.element("sourcepath",objects()))sourcepath <- "./src"
sourcepath <- "./src"

### --------------------- required STILT subroutines ---------------------- ###
stiltrpath <- file.path(sourcepath, "stiltR/")

# basic subroutines
source(paste(stiltrpath, "unix.r", sep=""))
source(paste(stiltrpath, "unix.shell.r", sep=""))
source(paste(stiltrpath, "getr.r", sep=""))
source(paste(stiltrpath, "assignr.r", sep=""))
source(paste(stiltrpath, "existsr.r", sep=""))
source(paste(stiltrpath, "distance.r", sep=""))

# time string, date...related --
source(paste(stiltrpath, "julian.r", sep=""))
source(paste(stiltrpath, "weekdayhr.r", sep=""))
source(paste(stiltrpath, "day.of.week.r", sep=""))
source(paste(stiltrpath, "month.day.year.r", sep=""))

# basic subroutine to generate STILT trajec
source(paste(stiltrpath, "Trajec.r", sep=""))
source(paste(stiltrpath, "getgridp.r", sep=""))

# basic subroutine to generate STILT footprint
source(paste(stiltrpath, "Trajecfoot.r", sep=""))

# updated version for generating multiple receptors in one hymodelc
source(paste(stiltrpath, "Trajecmulti.r", sep=""))

# basic subroutine to get ground hgt, wind info at receptor lat/lon
source(paste(stiltrpath, "trajwind.r", sep=""))


### -------------------- main subroutines for XSTILT ---------------------- ###
xpath <- file.path(sourcepath, "xstiltR/")

# subroutines to get receptor info, ground height and meteo files
source(paste(xpath, "ident.to.info.r", sep=""))
source(paste(xpath, "get.grdhgt.r", sep=""))
source(paste(xpath, "find.metfile.r", sep=""))
source(paste(xpath, "get.SIGUVERR.r", sep=""))
source(paste(xpath, "grab.oco2.r", sep=""))
source(paste(xpath, "area.r", sep="")) # calculate the area given lat/lon info

source(paste(xpath, "ggplot.forward.trajec.r", sep=""))
source(paste(xpath, "ggplot.map.r", sep=""))
source(paste(xpath, "default.col.r", sep=""))

# subroutines to generate backward trajec
source(paste(xpath, "get.more.namelist.r", sep=""))
source(paste(xpath, "run.backward.trajec.r", sep=""))
source(paste(xpath, "run.forward.trajec.r", sep=""))
source(paste(xpath, "find.create.dir.r", sep=""))

# subroutines that allocate receptors and allow running trajec at the same time
source(paste(xpath, "write.bash.r", sep=""))
source(paste(xpath, "allocate.recp.r", sep=""))

# topmost subroutine for XCO2 simulation
source(paste(xpath, "sim.xco2.r", sep=""))

# get satellite profiles given receptor lat/lon
source(paste(xpath, "OCO2.get.oco2info.r", sep=""))

# get satellite tracks given a city
source(paste(xpath, "OCO2.find.overpass.r", sep=""))

# X-STILT column footprint related --
source(paste(xpath, "OCO2.get.foot.r", sep=""))
source(paste(xpath, "OCO2.get.weight.funcv2.r", sep=""))
source(paste(xpath, "OCO2.weight.trajecfootv2.r", sep=""))

# X-STILT XCO2 simulation related --
source(paste(xpath, "convert.tif2nc.odiac.r", sep=""))
source(paste(xpath, "OCO2.odiac.anthro.r", sep=""))
source(paste(xpath, "OCO2.ctnrt.bio.ocean.r", sep=""))
source(paste(xpath, "OCO2.ctnrt.traj.edp.r", sep=""))
source(paste(xpath, "OCO2.apriori.r", sep=""))

# X-STILT error analysis related --
#source(paste(xpath, "ODIAC.trajfoot.r", sep=""))
#source(paste(xpath, "ODIAC.uncert.subroutine.r", sep=""))
#source(paste(xpath, "CT-NRT.bio.trajfoot.r", sep=""))
#source(paste(xpath, "CT-NRT.bg.trajfoot.r", sep=""))

# subroutines for bootstrap --
#source(paste(xpath, "OCO2.bootstrap.level.r", sep=""))
#source(paste(xpath, "OCO2.bootstrap.traj.r", sep=""))


# plotting
source(paste(xpath, "ggplot.sim.obs.r", sep=""))
