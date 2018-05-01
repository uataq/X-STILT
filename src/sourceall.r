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

# Load initial STILT subroutines
stiltrpath <- file.path(sourcepath,"stiltR/")

source(paste(stiltrpath, "assignr.r", sep=""))
source(paste(stiltrpath, "existsr.r", sep=""))
source(paste(stiltrpath, "distance.r", sep=""))
source(paste(stiltrpath, "getr.r", sep=""))
source(paste(stiltrpath, "Trajec.r", sep=""))
source(paste(stiltrpath, "Trajecfoot.r", sep=""))
source(paste(stiltrpath, "Trajecflux.r", sep=""))
source(paste(stiltrpath, "Trajecmod.r", sep=""))
source(paste(stiltrpath, "Trajecmulti.r", sep=""))  # updated version for generating multiple receptors in one hymodelc
source(paste(stiltrpath, "Trajecvprm.r", sep=""))
source(paste(stiltrpath, "trajwind.r", sep=""))
source(paste(stiltrpath, "weekdayhr.r", sep=""))
source(paste(stiltrpath, "unix.r", sep=""))
source(paste(stiltrpath, "unix.shell.r", sep=""))

# main subroutines for XSTILT
xpath <- file.path(sourcepath,"xstiltR/")

source(paste(xpath, "ident.to.info.r", sep=""))
source(paste(xpath, "get.grdhgt.r", sep=""))
source(paste(xpath, "find.metfile.r", sep=""))
source(paste(xpath, "run.backward.trajec.r", sep=""))  # run backward trajec


source(paste(xpath, "OCO2.get.oco2info.r", sep=""))  # get satellite profiles given receptor lat/lon
source(paste(xpath, "OCO2.find.overpass.r", sep=""))
source(paste(xpath, "OCO2.get.xco2.r", sep=""))
source(paste(xpath, "OCO2.get.weight.funcv2.r", sep=""))
source(paste(xpath, "OCO2.weight.trajecfootv2.r", sep=""))
source(paste(xpath, "OCO2.CT-NRT.oceanv2.r", sep=""))
source(paste(xpath, "OCO2.CT-NRT.biov2.r", sep=""))
source(paste(xpath, "OCO2.get.foot.r", sep=""))
source(paste(xpath, "OCO2.CT-NRT.background.r", sep=""))
source(paste(xpath, "OCO2.apriori.r", sep=""))
source(paste(xpath, "ODIAC.trajfoot.r", sep=""))
source(paste(xpath, "ODIAC.subroutine.r", sep=""))
source(paste(xpath, "ODIAC.uncert.subroutine.r", sep=""))
source(paste(xpath, "CT-NRT.bio.trajfoot.r", sep=""))
source(paste(xpath, "CT-NRT.bg.trajfoot.r", sep=""))

# subroutines for bootstrap
source(paste(xpath, "OCO2.bootstrap.level.r", sep=""))
source(paste(xpath, "OCO2.bootstrap.traj.r", sep=""))
