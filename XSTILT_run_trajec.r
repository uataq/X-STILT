# create a namelist for running X-STILT trajectories, by DW, 04/18/2018
# required subroutines:
# find.metfile(); run.backward.trajec(); run.forward.trajec()

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- "/uufs/chpc.utah.edu/common/home"
workdir <- file.path(homedir, "lin-group1/wde/github/XSTILT")
setwd(workdir) # move to working directory

source(file.path(workdir, "src/sourceall.r"))  # source all functions
library(ncdf4); library(ggplot2)

#------------------------------ STEP 1 --------------------------------------- #
#### CHOOSE CITIES, TRACKS AND OCO-2 LITE FILE VERSION ***
index <- 1
site  <- c("Riyadh", "Cairo", "PRD")[index]
met   <- c("1km", "gdas1", "gdas0p5")[3]  # customized WRF, 1 or 0.5deg GDAS
tt    <- 2                                # which track to model
oco2.version <- c("b7rb","b8r")[1]        # OCO-2 version

#### CHOOSE OVERPASSED TIMESTR ***
# input all track numbers to be modeled, can be YYYYMMDD OR YYYYMMDDHH ***
# track timestr can also be grabbed from another script
riyadh.timestr <- c("20141227", "20141229", "20150128", "20150817", "20151112",
                    "20151216", "20160115", "20160216", "20160725", "20161031")
cairo.timestr <- c("20150228", "20150318", "20160224"); prd.timestr<- "20150115"

# final timestr, YYYYMMDD and file strings for trajec
track.timestr <- c(riyadh.timestr[tt], cairo.timestr[tt], prd.timestr)[index]
filestr <- paste(substr(track.timestr,1,4), "x", substr(track.timestr,5,6), "x",
                 substr(track.timestr,7,8), sep="")

# spatial domains placing receptors, help select OCO-2 data ***
if(site == "Riyadh"){minlat <- 23; maxlat <- 26; minlon <- 45; maxlon <- 50}
if(site == "Cairo"){minlat <- 29; maxlat <- 31; minlon <- 30; maxlon <- 32}
if(site == "PRD"){minlat <- 22; maxlat <- 23; minlon <- 112; maxlon <- 115}

#------------------------------ STEP 2 --------------------------------------- #
#### TURN ON/OFF FLAGS for XSTILT setups ***
selTF <-F       # true for only use 1-day trajs.
nummodel <- 10	# copy for running trajwind() or trajec()
columnTF <- T   # whether a column receptor or fixed receptor
forwardTF<- F   # forward or backward traj

#### create output directory
# no need to change "outpath" below, it automatically creates one
# if there's no "output" directory, call find.create.dir() to create one-->
find.create.dir(path=workdir, dir="output", workdir=workdir)
outpath <- file.path(workdir, "output")  # output directory

# if no "Traj" directory for storing trajec .Rdata,
# call find.create.dir() to create one
find.create.dir(path=outpath, dir="Traj", workdir=workdir)
trajpath <- file.path(outpath, "Traj") # store trajec .RData file in this dir
find.create.dir(path=trajpath, dir="orig_traj", workdir=workdir)
trajpath <- file.path(trajpath, "orig_traj") # store trajec .RData file in this dir

# path for input data, OCO-2 Lite file
ocostr   <- paste("OCO-2/OCO2_lite_", oco2.version, sep="")
ocopath  <- file.path(homedir, "lin-group1/wde/STILT_input", ocostr)

#------------------------------ STEP 3 --------------------------------------- #
#### CHANGE METFILES ***
# path for the ARL format of WRF and GDAS, CANNOT BE TOO LONG ***
metpath <- file.path(homedir, "u0947337/")

# get metfile for generating backward trajectories
metfile <- find.metfile(timestr=track.timestr, nhrs=nhrs, metpath=metpath,
                        met=met, trajecTF=T) # do not change trajecTF=T !!!

# Where to run STILT, where Copy lies
rundir <- file.path(homedir, "lin-group1/wde/STILT_modeling/STILT_Exe/")

#### CHANGE X-STILT configurations ***
nhrs <- -72           # number of hours backward (-) or forward (+)
if(selTF)nhrs <- -24  # 1day if  selTF
overwrite <- F	      # T:rerun hymodelc, even if particle location object found
                      # F:re-use previously calculated particle location object
delt <- 2             # fixed timestep [min]; set =0 for dynamic timestep

#### Set model receptors, AGLs and particel numbers ***
# if release particles from fixed levels

# agl can be a vector, meaning releasing particles from several fixed level
# but if trying to represent air column, use columnTF=T, see below
agl <- 10
npar<- 3700

#### SET COLUMN RECEPTORS, using vector form
# will re-calculate "agl" and "npar" in run.backward.trajec()
# based on parameters shown below
if(columnTF){        # if release particles from a column

  # INPUT min, max heights and dh for STILT levels, in METERS
  minagl <- 0
  maxagl <- 6000
  zi <- 3000

  # vertical spacing below and above cutoff level, in METERS
  dh <- c(100, 500)

  # particle numbers per level
  dpar <- 100
}

# eyeball lat range for enhanced XCO2, or grab from available forward-time runs
# place denser receptors during this lat range (40pts in 1deg)
midlat <- c(24.5, 25)  # in deg North (+)

#------------------------------ STEP 4 --------------------------------------- #
#### NO NEED TO CHANGE ANYTHING LISTED BELOW -->
# create a namelist (data.frame) including all variables
# namelist required for generating trajec
namelist<-data.frame(met=met,site=site,timestr=track.timestr,filestr=filestr,
                     minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon,
                     ocopath=ocopath, outpath=trajpath, rundir=rundir,delt=delt,
                     metpath=metpath, metfile=metfile, selTF=selTF, nhrs=nhrs,
                     columnTF=columnTF, forwardTF=forwardTF, nummodel=nummodel,
                     overwrite=overwrite, minagl=minagl, maxagl=maxagl, zi=zi,
                     dpar=dpar, agl=agl,dh.lower=dh[1], dh.upper=dh[2],
                     midlat.left=midlat[1], midlat.right=midlat[2],
                     stringsAsFactors=F)

# store namelist to output directory
filename <- paste("trajlist_", site, "_", track.timestr, ".txt", sep="")
filename <- file.path(outpath, filename)
write.table(x=t(namelist), file=filename, sep="\n", col.names=F, quote=T)

#------------------------------ STEP 5 --------------------------------------- #
#### call run.backward.trajec to start running trajectories
# trajectories will be stored in ./output
# plotTF for whether plotting OCO-2 observed data if calling backward subroutine
if(forwardTF==F)back.info <- run.backward.trajec(namelist = namelist, plotTF=F)

# plotTF for whether plotting urban plume & obs XCO2 if calling forward function
if(forwardTF==T)forw.info <- run.forward.trajec(namelist = namelist, plotTF=F)


# end of running trajec
