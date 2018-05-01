# create a namelist for running X-STILT trajectories, by DW, 04/18/2018
# required functions: ident.to.info()

##### source all functions and load all libraries
# change working directory ***
workpath <- "/uufs/chpc.utah.edu/common/home/lin-group1/wde/github/XSTILT"
setwd(workpath) # move to working directory

source(file.path(workpath, "src/sourceall.r"))  # source all functions
library(ncdf4)
library(ggplot2)

#------------------------------ STEP 1 --------------------------------------- #
#### CHOOSE CITIES, TRACKS AND OCO-2 LITE FILE VERSION ***
index <- 1
site <- c("Riyadh", "Cairo", "PRD")[index]
met <- c("1km", "gdas1", "gdas0p5")[3]     # customized WRF, 1 or 0.5deg GDAS
tt <- 2  # which track to model
oco2.version <- c("b7rb","b8r")[1]  # OCO-2 version ***

#### CHOOSE OVERPASSED TIMESTR ***
# input all track numbers to be modeled, can be YYYYMMDD OR YYYYMMDDHH ***
riyadh.timestr <- c("2014122710", "2014122910", "2015012810", "2015081710",
                    "2015111210", "2015121610", "2016011510", "2016021610",
                    "2016072510")
cairo.timestr <- c("20150228", "20150318", "20160224")
prd.timestr <- "2015011505"

# final timestr, YYYYMMDD and file strings for trajec
track.timestr <- c(riyadh.timestr[tt], cairo.timestr[tt], prd.timestr)[index]
filestr <- paste(substr(track.timestr, 1, 4), "x", substr(track.timestr, 5, 6),
                 "x", substr(track.timestr, 7, 8), sep="")

# spatial domains placing receptors, help select OCO-2 data ***
if(site == "Riyadh"){minlat <- 23; maxlat <- 26; minlon <- 45; maxlon <- 50}
if(site == "Cairo"){minlat <- 29; maxlat <- 31; minlon <- 30; maxlon <- 32}
if(site == "PRD"){minlat <- 22; maxlat <- 23; minlon <- 112; maxlon <- 115}

#------------------------------ STEP 2 --------------------------------------- #
#### TURN ON/OFF FLAGS for XSTILT setups ***
selTF <-F       # true for only use 1-day trajs.
nummodel <- 10	# copy for running trajwind() or trajec()
columnTF <- T   # whether a column receptor or fixed receptor
forwardTF <- F  # forward or backward traj

#### INPUT/CHANGE PATHS ***
# if no "output" directory, create output directory
if(length(list.files(path=workpath,pattern="output"))==0)system("mkdir output")
outpath <- file.path(workpath, "output")  # store STILT trajec

# inpath for OCO2, ODIAC, CT..etc.
inpath <- "/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_input"

# OCO-2 path
ocopath <- file.path(inpath , paste("OCO-2/OCO2_lite_",oco2.version,sep=""))

# Where to run STILT, where Copy lies ***
rundir <- "/uufs/chpc.utah.edu/common/home/lin-group1/wde/STILT_modeling/STILT_Exe/"

#------------------------------ STEP 3 --------------------------------------- #
#### METFILES and X-STILT configurations ***
# path for the ARL format of WRF and GDAS, CANNOT BE TOO LONG ***
metpath <- "/uufs/chpc.utah.edu/common/home/u0947337/"

# get metfile for generating backward trajectories
metfile <- find.metfile(timestr=track.timestr, nhrs=nhrs, metpath=metpath,
                        met=met, trajecTF=T) # do not change trajecTF=T !!!

### numbers of hours for generating back trajec
nhrs <- -72
if(selTF)nhrs <- -24  # 1day if  selTF
overwrite <- F	      # T:rerun hymodelc, even if particle location object found
                      # F:re-use previously calculated particle location object
delt <- 2		          # fixed timestep [min]; set =0 for dynamic timestep

### Set model receptors, AGLs
# if release particles from fixed levels
agl <- 10
npar <- 3700

## SET COLUMN RECEPTORS, using vector form
# INPUT min, max heights and dh for STILT levels
# will re-calculate "agl" in subroutine based on parameters shown below
if(columnTF){    # release particles from a column
  minagl <- 0        # min, max agl heights, and cutoff level zi, in METERS
  maxagl <- 6000
  zi <- 3000
  dh <- c(100, 500)  # vertical spacing below and above cutoff level, in METERS
  dpar <- 100        # particle numbers per level
}

# eyeball lat range for enhanced XCO2, or grab from available forward-time runs
midlat <- c(24.5, 25)  # more receptors during this lat range

#------------------------------ STEP 4 --------------------------------------- #
#### create a namelist (data.frame) including all variables
# and pass to subroutines
namelist<-data.frame(met=met,site=site,timestr=track.timestr,filestr=filestr,
                     minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon,
                     ocopath=ocopath, outpath=outpath, rundir=rundir,
                     metpath=metpath, metfile=metfile, selTF=selTF, nhrs=nhrs,
                     columnTF=columnTF, forwardTF=forwardTF, nummodel=nummodel,
                     overwrite=overwrite,delt=delt,minagl=minagl,dh.lower=dh[1],
                     maxagl=maxagl,zi=zi,dpar=dpar,agl=agl,dh.upper=dh[2],
                     midlat.left=midlat[1], midlat.right=midlat[2],
                     stringsAsFactors=F)

# store traj.list as RData file
filename <- paste("trajec_list_", track.timestr, "_", site, ".txt", sep="")
filename <- file.path(workpath, filename)
write.table(x=t(namelist), file=filename, sep="\n", col.names=F, quote=T)

#------------------------------ STEP 5 --------------------------------------- #
### call run.backward.trajec to start running trajectories
source(file.path(workpath, "src/sourceall.r"))  # source all functions
if(forwardTF==F)back.info <- run.backward.trajec(namelist = namelist, plotTF=F)
if(forwardTF==T)back.info <- run.forward.trajec(namelist = namelist, plotTF=F)




# end of running trajec
