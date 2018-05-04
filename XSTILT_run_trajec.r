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

# spatial domains placing receptors and city center, help select OCO-2 data ***
# in form of "lat.lon <- c(minlat, maxlat, minlon, maxlon, city.lat, city.lon)"
# oco2.hr for overpass hour in UTC
if(site == "Riyadh"){lat.lon <- c(23, 26, 45, 50, 24.71, 46.75); oco2.hr <- 10}
if(site == "Cairo"){lat.lon <- c(29, 31, 30, 32, NA, NA); oco2.hr <- NA}
if(site == "PRD"){lat.lon <- c(22, 23, 112, 115, NA, NA); oco2.hr <- NA}

#------------------------------ STEP 2 --------------------------------------- #
#### TURN ON/OFF FLAGS for XSTILT setups ***
selTF <- F       # true for only use 1-day trajs.
columnTF <- F    # whether a column receptor or fixed receptor
forwardTF <- T   # forward or backward traj, if forward, release from a box
windowTF <- T    # whether to release particles every ?? minutes, "dhr" in hours
uncertTF <- T    # whether add wind error component to generate trajec
overwrite <- T	      # T:rerun hymodelc, even if particle location object found
                      # F:re-use previously calculated particle location object
delt <- 2             # fixed timestep [min]; set =0 for dynamic timestep

### change parameters according to above flags
nhrs <- -72           # number of hours backward (-) or forward (+)
if(selTF)nhrs <- -24  # 1day if  selTF

# copy for running trajwind() or trajec()
nummodel <- 1; dxyp <- NULL
if(forwardTF){
  dxyp <- 0.2*2       # +/-dxyp, +/-dxyp around the city center
  nummodel <- 997     # use TN's hymodelc to run forward box trajec
  nhrs <- 12          # number of hours backward (-) or forward (+)
}

## if allow for time window for releasing particles
dhr <- NULL; dtime <- NULL
if(windowTF){
  dhr <- 0.5     # release particle every 30 mins
  # FROM 10 hours ahead of sounding hour (-10), TO sounding hour (0)
  dtime <- seq(-10, 0, dhr)
  cat("Release particles", dtime[1], "hrs ahead & every", dhr * 60, "mins...\n")
}

#### obtaining wind errors, transport error component
# SD in wind errors, correlation timescale, horizontal and vertical lengthscales
if(uncertTF){
  cat("Using wind error component...\n")

  # intput correlation lengthscale (in m) and timescales (in mins) ***
  if(met == "gdas0p5"){TLuverr <- 1*60; zcoruverr <- 600; horcoruverr <- 40}
  if(met == "gdas1"){TLuverr <- 2.4*60; zcoruverr <- 700; horcoruverr <- 97}

  ## add errors, mainly siguverr, create a subroutine to compute siguverr
  # from model-data wind comparisons
  errpath <- file.path(homedir,"lin-group1/wde/STILT_input/wind_err")
  errpath <- file.path(errpath, tolower(site), tolower(met))

  # call get.SIGUVERR() to interpolate most near-field wind errors
  err.info <- get.SIGUVERR(site=site, timestr=track.timestr, gdaspath=errpath,
                           nfTF=FALSE,forwardTF=F)

  # if no wind error found, return NA for err.info
  # make a conservative assumption about the wind error, for the Middle East
  if(length(err.info)!=2){siguverr <- 1.8}

  met.rad <- err.info[[1]]
  siguverr <- as.numeric(err.info[[2]][1])
  u.bias <- as.numeric(err.info[[2]][2])
  v.bias <- as.numeric(err.info[[2]][3])

  cat(paste("SIGUVERR:", signif(siguverr,3), "m/s..\n"))
  cat(paste("u.bias:", signif(u.bias,3), "\n"))
  cat(paste("v.bias:", signif(u.bias,3), "\n"))

}else{

  cat("NOT using wind error component...\n")
  siguverr <- NULL; TLuverr <- NULL; horcoruverr <- NULL; zcoruverr <- NULL
}

#------------------------------ STEP 3 --------------------------------------- #
#### create output directory
# no need to change "outpath" below, it automatically creates one
# if there's no "output" directory, call find.create.dir() to create one-->
find.create.dir(path=workdir, dir="output", workdir=workdir)
outpath <- file.path(workdir, "output")  # output directory

# if no "Traj" directory for storing trajec .Rdata,
# call find.create.dir() to create one
find.create.dir(path=outpath, dir="Traj", workdir=workdir)
trajpath <- file.path(outpath, "Traj") # store trajec .RData file in this dir

if(forwardTF){
  # if running FORWARD trajec from city center, create "forw_traj"
  find.create.dir(path=trajpath, dir="forw_traj", workdir=workdir)
  trajpath <- file.path(trajpath, "forw_traj") # store trajec .RData file here
}else{
  # if running BACKWARD trajec from sounding lat/lon, create "orig_traj"
  find.create.dir(path=trajpath, dir="orig_traj", workdir=workdir)
  trajpath <- file.path(trajpath, "orig_traj") # store trajec .RData file here
}  # end if forwardTF

# path for input data, OCO-2 Lite file
ocostr   <- paste("OCO-2/OCO2_lite_", oco2.version, sep="")
ocopath  <- file.path(homedir, "lin-group1/wde/STILT_input", ocostr)

#### CHANGE METFILES ***
# path for the ARL format of WRF and GDAS, CANNOT BE TOO LONG ***
metpath <- file.path(homedir, "u0947337/")

# get metfile for generating backward trajectories
metfile <- find.metfile(timestr=track.timestr, nhrs=nhrs, metpath=metpath,
                        met=met, trajecTF=T) # do not change trajecTF=T !!!

# Where to run STILT, where Copy lies
rundir <- file.path(homedir, "lin-group1/wde/STILT_modeling/STILT_Exe/")

#------------------------------ STEP 4 --------------------------------------- #
#### Set model receptors, AGLs and particel numbers ***
### 1) if release particles from fixed levels
# for backward fixed-level runs OR forward fixed-level runs
# agl can be a vector, meaning releasing particles from several fixed level
# but if trying to represent air column, use columnTF=T, see below

agl <- 10
npar<- 1000            # par for each time window for forward box runs
minagl <- NA; maxagl <- NA; zi <- NA; dh <- NA; dpar <- NA

#### 2) SET COLUMN RECEPTORS, using vector form ***
 # if release particles from a column
if(columnTF){

  # INPUT min, max heights and dh for STILT levels, in METERS
  minagl <- 0; maxagl <- 6000; zi <- 3000

  # vertical spacing below and above cutoff level, in METERS
  dh <- c(100, 500)
  dpar <- 100           # particle numbers per level
  agl <- c(seq(minagl, zi, dh[1]), seq(zi+dh[2], maxagl, dh[2]))
  nlev <- length(agl)
  npar <- nlev * dpar   # total number of particles
}

## eyeball lat range for enhanced XCO2, or grab from available forward-time runs
# place denser receptors during this lat range (40pts in 1deg)
peak.lat <- c(24.5, 25)  # in deg North (+), required for backward run
bw.bg <- 1/20    # binwidth over small enhancements, e.g., every 20 pts in 1deg
bw.peak <- 1/40  # binwidth over large enhancements, during "peak.lat"

# recp.index: how to pick receptors from soundings
recp.index <- c(seq(lat.lon[1], peak.lat[1], bw.bg),
                seq(peak.lat[1], peak.lat[2], bw.peak),
                seq(peak.lat[1], lat.lon[2], bw.bg))

#------------------------------ STEP 5 --------------------------------------- #
#### !!! NO NEED TO CHANGE ANYTHING LISTED BELOW -->
# create a namelist including all variables
# namelist required for generating trajec
namelist <- list(met=met,site=site,timestr=track.timestr,filestr=filestr,
                 lat.lon=lat.lon, ocopath=ocopath, outpath=trajpath,
                 rundir=rundir,delt=delt, metpath=metpath, metfile=metfile,
                 selTF=selTF, nhrs=nhrs, columnTF=columnTF, forwardTF=forwardTF,
                 uncertTF=uncertTF, windowTF=windowTF, nummodel=nummodel,
                 overwrite=overwrite, minagl=minagl, maxagl=maxagl, zi=zi,
                 dpar=dpar, npar=npar, agl=agl, dh=dh, dxyp=dxyp, dtime=dtime,
                 oco2.hr=oco2.hr, siguverr=siguverr, TLuverr=TLuverr,
                 zcoruverr=zcoruverr, horcoruverr=horcoruverr,
                 recp.index=recp.index, stringsAsFactors=F)

# store namelist to output directory
filenm <- paste("trajlist_", site, "_", track.timestr, sep="")
if(forwardTF==F)filenm <- paste(filenm, "_backward.txt", sep="")
if(forwardTF==T)filenm <- paste(filenm, "_forward.txt", sep="")
filenm <- file.path(outpath, filenm)  # link to path
write.table(x=t(namelist), file=filenm, sep="\n", col.names=F, quote=T)

#------------------------------ STEP 7 --------------------------------------- #
#### call run.backward.trajec to start running trajectories
# trajectories will be stored in ./output
# plotTF for whether plotting OCO-2 observed data if calling backward subroutine
if(forwardTF==F){
  cat("Generating backward trajec...\n")
  back.info <- run.backward.trajec(namelist = namelist, plotTF=F)
}

# plotTF for whether plotting urban plume & obs XCO2 if calling forward function
# !!! if forward, release particles from a box around the city center
if(forwardTF==T){
  cat("Generating forward trajec...\n")
  source(file.path(workdir, "src/sourceall.r"))  # source all functions
  forw.info <- run.forward.trajec(namelist = namelist, plotTF=T)
}

# end of running trajec
