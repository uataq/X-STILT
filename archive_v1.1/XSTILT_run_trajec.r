# create a namelist for running X-STILT trajectories
# written by DW, 04/18/2018
# latest modification on 05/10/2018
# clean code on 06/12/2018

# required subroutines:
# find.metfile(); get.SIGUVERR(); find.create.dir(); get.more.namelist();
# allocate.recp(); write.bash(); run.backward.trajec(); run.forward.trajec()

#### source all functions and load all libraries
# CHANGE working directory
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT')
setwd(workdir) # move to working directory

# source all functions
source(file.path(workdir, 'src/sourceall.r'))
library(ncdf4); library(ggplot2)

#------------------------------ STEP 1 --------------------------------------- #
#### CHOOSE CITIES, TRACKS AND OCO-2 LITE FILE VERSION ***
site  <- 'Riyadh'
met   <- c('1km', 'gdas1', 'gdas0p5')[3]  # customized WRF, 1 or 0.5deg GDAS
tt    <- 2                                # which track to model
oco2.version <- c('b7rb', 'b8r')[1]       # OCO-2 version

#### CHOOSE OVERPASSED TIMESTR ***
# input all track numbers to be modeled, can be YYYYMMDD OR YYYYMMDDHH ***
# track timestr can also be grabbed from another script
riyadh.timestr <- c('20141227', '20141229', '20150128', '20150817', '20151112',
  '20151216', '20160115', '20160216', '20160725', '20161031')

# final timestr, YYYYMMDD and file strings for trajec
track.timestr <- riyadh.timestr[tt]
filestr <- paste0(substr(track.timestr, 1, 4), 'x', substr(track.timestr, 5, 6),
  'x', substr(track.timestr, 7, 8))

# spatial domains placing receptors and city center, help select OCO-2 data ***
# in form of 'lat.lon <- c(minlat, maxlat, minlon, maxlon, city.lat, city.lon)'
# oco2.hr for overpass hour in UTC
if (site == 'Riyadh'){lat.lon <- c(23, 26, 45, 50, 24.71, 46.75); oco2.hr <- 10}
track.timestr <- paste0(track.timestr, oco2.hr)


#------------------------------ STEP 2 --------------------------------------- #
#### TURN ON/OFF FLAGS for XSTILT setups ***
columnTF  <- T    # whether a column receptor or fixed receptor
forwardTF <- F    # forward or backward traj, if forward, release from a box
windowTF  <- F    # whether to release particles every ?? minutes, 'dhr' in hours
uncertTF  <- F    # whether add wind error component to generate trajec
overwrite <- T	  # T:rerun hymodelc, even if particle location object found
                  # F:re-use previously calculated particle location object
mpcTF <- F        # true for running trajec on multiple STILT copies
delt  <- 0        # fixed timestep [min]; set =0 for dynamic timestep

### change parameters according to above flags
nhrs <- -12      # number of hours backward (-) or forward (+)
nummodel <- 1    # copy for running trajwind() or trajec()
if (mpcTF) nummodel <- seq(11, 17)      # can be a vector denoting copy numbers

## if release particles from a box, +/-dxyp, +/-dxyp around the city center
dxyp <- NULL
if (forwardTF) {
  dxyp <- 0.2*2
  nummodel <- 997     # use TN's hymodelc to run forward box trajec
  nhrs <- 12          # number of hours backward (-) or forward (+)
}

## if allow for time window for releasing particles
dhr <- NULL; dtime <- NULL
if (windowTF) {
  dhr <- 0.5     # release particle every 30 mins
  # FROM 10 hours ahead of sounding hour (-10), TO sounding hour (0)
  dtime <- seq(-10, 0, dhr)
  cat('Release particles', dtime[1], 'hrs ahead & every', dhr * 60, 'mins...\n')
}

#### obtaining wind errors, transport error component
# SD in wind errors, correlation timescale, horizontal and vertical lengthscales
if (uncertTF) {
  cat('Using wind error component...\n')

  # intput correlation lengthscale (in m) and timescales (in mins) ***
  if (met == 'gdas0p5') {TLuverr <- 1*60; zcoruverr <- 600; horcoruverr <- 40}
  if (met == 'gdas1') {TLuverr <- 2.4*60; zcoruverr <- 700; horcoruverr <- 97}

  ## add errors, mainly siguverr, create a subroutine to compute siguverr
  # from model-data wind comparisons
  errpath <- file.path(homedir, 'lin-group5/wde/STILT_input/wind_err')
  errpath <- file.path(errpath, tolower(site), tolower(met))

  # call get.SIGUVERR() to interpolate most near-field wind errors
  err.info <- get.SIGUVERR(site = site, timestr = track.timestr,
                           gdaspath = errpath, nfTF = FALSE, forwardTF = F)
  if (length(err.info) != 2) {
    # if no wind error found, return NA for err.info
    # make a conservative assumption about the wind error, for the Middle East
    siguverr <- 1.8

  } else {
    met.rad <- err.info[[1]]
    siguverr <- as.numeric(err.info[[2]][1])
    u.bias <- as.numeric(err.info[[2]][2])
    v.bias <- as.numeric(err.info[[2]][3])

    cat(paste('SIGUVERR:', signif(siguverr,3), 'm/s..\n'))
    cat(paste('u.bias:', signif(u.bias,3), '\n'))
    cat(paste('v.bias:', signif(u.bias,3), '\n'))
  } # end if err.info

}else{

  cat('NOT using wind error component...\n')
  siguverr <- NULL; TLuverr <- NULL; horcoruverr <- NULL; zcoruverr <- NULL
}  # end if uncertTF


#------------------------------ STEP 3 --------------------------------------- #
#### create output directory
# no need to change 'outpath' below, it automatically creates one
# if there's no 'output' directory, call find.create.dir() to create one-->
find.create.dir(path = workdir, dir = 'output', workdir = workdir)
outpath <- file.path(workdir, 'output')   # output directory

# if no 'Traj' directory for storing trajec .Rdata,
# call find.create.dir() to create one
find.create.dir(path = outpath, dir = 'Traj', workdir = workdir)
trajpath <- file.path(outpath, 'Traj') # store trajec .RData file in this dir

if (forwardTF) {
  # if running FORWARD trajec from city center, create 'forw_traj'
  find.create.dir(path = trajpath, dir='forw_traj', workdir = workdir)
  trajpath <- file.path(trajpath, 'forw_traj') # store trajec .RData file here

}else{
  # if running BACKWARD trajec from sounding lat/lon, create 'orig_traj'
  find.create.dir(path = trajpath, dir='orig_traj', workdir = workdir)
  trajpath <- file.path(trajpath, 'orig_traj') # store trajec .RData file here
}  # end if forwardTF

## change paths--
# path for input data, OCO-2 Lite file
oco2.str   <- paste0('OCO-2/L2/OCO2_lite_', oco2.version)
oco2.path  <- file.path(homedir, 'lin-group5/wde/input_data', oco2.str)

# path for the ARL format of WRF and GDAS, CANNOT BE TOO LONG
metpath <- file.path(homedir, 'u0947337', met)

# get metfile for generating backward trajectories
# met.format: met file convention, e.g., '%Y%m%d_gdas0p5'
met.format <- paste0('%Y%m%d_', met)
metfile <- find.metfile(timestr = track.timestr, nhrs = nhrs, metpath = metpath,
  met.format = met.format)

# Where to run STILT, where Copy lies
rundir <- file.path(homedir, 'lin-group5/wde/STILT_modeling/STILT_Exe/')


#------------------------------ STEP 4 --------------------------------------- #
#### Set model receptors, AGLs and particel numbers ***
if (columnTF == T) {
  # 1) SET COLUMN RECEPTORS, if release particles from a column
  # INPUT min, max heights and dh for STILT levels, in METERS
  minagl <- 0
  maxagl <- 6000
  zi     <- 3000            # cutoff level, determined by mean PBL
  dh     <- c(100, 500)     # vertical spacing below and above cutoff level [m]
  dpar   <- 100             # particle numbers per level
  agl    <- c(seq(minagl, zi, dh[1]), seq(zi+dh[2], maxagl, dh[2]))
  nlev   <- length(agl)     # number of levels
  npar   <- nlev * dpar     # total number of particles

} else {
  # 2) if release particles from fixed levels
  # for backward fixed-level runs OR FORWARD fixed-level runs
  # if agl is a vector, meaning releasing particles from several fixed level
  # but if trying to represent air column, use columnTF = T (see below)
  agl    <- 10
  npar   <- 1000
  minagl <- NA; maxagl <- NA; zi <- NA; dh <- NA; dpar <- NA
}

## eyeball lat range for enhanced XCO2, or grab from available forward-time runs
# whether select receptors; or simulate all soundings
filterTF <- T

# place denser receptors during this lat range (40pts in 1deg)
peak.lat <- c(24.5, 25)  # in deg North (+), required for backward run
bw.bg   <- 1/20   # binwidth over small enhancements, e.g., every 20 pts in 1deg
bw.peak <- 1/40   # binwidth over large enhancements, during 'peak.lat'

# recp.index: how to pick receptors from soundings
recp.index <- c(seq(lat.lon[1], peak.lat[1], bw.bg),
  seq(peak.lat[1], peak.lat[2], bw.peak), seq(peak.lat[1], lat.lon[2], bw.bg))

# or overwrite recp.index with the latitude you wanted to model
recp.index <- 24.5444

#------------------------------ STEP 5 --------------------------------------- #
#### !!! NO NEED TO CHANGE ANYTHING LISTED BELOW -->
# create a namelist including all variables
# namelist required for generating trajec
namelist <- list(
  agl = agl, columnTF = columnTF, delt = delt, dh = dh, dpar = dpar,
  dtime = dtime, dxyp = dxyp, filestr = filestr, filterTF = filterTF,
  forwardTF = forwardTF, homedir = homedir, horcoruverr = horcoruverr,
  lat.lon = lat.lon, maxagl = maxagl, minagl = minagl, met = met,
  metfile = metfile, metpath = metpath, nhrs = nhrs, npar = npar,
  nummodel = nummodel, oco2.hr = oco2.hr, oco2.path = oco2.path,
  oco2.version = oco2.version, outpath = trajpath, overwrite = overwrite,
  recp.index = recp.index, rundir = rundir, siguverr = siguverr, site = site,
  timestr = track.timestr, TLuverr = TLuverr, uncertTF = uncertTF,
  windowTF = windowTF, workdir = workdir, zcoruverr = zcoruverr, zi = zi)

#------------------------------ STEP 6 --------------------------------------- #
#### call get.more.namelist() to get more info about receptors
# trajectories will be stored in ./output

if (forwardTF == F) {
  # further read OCO-2 data
  # then get receptor info and other info for running trajec
  # plotTF for whether plotting OCO-2 observed data
  namelist <- get.more.namelist(namelist = namelist, plotTF = F)

  ### prepare for creating run files, each X-STILT CONTROL r scripts,
  # trajlist & bash files that are needed to generate trajec on computing nodes
  cat('Checking and distributing copies...\n')

  # automatically creates RUN, CONTL, TRAJLIST and update namelist
  namelist$rmTF <- F  # whether remove previous files
  namelist <- allocate.recp(namelist)

  ## set up job info and write bash script
  job.time   <- '12:00:00'
  num.nodes  <- 1
  num.cores  <- length(namelist$nummodel)
  account    <- 'lin-kp'
  partition  <- account
  email      <- 'dien.wu@utah.edu'
  email.type <- 'FAIL,BEGIN,END'

  # write bash file --> config --> run_XSTILT and update namelist:
  bash.file <- write.bash(namelist, job.time, num.nodes, num.cores, account,
    partition, email, email.type)
  cat('Start generating backward trajec...\n')
  system(paste0('sbatch ', bash.file))
  q()
} # if backward


#------------------------------ STEP 7 --------------------------------------- #
# if for forward time runs
# plotTF for whether plotting urban plume & obs XCO2 if calling forward function
# !!! if forward, release particles from a box around the city center

if (forwardTF == T) {
  # store namelist to output directory
  filenm <- paste0('trajlist_', site, '_', track.timestr, '_forward.txt')
  filenm <- file.path(outpath, filenm)  # link to path
  write.table(x = t(namelist), file = filenm, sep = '\n', col.names = F, quote = T)

  cat('Start generating forward trajec...\n')
  forw.info <- run.forward.trajec(namelist = namelist, plotTF = T)
}

## end of running trajec
