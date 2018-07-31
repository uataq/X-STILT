# create a namelist for running X-STILT trajectories
# written by DW, 04/18/2018

# --------------------------------- updates ---------------------------------
# now build upon Ben's STILT-R version 2 codes, DW, 05/23/2018
# !!! need to clear up codes for forward run, based on Ben's parallel computing
# Add plotting scripts for footprint and XCO2 enhancements, DW, 06/28/2018
#
# Add STILTv1 dependencies and work well with existing framework and add dmassTF
# for mass violation corrections, DW, 07/17/2018
#
# Complement horizontal transport error module with flag 'run_hor_err'
#   with additional subroutines in /src/oco2-xstilt/trans_err, DW, 07/23/2018
#
# Complement vertical trans error module, with flag 'run_ver_err',
# ziscale when scaling PBL heights and break trans error part (step 3) into
# horizontal wind error and vertical trans error (via PBL perturb), DW, 07/25/2018
#
# Complement forward runs from a box around the city, DW, 07/27/2018
# this requires the AER_NOAA_branch modified hymodelc
#
# remove forwardTF and background estimates out of this code,
# for forward run, see 'compute_bg_xco2.r', DW, 07/30/2018
#
# for generating trajec with horizontal error compoennt,
# use the same lat/lon from original trajec, DW, 07/31/2018
# -----------------------------------------------------------------------------

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT/latest_v2.0')
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions

#------------------------------ STEP 1 --------------------------------------- #
# insert target city
site <- 'Riyadh'

# OCO-2 version, path
oco2.ver <- c('b7rb', 'b8r')[2]  # OCO-2 version
input.path <- file.path(homedir, 'lin-group5/wde/input_data')
output.path <-file.path(homedir, 'lin-group5/wde/github/result')

oco2.path <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))
sif.path <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_SIF_', oco2.ver))
txtpath <- file.path(output.path, 'oco2_overpass')

# date range for searching OCO-2 tracks, min, max YYYYMMDD
date.range <- c('20140101', '20181231')

# vector of examined region, c(minlon, maxlon, minlat, maxlat, citylon, citylat)
lon.lat <- get.lon.lat(site)  # can be NULL, default will be given in site.info()

# 'thred.count' for at least how many soundings needed per 1deg lat range
# -> calculate a total thred on total # of soundings given 'lon.lat'
thred.count.per.deg <- 100  # number of soundings per degree
thred.count.per.deg.urban <- 50

# whether to re-search OCO-2 overpasses and output in txtfile
# if FALSE, read timestr from existing txt file;
# always TRUE, if doing first simulation for a new site
searchTF <- F

# whether search for overpasses over urban region,
# defined as city.lat +/- dlat, city.lon +/- dlon
urbanTF <- T; dlon <- 0.5; dlat <- 0.5

# call get.site.info() to get lon.lat and OCO2 overpasses info
# PLEASE add lat lon info in 'get.site.track'
site.info <- get.site.track(site, oco2.ver, oco2.path, searchTF,
  date.range, thred.count.per.deg, lon.lat, urbanTF, dlon, dlat,
  thred.count.per.deg.urban, txtpath)

# get coordinate info and OCO2 track info from result 'site.info'
lon.lat <- site.info$lon.lat
oco2.track <- site.info$oco2.track
print(lon.lat)

# one can further subset 'oco2.track' based on sounding # over near city center
# one can further subset 'oco2.track' based on data quality
# see columns 'qf.count' or 'wl.count' in 'oco2.track'
# e.g., choose overpasses that have 100 soundings with QF == 0, & get reordered
if (urbanTF) oco2.track <- oco2.track %>% filter(tot.urban.count > 200)
if (oco2.ver == 'b7rb') oco2.track <- oco2.track %>% filter(qf.urban.count > 80)
if (oco2.ver == 'b8r') oco2.track <- oco2.track %>% filter(wl.urban.count > 100)

# finally narrow down and get timestr
all.timestr <- oco2.track$timestr[c(2, 3, 8, 9, 10)]

# once you have all timestr, you can choose whether to plot them on maps
# this helps you choose which overpass to simulate first, see 'tt' below
plotTF <- F
if (plotTF) {
  for(t in 1:length(all.timestr)){
  ggmap.obs.xco2(site, timestr = all.timestr[t], oco2.path, lon.lat, workdir,
    plotdir = file.path(workdir, 'plot/ggmap', site))
  ggmap.obs.sif(site, timestr = all.timestr[t], sif.path, lon.lat, workdir,
    plotdir = file.path(workdir, 'plot/ggmap', site))
  }
}

# *** NOW choose the timestr that you would like to work on...
tt <- 2
timestr <- all.timestr[tt]

cat(paste('Working on:', timestr, 'for city/region:', site, '...\n\n'))
cat('Done with choosing cities & overpasses...\n')


#------------------------------ STEP 2 --------------------------------------- #
#### Whether forward/backward, release from a column or a box
columnTF    <- T    # whether a column receptor or fixed receptor
stilt.ver   <- 2    # STILT versions (call different footprint algorithms)
delt        <- 2    # fixed timestep [min]; set = 0 for dynamic timestep
nhrs        <- -72  # number of hours backward (-) or forward (+)

run_trajec  <- T    # run trajec, will overwrite existing 'out'
run_foot    <- F    # run footprint, trajec requires
run_sim     <- F    # calculate simulated XCO2.ff, see STEP 8

run_hor_err <- T    # run traj with hor wind errors/calc trans error (see STEP3)
run_ver_err <- F    # run traj with mixed layer height scaling (see STEP3)
if (run_hor_err) err.level <- c('lower', 'upper')[1]  # which part to be modeled

if (run_trajec) cat('Need to generate trajec...\n')
if (run_foot)   cat('Need to generate footprint...\n\n')

# change to STILT-R version 2 conventions, see validate_varsiwant()
varstrajec <- c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
                'dmas', 'mlht', 'pres', 'sigw', 'tlgr', 'dens')
cat('Done with choosing forward box OR backward column runs...\n')

#------------------------------ STEP 3 --------------------------------------- #
# select receptors --
#### Set model receptors, AGLs and particel numbers ***
# agl can be a vector, meaning releasing particles from several fixed level
# but if trying to represent air column, use columnTF=T, see below

### 1) if release particles from fixed levels
agl    <- 10         # in mAGL
numpar <- 1000       # par for each time window for forward box runs

### 2) SET COLUMN RECEPTORS as a list, if release particles from a column
if (columnTF) {
  # min, middle, max heights for STILT levels, in METERS
  minagl <- 0
  midagl <- 3000  # cutoff level
  maxagl <- 6000
  dh   <- c(100, 500)   # vertical spacing below and above 'midagl', in METERs
  dpar <- c(10, 50, 100, 200, 2500)[3]  # part # per level, 2500 for Brute Force

  # compute the agl list
  agl.lower <- seq(minagl, midagl, dh[1])
  agl.upper <- seq(midagl+dh[2], maxagl, dh[2])

  # seperate release levels, when generating trajec with wind error
  if (run_hor_err & err.level == 'lower') agl <- list(agl.lower)
  if (run_hor_err & err.level == 'upper') agl <- list(agl.upper)
  if (run_hor_err == F) agl <- list(c(agl.lower, agl.upper))
  numpar <- length(unlist(agl)) * dpar   # total number of particles
} # end if columnTF

## eyeball lat range for enhanced XCO2, or grab from available forward-time runs
# place denser receptors during this lat range (40pts in 1deg)
selTF <- T  # whether select receptors; or simulate all soundings
if (selTF) {
  # lat range in deg N for placing denser receptors, required for backward run
  if (site == 'Riyadh')    peak.lat <- c(24, 25)
  if (site == 'Jerusalem') peak.lat <- c(31.75, 32.25)
  if (site == 'Cairo')     peak.lat <- c(29.5, 30.5)
  if (site == 'Phoenix')   peak.lat <- c(33.0, 34.0)
  if (site == 'Baghdad')   peak.lat <- c(32.5, 33.5)

  # number of points to aggregate within 1deg over small/large enhancements,
  # i.e., over background/enhancements, binwidth will be 1deg/num
  num.bg   <- 20   # e.g., every 20 pts in 1 deg
  num.peak <- 40   # e.g., every 40 pts in 1 deg

  # recp.indx: how to pick receptors from all screened soundings (QF = 0)
  recp.indx <- c(seq(lon.lat[3],  peak.lat[1], 1/num.bg),
    seq(peak.lat[1], peak.lat[2], 1/num.peak),
    seq(peak.lat[1], lon.lat[4],  1/num.bg))
} else {
  recp.indx <- NULL
}

# whether to subset receptors when debugging
recp.num <- NULL     # can be a number for max num of receptors
find.lat <- NULL     # for debug or test, model one sounding

# for generating trajec with horizontal error compoennt,
# use the same lat.lon from original trajec, DW, 07/31/2018
if (run_hor_err)
  trajpath <- file.path(homedir, 'lin-group5/wde/github/stilt/workdir',
    paste0('out_', timestr, '_72hrs_v8'), 'by-id')

# select satellite soundings, plotTF for whether plotting OCO-2 observed XCO2
recp.info <- get.recp.info(timestr = timestr, oco2.path = oco2.path,
  oco2.ver = oco2.ver, lon.lat = lon.lat, selTF = selTF,
  recp.indx = recp.indx, recp.num = recp.num, find.lat = find.lat, agl = agl,
  plotTF = F, run_trajec = run_trajec, run_hor_err = run_hor_err,
  trajpath = trajpath, stilt.ver = stilt.ver)

nrecp <- nrow(recp.info)
print(nrecp)
cat(paste('Done with receptor setup...\n'))

#------------------------------ STEP 4 --------------------------------------- #
# path for the ARL format of WRF and GDAS
# simulation_step() will find corresponding met files
met        <- c('1km', 'gdas1', 'gdas0p5')[3]  # choose met fields
met.path   <- file.path(homedir, 'u0947337', met)
met.num    <- 1                                # min number of met files needed
met.format <- '%Y%m%d_gdas0p5'                 # met file name convention

# one can link to other direcetory that store trajec,
# but need to have the same directory structure, including by-id, footprint...
outdir     <- file.path(workdir, 'out')  # path for storing trajec, foot

#### Whether obtaining wind errors, transport error component
# require wind error comparisons stored in txtfile *****
if (run_hor_err) {
  cat('+++ horizontal wind error component +++\n')

  # intput correlation lengthscale (in m) and timescales (in mins)
  # correlation timescale, horizontal and vertical lengthscales
  if (met == 'gdas0p5') {TLuverr <- 1*60; zcoruverr <- 600; horcoruverr <- 40}
  if (met == 'gdas1') {TLuverr <- 2.4*60; zcoruverr <- 700; horcoruverr <- 97}

  ## add errors, mainly siguverr, create a subroutine to compute siguverr
  # from model-data wind comparisons
  err.path <- file.path(homedir, 'lin-group5/wde/input_data/wind_err')
  err.path <- file.path(err.path, tolower(site), tolower(met))

  # call get.SIGUVERR() to interpolate most near-field wind errors
  err.info <- get.siguverr(site, timestr, errpath = err.path, nfTF = F,
    forwardTF = F, lon.lat, nhrs, agl)

  if (is.null(err.info)) {
    cat('no wind error found; make consevative assumption of siguverr...\n')
    # make a conservative assumption about the wind error, for the Middle East
    siguverr <- 1.8  # < 2 m/s for GDAS 1deg, based on Wu et al., GMDD

  } else {
    met.rad  <- err.info[[1]]
    siguverr <- as.numeric(err.info[[2]][1])    # SD in wind errors
    u.bias   <- as.numeric(err.info[[2]][2])
    v.bias   <- as.numeric(err.info[[2]][3])
    cat(paste('u.bias:', signif(u.bias,3), 'm/s; v.bias:', signif(v.bias,3), 'm/s\n'))

  }  # end if is.null(err.info)
  cat(paste('SIGUVERR:', signif(siguverr, 3), 'm/s..\n'))

} else {  # if no wine error component used
  cat('NO horizontal wind error component for generating trajec...\n')
  siguverr    <- NA
  TLuverr     <- NA
  horcoruverr <- NA
  zcoruverr   <- NA
}  # end if run_hor_err

# no error covariance on PBL heights used for now
# but one can assign values as below
sigzierr    <- NA
TLzierr     <- NA
horcorzierr <- NA

### Besides horizontal wind error, do we want to account for PBL height?
# add vertical trans error via ziscale *****
if (run_ver_err) {
  zicontroltf <- 1              # 0 for FALSE; 1 for scaling, TRUE
  ziscale     <- rep(list(rep(0.8, 24)), nrecp)  # create as list
  # 1st # for scaling factor; 2nd # for # of hours (always use abs())
  cat(paste('+++ Mixed layer height scaling of', ziscale[[1]][1],
    'when generating trajec +++\n'))

} else {
  cat('NO Mixed layer height scaling ...\n')
  zicontroltf <- 0
  ziscale     <- 1.0
} # end if run_ver_err
cat('Done with choosing met & inputting wind errors...\n')

#------------------------------ STEP 5 --------------------------------------- #
#### Settings for generating footprint maps, no need for getting background
## SET spatial domains and resolution for calculating footprints
foot.res <- 1/120  # footprint resolution, 1km for ODIAC

# these variables will determine resoluation and spatial domain of footprint
# 20x20 degree domain around the city center
foot.info <- data.frame(
  xmn = round(lon.lat[5]) - 10, xmx = round(lon.lat[5]) + 10,
  ymn = round(lon.lat[6]) - 10, ymx = round(lon.lat[6]) + 10,
  xres = foot.res, yres = foot.res
)
# OR customize foot domain, in deg E and deg N
foot.info <- data.frame(xmn = 30, xmx = 50, ymn = 15, ymx = 35, xres = foot.res,
  yres = foot.res)
print(foot.info)

## whether weighted footprint by AK and PW for column simulations
if (columnTF) {
  ak.wgt <- T;  pwf.wgt <- T
} else {  # no weighting needed for fixed receptor simulations
  ak.wgt <- NA; pwf.wgt <- NA
}

# 1st-order correction on dmass for footprint if using STILTv1
dmassTF <- F

# other footprint parameters if using STILTv2
hnf_plume      <- T  # whether turn on hyper near-field (HNP) for mising hgts
smooth_factor  <- 1  # Gaussian smooth factor, 0 to disable
time_integrate <- T  # whether integrate footprint along time
projection     <- '+proj=longlat'
cat('Done with footprint setup...\n')


#------------------------------ STEP 6 --------------------------------------- #
#### !!! NO NEED TO CHANGE ANYTHING LISTED BELOW -->
# create a namelist including all variables
# namelist required for generating trajec
namelist <- list(agl = agl, ak.wgt = ak.wgt, delt = delt, dmassTF = dmassTF,
  dpar = dpar, foot.info = foot.info, hnf_plume = hnf_plume, homedir = homedir,
  horcoruverr = horcoruverr, horcorzierr = horcorzierr, lon.lat = lon.lat,
  met = met, met.format = met.format, met.num = met.num, met.path = met.path,
  nhrs = nhrs, numpar = numpar, outdir = outdir, oco2.path = oco2.path,
  projection = projection, pwf.wgt = pwf.wgt, recp.info = recp.info,
  run_foot = run_foot, run_sim = run_sim, run_trajec = run_trajec,
  run_hor_err = run_hor_err, run_ver_err = run_ver_err, siguverr = siguverr,
  sigzierr = sigzierr, site = site, smooth_factor = smooth_factor,
  stilt.ver = stilt.ver, time_integrate = time_integrate, timestr = timestr,
  TLuverr = TLuverr, TLzierr = TLzierr, varstrajec = varstrajec,
  windowTF = windowTF, workdir = workdir, zicontroltf = zicontroltf,
  ziscale = ziscale, zcoruverr = zcoruverr)


#------------------------------ STEP 7 --------------------------------------- #
## run trajec/footprint/simulations for backward simulations --
if (run_trajec | run_foot) {    ## if running trajec or footprint

  ## use Ben's algorithm for parallel simulation settings
  n_nodes  <- 5
  n_cores  <- ceiling(nrecp/n_nodes)
  job.time <- '24:00:00'
  slurm    <- n_nodes > 1
  namelist$slurm_options <- list(time = job.time, account = 'lin-kp',
    partition = 'lin-kp')

  # time allowed for running hymodelc before forced terminations
  timeout  <- 2 * 60 * 60  # in sec
  namelist <- c(namelist, n_cores = n_cores, n_nodes = n_nodes, slurm = slurm,
    timeout = timeout)
  cat('Done with creating namelist...\n')

  # call run_stilt_mod(), start running trajec and foot
  run.stiltv2(namelist = namelist)

} else if (run_sim) {    ## if running trajec or footprint

  # requires trajec and footprints ready for running this sections:
  # Simulate XCO2.ff using ODIAC emissions, DW, 06/04/2018
  # grab footprint info
  foot.path <- file.path(workdir, 'out/by-id')
  foot.file <- file.path(foot.path, list.files(foot.path, pattern = 'foot.nc',
    recursive = T))

  tmp.foot <- raster(foot.file[1]); foot.extent <- extent(tmp.foot)
  cat('Done reading footprint..\n')

  # txt file name for outputting model results
  txtfile <- file.path(workdir, paste0(timestr, '_', site, '_XCO2ff_',
    abs(nhrs), 'hrs_', dpar, 'dpar_sf', smooth_factor, '.txt'))

  # before simulations, subset emissions and convert tif format to nc format
  vname <- '2017'; tiff.path <- file.path(homedir,
    paste0('lin-group2/group_data/ODIAC/ODIAC', vname),
    substr(timestr, 1,4))  # tif file from ODIAC website

  # call tif2nc.odiacv2() to subset and get emiss file name
  # 'store.path' is the path for outputting emissions
  cat('Start reading and subsetting emissions that match foot...\n')
  emiss.file <- tif2nc.odiacv3(site, timestr, vname, workdir, foot.extent,
    store.path = file.path(workdir, 'in', 'ODIAC'), tiff.path, gzTF = F)

  # reformatted ODIAC emissions file name should include 'YYYYMM'
  # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
  cat('Start XCO2.ff simulations...\n')
  receptor <- foot.odiacv3(foot.file, emiss.file, workdir, txtfile, lon.lat,
    plotTF = F)

  #ff2 <- read.table(txtfile, sep = ',', header = T)
  #l1 <- ggplot() + geom_point(data = receptor, aes(lat, xco2.ff), colour = 'red')

  ### add latitude integrations--
  library(zoo)
  auc <- diff(receptor$lat) * rollmean(receptor$xco2.ff, 2)
  xco2.ff.int <- sum(auc[auc > 0])
  cat(paste('Lat-integrated XCO2.ff:', signif(xco2.ff.int, 3), 'ppm\n'))
} # end if run traj/foot/sim



#------------------------------ STEP 8 --------------------------------------- #
### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
# requires two sets of trajectories before running the following section:
if (run_hor_err) {

  # for ffco2 emission path and files
  vname <- '2017'
  tiff.path <- file.path(homedir, paste0('lin-group2/group_data/ODIAC/ODIAC',
    vname), substr(timestr, 1,4))  # tif file from ODIAC website
  foot.extent <- extent(as.numeric(foot.info[1:4]))
  emiss.file <- tif2nc.odiacv3(site, timestr, vname, workdir, foot.extent,
    store.path = file.path(workdir, 'in', 'ODIAC'), tiff.path, gzTF = F)

  # load two paths for trajec before and after perturbations
  traj.path1 <- file.path(homedir, 'lin-group5/wde/github/cp_trajecfoot/out_2014122910_100dpar/by-id')
  traj.path2 <- file.path(workdir, 'out/by-id')

  # add CT paths and files, DW, 07/28/2018
  ct.ver <- '2016-1'; if (timestr >= 20160101) ct.version <- '2017'
  ct.path <- '/uufs/chpc.utah.edu/common/home/lin-group2/group_data/CT-NRT'
  ctflux.path <- file.path(ct.path, paste0('v', ct.ver), 'fluxes/optimized')
  ctmole.path <- file.path(ct.path, paste0('v', ct.ver), 'molefractions/co2_total')

  # path for storing CO2 statistics
  store.path <- file.path(output.path, 'trans_err')

  # txt file name for outputting model results
  txtfile <- file.path(workdir, paste0(timestr, '_', site, '_trans_err_',
    met, '.txt'))

  namelist <- c(namelist, emiss.file = emiss.file, traj.path1 = traj.path1,
    traj.path2 = traj.path2, ct.ver = ct.ver, ctflux.path = ctflux.path,
    ctmole.path = ctmole.path, store.path = store.path, txtfile = txtfile)

  # call cal.trans.err to estimate trans err
  trans.err <- cal.trans.err(namelist)
} # end if run_hor_err

# end of script
