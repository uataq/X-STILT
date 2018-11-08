#' create a namelist for running X-STILT trajectories
#' @author Dien Wu, 04/18/2018
#' -----------------------------------------------------------------------------

#' @flags: 
#' 1. run_trajec, run_foot: run trajectory or footprint
#' 1) if run_hor_err == T or run_ver_err == T
#'    X-STILT generates trajectory with perturbations from wind or PBL errors
#'    see details in STEP 4. 
#' 2) if run_emiss_err is turned on
#'    X-STILT generates footprint with horizontal resolution of 0.1deg, 
#'    rather than 1km, because the absolute emission error has res of 0.1deg. 
#'    see details in STEP 8. 

#' 2. run_sim: run simulations (requires at least trajec done)
#' 1) if no other flags is turned on, 
#'    X-STILT models the XCO2 enhancements based on FFCO2 emission from 1km ODIAC
#' 2) if run_hor_err is turned on, 
#'    X-STILT models XCO2 error due to hortizontal transport error, See STEP 8. 
#' -----------------------------------------------------------------------------

#' @updates:
#' now build upon Ben's STILT-R version 2 codes, DW, 05/23/2018
#' !!! need to clear up codes for forward run, based on Ben's parallel computing
#' Add plotting scripts for footprint and XCO2 enhancements, DW, 06/28/2018
#
#' Add STILTv1 dependencies and work well with existing framework;
#'   add dmassTF for mass violation corrections, DW, 07/17/2018
#
#' Add horizontal transport error module with flag 'run_hor_err'
#'   with additional subroutines in /src/oco2-xstilt/trans_err, DW, 07/23/2018
#
#' Add vertical trans error module, with flag 'run_ver_err', add ziscale
#'   when scaling PBL heights and break trans error part (step 3) into
#'   horizontal wind error and vertical trans error (via PBL perturb), 
#'   DW, 07/25/2018
#
#' simplify main code and use 'get.uverr()' to get horizontal wind error stat,
#'   DW, 08/31/2018
#
#' add emission uncertainty, DW, 09/06/2018
#' update code according to v9 changes, DW, 10/19/2018
#
#' Optimize horizontal trans error as parallel computing to save time, DW, 10/21/2018
#' -----------------------------------------------------------------------------

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group7/wde/github/XSTILT') # current dir
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions

# insert your API for the use of ggplot and ggmap
api.key <- ''
register_google(key = api.key)


#------------------------------ STEP 1 --------------------------------------- #
# input dlat, dlon to get spatial domain around city center
site <- 'Riyadh'   # choose a city
lon.lat <- get.lon.lat(site = site, dlon = 1, dlat = 1.5)

# required paths
oco2.ver   <- c('b7rb', 'b8r', 'b9r')[1]  # OCO-2 version
input.path <- file.path(homedir, 'lin-group7/wde/input_data')
oco2.path  <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))
sif.path   <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_SIF_', oco2.ver))
raob.path  <- file.path(input.path, 'RAOB/middle.east/riyadh')  # radiosonde

# emission prior, ODIAC version and paths, original from ODIAC website in tiff
vname <- c('2016', '2017')[2]  # ODIAC version
tiff.path <- file.path(input.path, 'ODIAC', paste0('ODIAC', vname))  

# path for storing plot or CO2 statistics
store.path <- file.path(workdir, gsub(' ', '', lon.lat$regid), site)
dir.create(store.path, showWarnings = F, recursive = T)

# whether search for overpasses over urban region,
# defined as city.lat +/- dlat, city.lon +/- dlon
urbanTF <- T; dlon.urban <- 0.5; dlat.urban <- 0.5

# path for storing overpass sampling info by 'get,site.track'
txt.path <- file.path(input.path, 'OCO-2/overpass_city') 

# call get.site.info() to get lon.lat and OCO2 overpasses info
oco2.track <- get.site.track(site, oco2.ver, oco2.path, searchTF = F, 
                             date.range = c('20140101', '20181231'), 
                             thred.count.per.deg = 100, lon.lat, 
                             urbanTF, dlon.urban, dlat.urban, 
                             thred.count.per.deg.urban = 50, txt.path)

# one can further subset 'oco2.track' based on sounding # or data quality
# over entire lon.lat domain or near city center
# see columns 'qf.count' or 'wl.count' in 'oco2.track'
oco2.track <- oco2.track %>% filter(qf.urban.count > 50)

rmTF <- F  # whether to remove summtertime tracks
if (rmTF) {
   if (lon.lat$citylat > 0) {   
    # northern Hemi
    oco2.track <- oco2.track %>% filter(substr(timestr, 5, 6) < '05' | 
                                        substr(timestr, 5, 6) > '08')
  } else {  
    # southern Hemi
    oco2.track <- oco2.track %>% filter(substr(timestr, 5, 6) < '12' & 
                                        substr(timestr, 5, 6) > '03')
  } # end if Hemisphere
}  # end if rmTF

# finally narrow down and get timestr
all.timestr <- oco2.track$timestr[c(2, 3, 5, 6, 7)]         # v7, Riyadh
print(all.timestr)

# once you have all timestr, you can choose whether to plot them on maps
# this helps you choose which overpass to simulate first, see 'tt' below
plotTF <- F
if (plotTF) {
  plotdir <- file.path(store.path, 'plot')
  dir.create(plotdir, showWarnings = F, recursive = T)

  # zoom scales 1 to 10, larger the value, more zoomed in
  # qfTF = T for only plotting data with QF = 0
  for (t in 1:length(all.timestr)) {
    ggmap.obs.xco2(site, all.timestr[t], oco2.ver, oco2.path, lon.lat, workdir, 
                   plotdir, zoom = 8, qfTF = F, box.dlat = dlat.urban, 
                   box.dlon = dlon.urban)
    ggmap.obs.sif(site, all.timestr[t], sif.path, lon.lat, workdir, plotdir, 
                   zoom = 8, box.dlon = dlon.urban, box.dlat = dlat.urban)
  } # end for t
}  # end if plotTF


# *** NOW choose the timestr that you would like to work on...
tt <- 2
timestr <- all.timestr[tt]
cat(paste('Working on:', timestr, 'for city/region:', site, '...\n\n'))
cat('Done with choosing cities & overpasses...\n')


#------------------------------ STEP 2 --------------------------------------- #
# T:rerun hymodelc, even if particle location object found
# F:re-use previously calculated particle location object
run_trajec <- T      # whether to generate trajec
run_foot   <- T      # whether to generate footprint
columnTF   <- T      # whether a column receptor or fixed receptor
if (run_trajec) cat('Need to generate trajec...\n')
if (run_foot)   cat('Need to generate footprint...\n\n')

# whether to perform XCO2 and its error simulations
run_sim       <- F    # whether to run analysis, see details in STEP 8
run_hor_err   <- F    # T: set parameters in STEP 4
run_ver_err   <- F    # T: set parameters in STEP 4
run_emiss_err <- F    # T: get XCO2 error due to prior emiss err, see STEP 8

stilt.ver <- 2    # STILT versions (call different footprint algorithms)
delt      <- 2    # fixed timestep [min]; set = 0 for dynamic timestep
nhrs      <- -72  # number of hours backward (-) or forward (+)

# where to grab or store trajec, DW, 07/31/2018
# ourput directory for storing traj, default convention
# one can change the path
outdir <- file.path(store.path, paste('out', timestr, site, sep = '_'))

# change to Ben's definitions,  see validate_varsiwant()
varstrajec <- c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
                'dmas', 'mlht', 'temp', 'pres', 'sigw', 'tlgr', 'dens')
cat('Done with setting flags...\n')


#------------------------------ STEP 3 --------------------------------------- #
# select receptors --
#### Set model receptors, AGLs and particel numbers ***
# for backward fixed-level runs OR forward fixed-level runs
# agl can be a vector, meaning releasing particles from several fixed level
# but if trying to represent air column, use columnTF=T, see below

### 1) if release particles from fixed levels
agl    <- c(10, 500, 2000, 3000)[4]        # in mAGL
numpar <- 100       # par for each time window for forward box runs
dpar   <- numpar

### 2) SET COLUMN RECEPTORS as a list, if release particles from a column
if (columnTF) {

  # min, middle, max heights for STILT levels, in METERS
  minagl <- 0; midagl <- 3000; maxagl <- 6000

  # vertical spacing below and above cutoff level 'midagl', in METERS
  dh <- c(100, 500)

  # particle numbers per level, 2500 for the Brute Force test
  dpar <- c(10, 50, 100, 200, 2500)[3]

  # compute the agl list
  agl  <- list(c(seq(minagl, midagl, dh[1]), seq(midagl+dh[2], maxagl, dh[2])))
  nlev <- length(unlist(agl))
  numpar <- nlev * dpar   # total number of particles
}  # end if columnTF

## eyeball lat range for enhanced XCO2, or grab from available forward-time runs
# place denser receptors during this lat range (40pts in 1deg)
selTF <- T  # whether to select receptors; or simulate all soundings
if (selTF) {
  # lat range in deg N for placing denser receptors, required for backward run
  peak.lat <- c(lon.lat$citylat - dlat.urban, lon.lat$citylat + dlat.urban)

  # number of points to aggregate within 1deg over small/large enhancements,
  # i.e., over background/enhancements, binwidth will be 1deg/num
  num.bg   <- 20   # e.g., every 20 pts in 1 deg
  num.peak <- 40   # e.g., every 40 pts in 1 deg

  # recp.indx: how to pick receptors from all screened soundings (QF = 0)
  recp.indx <- c(seq(lon.lat$minlat,  peak.lat[1], 1/num.bg),
                 seq(peak.lat[1], peak.lat[2], 1/num.peak),
                 seq(peak.lat[1], lon.lat$maxlat,  1/num.bg))

} else { recp.indx <- NULL }

# whether to subset receptors when debugging; if no subset, insert NULL
recp.num <- NULL     # can be a number for max num of receptors
find.lat <- NULL     # for debug or test, model one sounding

# select satellite soundings, plotTF for whether plotting OCO-2 observed XCO2
# data.filter: use WL or QF to filter data, default is QF = 0
recp.info <- get.recp.info(timestr, oco2.ver, oco2.path, lon.lat, selTF, 
                           recp.indx, recp.num, find.lat, agl, plotTF = F, 
                           trajpath = file.path(outdir, 'by-id'), run_trajec,
                           stilt.ver, data.filter = c('QF', 0))
nrecp <- nrow(recp.info)
cat(paste('Done with receptor setup...total', nrecp, 'receptors..\n'))


#------------------------------ STEP 4 --------------------------------------- #
## path for the ARL format of WRF and GDAS
# simulation_step() will find corresponding met files
met        <- 'gdas0p5' # choose met fields
met.path   <- file.path(homedir, 'u0947337', met)  # path of met fields
met.format <- paste0('%Y%m%d_', met)           # met file name convention
met.num    <- 1                                # min number of met files needed

## get horizontal transport error component if run_hor_err = T
# path for outputting wind error stats
err.path <- file.path(store.path, 'wind_err')  
hor.err  <- get.uverr(run_hor_err, site, timestr, workdir, overwrite = F,
                      raob.path, raob.format = 'fsl', nhrs, met, met.path, 
                      met.format, lon.lat, agl, err.path)

## get vertical transport error component if run_ver_err = T
# set zisf = 1 if run_ver_err = F
zisf <- c(0.6, 0.8, 1.0, 1.2, 1.4)[3]; if (!run_ver_err) zisf <- 1.0
pbl.err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = zisf)

# if calculating XCO2 error due to emiss error, need files for EDGAR and FFDAS
# rearrange by DW, 10/21/2018
if (run_emiss_err) { 
  edgar.file <- file.path(homedir,
    'lin-group2/group_data/EDGAR/CO2_2000-2008/v42_CO2_2008_TOT.0.1x0.1.nc')
  ffdas.path <- file.path(input.path, 'anthro_inventories/FFDAS/')
  ffdas.file <- list.files(path = ffdas.path, pattern = 'totals')
  ffdas.file <- file.path(ffdas.path, ffdas.file[grep('2008', ffdas.file)])
} else { edgar.file = NULL; ffdas.file = NULL } 
# end if run_emiss_err

# if calculating XCO2 error due to horizontal wind error, need biospheric input
# add CT paths and files, DW, 07/28/2018
if (run_hor_err) {
  ct.ver  <- ifelse(substr(timestr, 1, 4) >= '2016', 'v2017', 'v2016-1')
  ct.path <- file.path(input.path, 'CT-NRT', ct.ver)
  ctflux.path <- file.path(ct.path, 'fluxes/optimized')
  ctmole.path <- file.path(ct.path, 'molefractions/co2_total')
} # end of run_hor_err

cat('Done with choosing met & inputting parameters for error estimates...\n')


#------------------------------ STEP 5 --------------------------------------- #
#### Settings for generating footprint maps
## SET spatial domains and resolution for calculating footprints
foot.res <- 1/120  # footprint resolution, 1km for ODIAC
if (run_emiss_err) foot.res <- 1/10  # for emiss error, generate 0.1deg foot
#foot.res <- 1  # for generating foot that matches CarbonTracker

# these variables will determine resoluation and spatial domain of footprint
# 20x20 degree domain around the city center
foot.info <- data.frame(xmn = round(lon.lat$minlon) - 10, 
                        xmx = round(lon.lat$minlon) + 10,
                        ymn = round(lon.lat$minlat) - 10, 
                        ymx = round(lon.lat$minlat) + 10,
                        xres = foot.res, yres = foot.res)

### OR customize foot domain, in deg E and deg N
#foot.info <- data.frame(xmn = 30, xmx = 50, ymn = 15, ymx = 35, 
#  xres = foot.res, yres = foot.res)
foot.extent <- extent(foot.info$xmn, foot.info$xmx, foot.info$ymn, foot.info$ymx)
print(foot.info)

## whether weighted footprint by AK and PW for column simulations
if (columnTF) {
  ak.wgt <- T; pwf.wgt  <- T
} else {  # no weighting needed for fixed receptor simulations
  ak.wgt <- NA; pwf.wgt <- NA
}

dmassTF <- F        # 1st-order correction on dmass for footprint using STILTv1
# other footprint parameters using STILTv2:
hnf_plume      <- T  # whether turn on hyper near-field (HNP) for mising hgts
smooth_factor  <- 1  # Gaussian smooth factor, 0 to disable
time_integrate <- T  # whether integrate footprint along time, T, no hourly foot
projection     <- '+proj=longlat'
cat('Done with footprint setup...\n\n')


#------------------------------ STEP 7 --------------------------------------- #
## if running trajec or footprint
if (run_trajec | run_foot) {

  ## use SLURM for parallel simulation settings
  n_nodes  <- 8
  n_cores  <- ceiling(nrecp/n_nodes)

  # time allowed for running hymodelc before forced terminations
  timeout  <- 24 * 60 * 60  # in sec
  job.time <- '24:00:00'    # total job time
  slurm    <- n_nodes > 0
  slurm_options <- list(time = job.time, account = 'lin-kp', partition = 'lin-kp')

  # create a namelist including all variables
  # namelist required for generating trajec
  namelist <- list(agl = agl, ak.wgt = ak.wgt, delt = delt, dmassTF = dmassTF,
                   dpar = dpar, foot.info = foot.info, hnf_plume = hnf_plume, 
                   hor.err = hor.err, lon.lat = lon.lat, met = met, 
                   met.format = met.format, met.num = met.num, 
                   met.path = met.path, n_cores = n_cores, nhrs = nhrs,
                   n_nodes = n_nodes, numpar = numpar, outdir = outdir, 
                   oco2.path = oco2.path, pbl.err = pbl.err, 
                   projection = projection, pwf.wgt = pwf.wgt, 
                   recp.info = recp.info, run_foot = run_foot, 
                   run_hor_err = run_hor_err, run_trajec = run_trajec, site = site, 
                   slurm = slurm, slurm_options = slurm_options, 
                   smooth_factor = smooth_factor, stilt.ver = stilt.ver, 
                   time_integrate = time_integrate, timeout = timeout, 
                   timestr = timestr, varstrajec = varstrajec, workdir = workdir)         
  cat('Done with creating namelist...\n')

  # call run_stiltv2() to start running trajec and foot
  run.stiltv2(namelist)
} # end if run trajec or foot


#------------------------------ STEP 8 --------------------------------------- #
## if calculating XCO2 and its error (need trajec and footprint ready)
if (run_trajec * run_foot == F) {

  #------------------------  Horizontal trans error -------------------------- #
  ### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
  # requires two sets of trajectories before running the following section:
  if (run_hor_err) { # this does not need footprint

    # call function cal.trans.err() to estimate trans err
    if (run_sim) {
      library(zoo)
      cat('Start simulations of XCO2 error due to horizontal trans err...\n')
      # get actual ppm error, need to have error statistics ready (see above)
      receptor <- cal.trans.err(site, timestr, workdir, outdir, store.path, met)

      # add latitudinally integrated horizontal trans error, DW, 10/22/2018
      # need to run rolling mean on XCO2 hor trans error (** operate on variance)
      if (!is.null(receptor)) {
        w <- abs(diff(receptor$lat))                     # diff in lat as weights
        x <- as.data.frame(rollmean(receptor, 2))        # rollmean on all results
        x$sd <- sqrt(rollmean(receptor$xco2.trans^2, 2)) # replace SD

        # error variance covariance matrix and final integrated SD in ppm
        err.var.cov <- cal.var.cov(sd = x$sd, L = 40000, w, x, type = 'hor')
        tot.hor.err <- sqrt(sum(err.var.cov, na.rm = T)) 
      }  # end if

    } else {  # run horizontal trans error statistics ------------------------ #

      ## use SLURM for parallel simulation settings
      # time (in sec) allowed for running hymodelc before forced terminations
      n_nodes  <- 6
      n_cores  <- ceiling(nrecp / n_nodes)
      slurm    <- n_nodes > 1
      job.time <- '01:00:00'    # total job time
      slurm_options <- list(time = job.time, account = 'lin-kp', partition = 'lin-kp')

      # need ODIAC emissions to calculate trans error
      emiss.file <- tif2nc.odiacv3(site, timestr, vname, workdir, foot.extent,
                                   tiff.path, gzTF = F)

      # get error statistics via calling 'cal.trajfoot.stat()'
      stilt_apply(X = 1:nrecp, FUN = cal.trajfoot.stat, slurm = slurm, 
                  slurm_options = slurm_options, n_nodes = n_nodes, 
                  n_cores = n_cores, workdir = workdir, outdir = outdir, 
                  emiss.file = emiss.file, met = met, dpar = dpar, 
                  ct.ver = ct.ver, ctflux.path = ctflux.path, 
                  ctmole.path = ctmole.path, r_run_time = recp.info$run_time, 
                  r_lati = recp.info$lati, r_long = recp.info$long, 
                  r_zagl = recp.info$zagl)
      q('no')
    } # end if run_sim
  
  } else {
  
    #--------------------- XCO2 or emiss error simulations ---------------------- #
    # requires trajec and footprints ready for running things below, DW, 06/04/2018  
    # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
    cat('Start simulations of XCO2.ff or its error due to emiss err...\n')
    receptor <- run.xco2.sim(site, timestr, vname, tiff.path, outdir, foot.res, 
                             workdir, store.path, nhrs, dpar, smooth_factor, 
                             zisf, oco2.ver, met, lon.lat, run_emiss_err, 
                             edgar.file, ffdas.file, plotTF = F, writeTF = T)
    
    # add latitudinally integrated emiss error, DW, 10/22/2018
    # need to run rolling mean on XCO2 emiss err (** operate on variance)
    if (!is.null(receptor)) {
      if (run_emiss_err) {
        w <- abs(diff(receptor$lat))                      # diff in lat as weights
        x <- as.data.frame(rollmean(receptor, 2))         # rollmean all results
        x$sd <- sqrt(rollmean(receptor$xco2.ff.err^2, 2)) # replace with correct SD

        # L for error covariance of horizontal transport error (ppm) in meters, 
        # not error covariance of wind error (m/s)
        err.var.cov <- cal.var.cov(sd = x$sd, L = 40000, w, x, type = 'hor')
        tot.hor.err <- sqrt(sum(err.var.cov, na.rm = T))  # final integrated SD in ppm

      } else {
        # add lat-int XCO2 signals
        w <- abs(diff(receptor$lat))                    # diff in lat as weights
        x <- as.data.frame(rollmean(receptor, 2))  
        tot.sig <- sum(w * x$xco2.ff)
      } # end if run_emiss_err
    } # end if is.null

  } # end if run_hor_err
} # end if run_sim

##### end of script
