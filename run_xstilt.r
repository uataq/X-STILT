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
#'    see details in STEP 7. 

#' 2. run_sim: run simulations (requires at least trajec done)
#' 1) if no other flags is turned on, 
#'    X-STILT models the XCO2 enhancements based on FFCO2 emission from 1km ODIAC
#' 2) if run_hor_err is turned on, 
#'    X-STILT models XCO2 error due to hortizontal transport error, See STEP 7. 
#' -----------------------------------------------------------------------------

#' @updates:
#' now build upon Ben Fasoli's STILT-R version 2 codes, DW, 05/23/2018
#' !!! need to clear up codes for forward run, based on Ben's parallel computing
#' Add plotting scripts for footprint and XCO2 enhancements, DW, 06/28/2018
#
#' Add horizontal transport error module with flag 'run_hor_err'
#'   with additional subroutines in /src/oco2-xstilt/trans_err, DW, 07/23/2018
#
#' Add vertical trans error module, with flag 'run_ver_err', add ziscale
#'   when scaling PBL heights and break trans error part (step 4) into
#'   horizontal wind error and vertical trans error (via PBL perturb), 
#'   DW, 07/25/2018
#
#' simplify main code and use 'get.uverr()' to get horizontal wind error stat,
#'   DW, 08/31/2018
#
#' add emission uncertainty, DW, 09/06/2018
#' update code according to v9 changes, DW, 10/19/2018
#
#' Optimize horizontal trans error using parallel computing, DW, 10/21/2018
#'
#' --------------- Reconstruct X-STILT codes, DW, BF, 01/25-29/2019 -----------
#' remove STILTv1 (Lin et al., 2003) from X-STILT;
#' separate its modification from STILT; 
#' now STILTv2 (Fasoli et al., 2018) is a SUBMODULE of X-STILT;
#' customize before_trajec and before_footprint functions to mutate `output`. 
#'
#' move the section of estimating trajec-level CO2 from main script to 
#'    before_footprint_xstilt.r. So, each core will work on calculating the 
#'    required trans error statistics (since it takes time).  DW, 01/29/2019
#' to simplify the main script, remove calculations of lat-int signals or err, 
#'    DW, 01/29/2019 
#' -----------------------------------------------------------------------------


#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group7/wde/X-STILT') # current dir
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions

# Please insert your API for the use of ggplot and ggmap
api.key <- ''
register_google(key = api.key)


#------------------------------ STEP 1 --------------------------------------- #
### 1) input dlat, dlon to get spatial domain around city center
site <- 'Riyadh'   # choose a city
lon.lat <- get.lon.lat(site = site, dlon = 1, dlat = 1.5)


### 2) required paths for input datasets
# e.g., OCO2 XCO2, SIF, NOAA RAOB, ODAIC emission (1km tiff files)
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
oco2.ver    <- c('b7rb', 'b8r', 'b9r')[3]           # OCO-2 version
oco2.path   <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))
sif.path    <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_SIF_', oco2.ver))
raob.path   <- file.path(input.path, 'RAOB', site)  # NOAA radiosonde
odiac.vname <- c('2016', '2017', '2018')[3]         # ODIAC version
tiff.path   <- file.path(input.path, 'ODIAC', paste0('ODIAC', odiac.vname))  

## path for storing output or anything related to trans err
store.path <- file.path(homedir, 'lin-group7/wde/output', site)
err.path   <- file.path(store.path, 'wind_err')  
dir.create(store.path, showWarnings = F, recursive = T)
dir.create(err.path, showWarnings = F, recursive = T)


### 3) call get.site.track() to get lon.lat and OCO2 overpasses info
# txt.path, path for storing output from 'get.site.track()'
txt.path   <- file.path(input.path, 'OCO-2/overpass_city') 

# whether to search for overpasses over urban region,
# defined as city.lat +/- dlat, city.lon +/- dlon
urbanTF <- T; dlon.urban <- 0.5; dlat.urban <- 0.5

# rmTF = T, for removing overpasses during hemispheric growing seasons
oco2.track <- get.site.track(site, oco2.ver, oco2.path, searchTF = F, 
                             date.range = c('20140101', '20181231'), 
                             thred.count.per.deg = 100, lon.lat, 
                             urbanTF, dlon.urban, dlat.urban, 
                             thred.count.per.deg.urban = 50, txt.path, rmTF = T)

# one can further subset 'oco2.track' based on sounding # or data quality
# over entire lon.lat domain or near city center
# see columns 'qf.count' or 'wl.count' in 'oco2.track'
oco2.track <- oco2.track %>% filter(qf.urban.count > 50)


### 4) finally narrow down and get timestr
all.timestr <- oco2.track$timestr[c(2, 3, 8, 9, 10)]         # v7, Riyadh
print(all.timestr)

# whether to plot them on maps, plotTF = T/F,
# this helps you choose which overpass to simulate, see 'tt' below
ggmap.obs.info(plotTF = F, store.path, all.timestr, oco2.ver, oco2.path, lon.lat, 
               workdir, dlat.urban, dlon.urban)

### 5) *** NOW choose the timestr that you'd like to work on...
tt <- 5
timestr <- all.timestr[tt]
cat(paste('Working on:', timestr, 'for city/region:', site, '...\n\n'))
cat('Done with choosing cities & overpasses...\n')


#------------------------------ STEP 2 --------------------------------------- #
# T:rerun hymodelc, even if particle location object found
# F:re-use previously calculated particle location object
run_trajec <- F     # whether to generate trajec, see STEP 6
run_foot   <- F     # whether to generate footprint, see STEP 6
if (run_trajec) cat('Need to generate trajec...\n')
if (run_foot)   cat('Need to generate footprint...\n\n')

# if columnTF == F, no need to get ground height or 
#                   perform vertical weighting on trajec-level footprint
columnTF <- T       # column receptor (T) or fixed receptor

# whether to perform XCO2 and its error simulations
run_hor_err   <- T  # T: set error parameters in STEP 4 and set run_foot == T
run_ver_err   <- F  # T: set error parameters in STEP 4
run_emiss_err <- F  # T: get XCO2 error due to prior emiss err, see STEP 7
run_sim       <- T  # T: do analysis with existing trajec/foot, see STEP 7

delt <- 2           # fixed timestep [min]; set = 0 for dynamic timestep
nhrs <- -72         # number of hours backward (-) or forward (+)

# path to grab or store trajec, foot and potential trans err stat DW, 07/31/2018
# ourput directory for storing traj with default convention;
# store traj with wind err in a separate directory if run_hor_err = T
outdir <- file.path(store.path, paste('out', timestr, site, sep = '_'))
if (run_hor_err) outdir <- file.path(err.path, 
                                     paste('out_err', timestr, site, sep = '_'))

# change to Ben's definitions,  see validate_varsiwant()
varstrajec <- c('time', 'indx', 'lati', 'long', 'zagl', 'zsfc', 'foot', 'samt',
                'dmas', 'mlht', 'temp', 'pres', 'sigw', 'tlgr', 'dens')
cat('Done with setting flags...\n')


#------------------------------ STEP 3 --------------------------------------- #
# select receptors --
### 1) Set model receptors, AGLs and particel numbers ***
# for backward fixed-level runs OR forward fixed-level runs
# agl can be a vector, meaning releasing particles from several fixed level
# but if trying to represent air column, use columnTF=T, see below

# if release particles from fixed levels
agl    <- c(10, 500, 2000, 3000)[4]        # in mAGL
numpar <- 100       # par for each time window for forward box runs
dpar   <- numpar

# SET COLUMN RECEPTORS as a list, if release particles from a column
if (columnTF) {

  # min, middle, max heights for STILT levels, in METERS
  minagl <- 0
  midagl <- 3000         # cut-off level for different vertical spacing
  maxagl <- 6000
  dh     <- c(100, 500)  # vertical spacing below and above 'midagl', in METERS
  dpar   <- c(10, 50, 100, 200, 2500)[3]   # particle numbers per level
  agl    <- c(seq(minagl, midagl, dh[1]), seq(midagl + dh[2], maxagl, dh[2]))
  nlev   <- length(agl)
  numpar <- nlev * dpar   # total number of particles
}  # end if columnTF


### 2) place denser receptors within lat range with high XCO2
# whether to select receptors; or simulate all soundings
selTF <- T  
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
recp.num <- NULL     # can be a number for max num of receptors or a vector, 1:20
find.lat <- NULL     # for debug or test, model one sounding


### 3) select satellite soundings, plotTF for whether plotting OCO-2 observed XCO2
recp.info <- get.recp.info(timestr, oco2.ver, oco2.path, lon.lat, selTF, 
                           recp.indx, recp.num, find.lat, agl, plotTF = F, 
                           trajpath = file.path(outdir, 'by-id'), run_trajec)
nrecp <- nrow(recp.info)
cat(paste('Done with receptor setup...total', nrecp, 'receptors..\n'))


#------------------------------ STEP 4 --------------------------------------- #
### 1) path for the ARL format of WRF and GDAS
# simulation_step() will find corresponding met files
met        <- 'gdas0p5'                             # choose met fields
met.path   <- file.path(homedir, 'u0947337', met)  # path of met fields
met.format <- paste0('%Y%m%d_', met)           # met file name convention
met.num    <- 1                                # min number of met files needed


### 2) get horizontal transport error component if run_hor_err = T
# path for outputting wind error stats
hor.err  <- get.uverr(run_hor_err, site, timestr, workdir, overwrite = F,
                      raob.path, raob.format = 'fsl', nhrs, met, met.path, 
                      met.format, lon.lat, agl, err.path)

### if calculating XCO2 error due to horizontal wind error, need biospheric input
# add CT paths and files, DW, 07/28/2018
if (run_hor_err) {
  ct.ver  <- ifelse(substr(timestr, 1, 4) >= '2016', 'v2017', 'v2016-1')
  ct.path <- file.path(input.path, 'CT-NRT', ct.ver)
  ctflux.path <- file.path(ct.path, 'fluxes/optimized')
  ctmole.path <- file.path(ct.path, 'molefractions/co2_total')
} else { ct.ver <- NA; ctflux.path <- NA; ctmole.path <- NA } # end if run_hor_err


### 3) get vertical transport error component if run_ver_err = T
# set zisf = 1 if run_ver_err = F
zisf <- c(0.6, 0.8, 1.0, 1.2, 1.4)[3]; if (!run_ver_err) zisf <- 1.0
pbl.err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = zisf)


### 4) if calculating XCO2 error due to emiss error, need EDGAR and FFDAS files
# rearrange by DW, 10/21/2018
if (run_emiss_err) { 
  edgar.file <- file.path(homedir,
    'lin-group2/group_data/EDGAR/CO2_2000-2008/v42_CO2_2008_TOT.0.1x0.1.nc')
  ffdas.path <- file.path(input.path, 'anthro_inventories/FFDAS/')
  ffdas.file <- list.files(path = ffdas.path, pattern = 'totals')
  ffdas.file <- file.path(ffdas.path, ffdas.file[grep('2008', ffdas.file)])
} else { edgar.file = NA; ffdas.file = NA }  # end if run_emiss_err

cat('Done with choosing met & inputting parameters for error estimates...\n')


#------------------------------ STEP 5 --------------------------------------- #
#### Settings for generating footprint maps
### 1) SET spatial domains and resolution for calculating footprints
foot.res <- 1/120  # footprint resolution, 1km for ODIAC
if (run_emiss_err) foot.res <- 1/10  # for emiss error, generate 0.1deg foot
#foot.res <- 1     # for generating foot that matches CarbonTracker

# these variables will determine resoluation and spatial domain of footprint
# 20x20 degree domain around the city center
# one can also customize the data.frame of `foot.info`
foot.info <- data.frame(xmn = round(lon.lat$citylon) - 10, 
                        xmx = round(lon.lat$citylon) + 10,
                        ymn = round(lon.lat$citylat) - 10, 
                        ymx = round(lon.lat$citylat) + 10,
                        xres = foot.res, yres = foot.res); print(foot.info)
foot.extent <- extent(foot.info$xmn, foot.info$xmx, foot.info$ymn, foot.info$ymx)


### 2) whether weighted footprint by AK and PW for column simulations (X-STILT)
# NA: no weighting performed for fixed receptor simulations
ak.wgt  <- ifelse(columnTF, TRUE, NA)
pwf.wgt <- ifelse(columnTF, TRUE, NA)

# whether to overwrite existing wgttraj file
# if false, read from existing wgttraj.rds
overwrite_wgttraj <- TRUE  


### 3) other footprint parameters using STILTv2 (Fasoli et al., 2018)
hnf_plume      <- T  # whether turn on hyper near-field (HNP) for mising hgts
smooth_factor  <- 1  # Gaussian smooth factor, 0 to disable
time_integrate <- T  # whether integrate footprint along time, T, no hourly foot
projection     <- '+proj=longlat'
cat('Done with footprint setup...\n\n')


#------------------------------ STEP 6 --------------------------------------- #
## if running trajec or footprint
if (run_trajec | run_foot) {

  ## use SLURM for parallel simulation settings
  n_nodes  <- 6
  n_cores  <- ceiling(nrecp/n_nodes)

  # time allowed for running hymodelc before forced terminations
  timeout  <- 24 * 60 * 60  # in sec
  job.time <- '24:00:00'    # total job time
  slurm    <- n_nodes > 0
  slurm_options <- list(time = job.time, account = 'lin-kp', partition = 'lin-kp')

  # if run_hor_err = T, require ODIAC or other flux grids to calculate trans error
  # trans error simulation is performed for each core, see before_footprint_xstilt()
  # DW, 01/29/2019
  if (run_hor_err) {
    emiss.file <- tif2nc.odiacv3(site, timestr, vname = odiac.vname, workdir, 
                                 foot.extent, tiff.path, gzTF = F)
  } else { emiss.file = NA }

  # create a namelist including all variables
  # namelist required for generating trajec
  namelist <- list(ak.wgt = ak.wgt, columnTF = columnTF, ct.ver = ct.ver, 
                   ctflux.path = ctflux.path, ctmole.path = ctmole.path, 
                   delt = delt, emiss.file = emiss.file, foot.info = foot.info, 
                   hnf_plume = hnf_plume, hor.err = hor.err, met = met, 
                   met.format = met.format, met.num = met.num, 
                   met.path = met.path, nhrs = nhrs, n_cores = n_cores,
                   n_nodes = n_nodes, numpar = numpar, outdir = outdir, 
                   oco2.path = oco2.path, overwrite_wgttraj = overwrite_wgttraj,
                   pbl.err = pbl.err, projection = projection, pwf.wgt = pwf.wgt, 
                   recp.info = recp.info, run_foot = run_foot, 
                   run_hor_err = run_hor_err, run_trajec = run_trajec, 
                   slurm = slurm, slurm_options = slurm_options, 
                   smooth_factor = smooth_factor, time_integrate = time_integrate, 
                   timeout = timeout, varstrajec = varstrajec, workdir = workdir)        
  cat('Done with creating namelist...\n')

  # call run_stiltv2() to start running trajec and foot
  run.xstilt(namelist)  # more variables are defined in run.xstilt()
} # end if run trajec or foot


#------------------------------ STEP 7 --------------------------------------- #
### calculate XCO2 concentration and its error (need trajec and footprint ready)
if (!run_trajec & !run_foot & run_sim) {
  library(zoo)

  #------------------------  Horizontal trans error -------------------------- #
  ### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
  # requires two sets of trajectories before running the following section:
  if (run_hor_err) { # this does not need footprint

    ## call function cal.trans.err() to estimate receptor-level trans err [ppm]
    # get actual ppm error, need to have error statistics ready
    # see cal.trajfoot.stat() in called before_footprint_xstilt.r for err stat
    cat('Start simulations of XCO2 error due to horizontal trans err...\n')
    receptor <- cal.trans.err(site, timestr, workdir, outdir, 
                              store.path = err.path, met)
    if (is.null(receptor)) stop('No results calculated, check cal.trans.err()\n')

  } else {
  
    #--------------------- XCO2 or emiss error simulations ---------------------- #
    # requires trajec and footprints ready for running things below, DW, 06/04/2018  
    # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
    cat('Start simulations of XCO2.ff or its error due to emiss err...\n')
    receptor <- run.xco2ff.sim(site, timestr, vname = odiac.vname, tiff.path, 
                               outdir, foot.res, workdir, store.path, nhrs, dpar,
                               smooth_factor, zisf, oco2.ver, met, lon.lat, 
                               run_emiss_err, edgar.file, ffdas.file, 
                               plotTF = F, writeTF = T)
    if (is.null(receptor)) stop('No results calculated, check run.xco2ff.sim()\n')
  } # end if run_hor_err
} # end if run_sim

##### end of script
