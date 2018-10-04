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
# Add horizontal transport error module with flag 'run_hor_err'
#   with additional subroutines in /src/oco2-xstilt/trans_err, DW, 07/23/2018
#
# Add vertical trans error module, with flag 'run_ver_err',
# ziscale when scaling PBL heights and break trans error part (step 3) into
# horizontal wind error and vertical trans error (via PBL perturb), DW, 07/25/2018
#
# simplify main code and use 'get.uverr()' to get horizontal wind error stat,
# DW, 08/31/2018
#
# add emission uncertainty, DW, 09/06/2018
# -----------------------------------------------------------------------------

#### source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT') # current dir
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions


#------------------------------ STEP 1 --------------------------------------- #
# input dlat, dlon to get spatial domain around city center
site <- 'Riyadh'   # choose a city

# please get a google API and insert in the "" as below
register_google(key = '')
lon.lat  <- get.lon.lat(site = site, dlon = 1, dlat = 2)

# required paths
oco2.ver   <- c('b7rb', 'b8r')[1]  # OCO-2 version
input.path <- file.path(homedir, 'lin-group5/wde/input_data')
oco2.path  <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_', oco2.ver))
sif.path   <- file.path(input.path, paste0('OCO-2/L2/OCO2_lite_SIF_', oco2.ver))
raob.path  <- file.path(input.path, 'RAOB/middle.east/riyadh')  # path for RAOB

output.path <- file.path(homedir, 'lin-group5/wde/github/result')
txt.path    <- file.path(output.path, 'oco2_overpass')
trans.path  <- file.path(output.path, 'trans_err')  # path for storing CO2 statistics

# date range for searching OCO-2 tracks, min, max YYYYMMDD
date.range  <- c('20140101', '20181231')

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
oco2.track <- get.site.track(site, oco2.ver, oco2.path, searchTF,
                             date.range, thred.count.per.deg, lon.lat, 
                             urbanTF, dlon, dlat, thred.count.per.deg.urban,
                             txt.path)

# one can further subset 'oco2.track' based on sounding # or data quality
# over entire lon.lat domain or near city center
# see columns 'qf.count' or 'wl.count' in 'oco2.track'
oco2.track <- oco2.track %>% filter(qf.urban.count > 50)

# remove summtertime tracks
oco2.track <- oco2.track %>%
  filter(substr(timestr, 5, 6) < '05' | substr(timestr, 5, 6) > '08')

# finally narrow down and get timestr
all.timestr <- oco2.track$timestr[c(2, 3, 5, 6, 7)]         # v7, Riyadh
print(all.timestr)

# once you have all timestr, you can choose whether to plot them on maps
# this helps you choose which overpass to simulate first, see 'tt' below
plotTF <- F
if (plotTF) {
  plotdir <- file.path(homedir, 'lin-group5/wde/github/stilt', site)
  dir.create(plotdir, showWarnings = F)

  for(t in 1:length(all.timestr)){
    ggmap.obs.xco2(site, all.timestr[t], oco2.path, lon.lat, workdir, plotdir)
    ggmap.obs.sif(site, all.timestr[t], sif.path, lon.lat, workdir, plotdir)
  }
}  # end if plotTF


# *** NOW choose the timestr that you would like to work on...
tt <- 2
timestr <- all.timestr[tt]
cat(paste('Working on:', timestr, 'for city/region:', site, '...\n\n'))
cat('Done with choosing cities & overpasses...\n')


#------------------------------ STEP 2 --------------------------------------- #
#### Whether forward/backward, release from a column or a box
columnTF   <- T      # whether a column receptor or fixed receptor
forwardTF  <- F      # forward or backward traj, if forward, release from a box

# T:rerun hymodelc, even if particle location object found
# F:re-use previously calculated particle location object
run_trajec <- T      # whether to generate trajec
run_foot   <- T      # whether to generate footprint
run_sim    <- T      # whether to calculate simulated XCO2.ff, see STEP 8
if (run_trajec) cat('Need to generate trajec...\n')
if (run_foot)   cat('Need to generate footprint...\n\n')

# whether to generate trajec with horizontal wind err component/to calc trans error
# OR generate trajec with PBL scaling
run_hor_err   <- F    # T: set parameters in STEP 3, call functions in STEP 9
run_ver_err   <- F    # T: set parameters in STEP 3, call functions in STEP 9
run_emiss_err <- F    # T: for getting XCO2 error due to prior emiss uncertainty

stilt.ver     <- 2    # STILT versions (call different footprint algorithms)
delt          <- 2    # fixed timestep [min]; set = 0 for dynamic timestep
nhrs          <- -72  # number of hours backward (-) or forward (+)

# where to grab or store trajec, DW, 07/31/2018
if (run_trajec) {
  # ourput directory for storing traj, default convention
  outdir <- file.path(workdir, paste('out', timestr, site, sep = '_'))
} else {
  # ourput directory that stores existing trajec
  outdir <- file.path(homedir, 'lin-group5/wde/github/stilt', 
                      site, paste0('out_', timestr, '_v7_zisf1p0'))  
}  # end if run_trajec


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
  peak.lat <- c(lon.lat$citylat - 0.5, lon.lat$citylat + 0.5)

  # number of points to aggregate within 1deg over small/large enhancements,
  # i.e., over background/enhancements, binwidth will be 1deg/num
  num.bg   <- 20   # e.g., every 20 pts in 1 deg
  num.peak <- 40   # e.g., every 40 pts in 1 deg

  # recp.indx: how to pick receptors from all screened soundings (QF = 0)
  recp.indx <- c(seq(lon.lat$minlat,  peak.lat[1], 1/num.bg),
                 seq(peak.lat[1], peak.lat[2], 1/num.peak),
                 seq(peak.lat[1], lon.lat$maxlat,  1/num.bg))

} else { recp.indx <- NULL }

# whether to subset receptors when debugging
recp.num <- NULL           # can be a number for max num of receptors
find.lat <- NULL           # for debug or test, model one sounding
data.filter <- c('QF', 0)  # use WL or QF to filter data

# select satellite soundings, plotTF for whether plotting OCO-2 observed XCO2
recp.info <- get.recp.info(timestr, oco2.path, lon.lat, selTF, recp.indx,
                           recp.num, find.lat, agl, plotTF = F, 
                           trajpath = file.path(outdir, 'by-id'), 
                           stilt.ver, data.filter)

nrecp <- nrow(recp.info)
cat(paste('Done with receptor setup...total', nrecp, 'receptors..\n'))


#------------------------------ STEP 4 --------------------------------------- #
## path for the ARL format of WRF and GDAS
# simulation_step() will find corresponding met files
met.indx   <- 3
met        <- c('hrrr', 'wrf', 'gdas0p5')[met.indx] # choose met fields
met.path   <- file.path(homedir, 'u0947337', met)

# met file name convention
met.format <- c('%Y%m%d.%Hz.hrrra', 'wrfout_', '%Y%m%d_gdas0p5')[met.indx]
met.num    <- 1     # min number of met files needed
met.files  <- NULL
if (met == 'wrf')
met.files  <- list.files(path = met.path, pattern = substr(timestr, 1, 6),
                         full.names = T)

## get horizontal transport error component if run_hor_err = T
hor.err <- get.uverr(run_hor_err, site, timestr, workdir, overwrite = F,
                     raob.path, raob.format = 'fsl', nhrs, met, met.path, 
                     met.format, met.files, lon.lat, agl)

## get vertical transport error component if run_ver_err = T
# set zisf = 1 if run_ver_err = F
zisf <- c(0.6, 0.8, 1.0, 1.2, 1.4)[3]
pbl.err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = zisf)
cat('Done with choosing met & inputting wind errors...\n')


#------------------------------ STEP 5 --------------------------------------- #
#### Settings for generating footprint maps
## SET spatial domains and resolution for calculating footprints
foot.res <- 1/120  # footprint resolution, 1km for ODIAC
#foot.res <- 1/10  # for emiss error, generate 0.1deg foot
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
print(foot.info)

## whether weighted footprint by AK and PW for column simulations
if (columnTF) {
  ak.wgt <- T; pwf.wgt  <- T
} else {  # no weighting needed for fixed receptor simulations
  ak.wgt <- NA; pwf.wgt <- NA
}

dmassTF <- F        # 1st-order correction on dmass for footprint using STILTv1

# other footprint parameters using STILTv2
hnf_plume      <- T  # whether turn on hyper near-field (HNP) for mising hgts
smooth_factor  <- 1  # Gaussian smooth factor, 0 to disable
time_integrate <- T  # whether integrate footprint along time, T, no hourly foot
projection     <- '+proj=longlat'
cat('Done with footprint setup...\n')


#------------------------------ STEP 6 --------------------------------------- #
#### !!! NO NEED TO CHANGE ANYTHING LISTED BELOW -->
# create a namelist including all variables
# namelist required for generating trajec
namelist <- list(agl = agl, ak.wgt = ak.wgt, delt = delt, dmassTF = dmassTF,
                 dpar = dpar, foot.info = foot.info, forwardTF = forwardTF,
                 hnf_plume = hnf_plume, hor.err = hor.err, lon.lat = lon.lat, 
                 met = met, met.format = met.format, met.num = met.num, 
                 met.path = met.path, nhrs = nhrs, numpar = numpar, 
                 outdir = outdir, oco2.path = oco2.path, pbl.err = pbl.err,
                 projection = projection, pwf.wgt = pwf.wgt, 
                 recp.info = recp.info, run_foot = run_foot, 
                 run_trajec = run_trajec, site = site,
                 smooth_factor = smooth_factor, stilt.ver = stilt.ver,
                 time_integrate = time_integrate, timestr = timestr, 
                 varstrajec = varstrajec, workdir = workdir)


#------------------------------ STEP 7 --------------------------------------- #

## if running trajec or footprint ----------------------
if (run_trajec | run_foot) {

  ## use SLURM for parallel simulation settings
  n_nodes  <- 6
  n_cores  <- ceiling(nrecp/n_nodes)
  job.time <- '24:00:00'
  slurm    <- n_nodes > 0
  namelist$slurm_options <- list(time = job.time, account = 'lin-kp',
                                 partition = 'lin-kp')

  # time allowed for running hymodelc before forced terminations
  timeout  <- 24 * 60 * 60  # in sec
  namelist <- c(namelist, n_cores = n_cores, n_nodes = n_nodes, slurm = slurm,
                timeout = timeout)
  cat('Done with creating namelist...\n')

  # call run_stiltv2(), start running trajec and foot
  run.stiltv2(namelist = namelist)

} else if (run_sim | run_emiss_err) {

  ## if running trajec or footprint --------------------
  # requires trajec and footprints ready for running this sections:
  # Simulate XCO2.ff using ODIAC emissions, DW, 06/04/2018
  # use footprint domain to crop emissions
  foot.path   <- file.path(outdir, 'by-id')
  foot.file   <- list.files(foot.path, pattern = 'foot.nc', recursive = T,
                            full.names = T)
  tmp.foot    <- raster(foot.file[1])
  foot.extent <- extent(tmp.foot)
  cat('Done reading footprint..\n')

  # before simulations, subset emissions and convert tif format to nc format
  vname <- '2017'
  tiff.path <- file.path(homedir,
    paste0('lin-group2/group_data/ODIAC/ODIAC', vname))  # tif file from ODIAC website

  # call tif2nc.odiacv2() to subset and get emiss file name
  # 'store.path' is the path for outputting emissions
  cat('Start reading and subsetting emissions that match foot...\n')
  emiss.path <- file.path(workdir, 'in')
  dir.create(emiss.path, showWarnings = F)  # create the dir to store emiss
  tiff.path <- file.path(tiff.path, substr(timestr, 1, 4))

  # get cropped ODIAC emission for the overpass month
  emiss.file <- tif2nc.odiacv3(site, timestr, vname, workdir, foot.extent,
                               store.path = emiss.path, tiff.path, gzTF = F)

  if (run_emiss_err) {  
   
    # calculating absolute emission errors
    tiff.path <- file.path(homedir,
      paste0('lin-group2/group_data/ODIAC/ODIAC', vname), '2008')

    # get the ODAIC emissions for year 2008 for conducting emission error
    odiac.file.2008 <- tif2nc.odiacv3(site, timestr = '20081229', vname, workdir,
                                      foot.extent, store.path = emiss.path, 
                                      tiff.path, gzTF = F)

    edgar.file <- file.path(homedir,
      'lin-group2/group_data/EDGAR/CO2_2000-2008/v42_CO2_2008_TOT.0.1x0.1.nc')

    ffdas.path <- file.path(input.path, 'anthro_inventories/FFDAS/')
    ffdas.file <- list.files(path = ffdas.path, pattern = 'totals')
    ffdas.file <- file.path(ffdas.path, ffdas.file[grep('2008', ffdas.file)])

    # update emiss.file with absolute emission uncertainty
    # that can further convolve with footprints
    source('r/dependencies.r') # source all functions
    emiss.file <- cal.emiss.err(site, timestr, odiac.file.2008, edgar.file,
                                ffdas.file, emiss.file, overwrite = F)

    # txt file name for outputting model results
    txtfile <- file.path(workdir, 
      paste0(timestr, '_', site, '_XCO2ff_emiss_err_', abs(nhrs), 'hrs_', dpar, 
             'dpar_sf', smooth_factor, '_zisf', zisf, '_', oco2.ver, '_', met,
             '.txt'))

  } else {  # or just grab ODIAC for emissions
   
    # txt file name for outputting model results
    txtfile <- file.path(workdir, 
      paste0(timestr, '_', site, '_XCO2ff_', abs(nhrs), 'hrs_', dpar, 'dpar_sf', 
             smooth_factor, '_zisf', zisf, '_', oco2.ver, '_', met, '.txt'))

  }  # end if emission err
  print(txtfile)

  # reformatted ODIAC emissions file name should include 'YYYYMM'
  # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
  cat('Start XCO2.ff simulations...\n')
  receptor <- foot.odiacv3(foot.file, emiss.file, workdir, txtfile, lon.lat,
                           plotTF = F, writeTF = T)

} # end if run traj/foot/sim


#------------------------------ STEP 9 --------------------------------------- #
### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
# requires two sets of trajectories before running the following section:
if (run_hor_err) {

  # for ffco2 emission path and files
  vname <- '2017'

  # tif file from ODIAC website
  tiff.path <- file.path(homedir, 
    paste0('lin-group2/group_data/ODIAC/ODIAC', vname), substr(timestr, 1,4))  

  # grab footprint domain
  foot.extent <- extent(as.numeric(foot.info[1:4]))
  emiss.file <- tif2nc.odiacv3(site, timestr, vname, workdir, foot.extent,
                               store.path = file.path(workdir, 'in', 'ODIAC'), 
                               tiff.path, gzTF = F)

  # load two paths for trajec before and after perturbations
  traj.path1 <- file.path(homedir, 
     'lin-group5/wde/github/cp_trajecfoot/out_2014122910_100dpar/by-id')
  traj.path2 <- file.path(workdir, 'out/by-id')

  # add CT paths and files, DW, 07/28/2018
  ct.ver <- '2016-1'; if (timestr >= 20160101) ct.version <- '2017'
  ct.path <- '/uufs/chpc.utah.edu/common/home/lin-group2/group_data/CT-NRT'
  ctflux.path <- file.path(ct.path, paste0('v', ct.ver), 'fluxes/optimized')
  ctmole.path <- file.path(ct.path, paste0('v', ct.ver), 'molefractions/co2_total')

  # txt file name for outputting model results
  txtfile <- file.path(workdir, 
    paste0(timestr, '_', site, '_trans_err_', met, '.txt'))

  namelist <- c(namelist, emiss.file = emiss.file, traj.path1 = traj.path1,
                traj.path2 = traj.path2, ct.ver = ct.ver, 
                ctflux.path = ctflux.path, ctmole.path = ctmole.path, 
                store.path = trans.path, txtfile = txtfile)

  # call cal.trans.err to estimate trans err
  trans.err <- cal.trans.err(namelist)
} # end if run_hor_err


##### end of script
