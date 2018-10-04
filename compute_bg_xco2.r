#' Compute background XCO2, given different methods:
#' M1. Trajec-endpoint (using CarbonTracker)
#' M2H. Regional daily median (based on Hakkareinen et al., 2016)
#' M2S. Localized normal statistics (based on Silva and Arellano, 2017)
#' M3. X-STILT overpass-specific (based on Wu et al., GMDD)
#' originated from 'create_namelist_oco2-xsilt.r'
#' @author Dien Wu, 04/18/2018

#' @updates:
#' add customized data filtering, DW, 08/20/2018
#' get overpass dates from existing back traj, DW, 08/23/2018
#' add background uncertainty (including spread sd + retrieval err), DW, 09/07/2018
#' add CT-derived background M1, DW, 09/14/2018


## source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
workdir <- file.path(homedir, 'lin-group5/wde/github/XSTILT') #current dir
setwd(workdir)   # move to working directory
source('r/dependencies.r') # source all functions

#------------------------------ STEP 1 --------------------------------------- #
# insert target city
site     <- 'Riyadh'

# please get a google API and insert in the "" as below
register_google(key = '')
lon.lat  <- get.lon.lat(site = site, dlon = 1, dlat = 2)

oco2.ver  <- c('b7rb', 'b8r')[2]  # OCO-2 version
oco2.path <- file.path(homedir, 'lin-group5/wde/input_data/OCO-2/L2', 
                       paste0('OCO2_lite_', oco2.ver))

# path for storing overpass info by 'get,site.track'
txtpath  <- file.path(homedir, 'lin-group5/wde/github/result/oco2_overpass')

# lon.lat: minlon, maxlon, minlat, maxlat, city.lon, city.lat
oco2.track <- get.site.track(site, oco2.ver, oco2.path, searchTF = F,
                             date.range = c('20140901', '20171231'), 
                             thred.count.per.deg = 200, lon.lat = lon.lat, 
                             urbanTF = T, dlon.urban = 0.5, dlat.urban = 0.5,
                             thred.count.per.deg.urban = 100, txtpath = txtpath)
oco2.track  <- oco2.track %>% filter(qf.urban.count > 80)
all.timestr <- oco2.track$timestr
print(all.timestr)

## choose which background method:
# M1. Trajec-endpoint (using CarbonTracker)
# M2H. Regional daily median (based on Hakkareinen et al., 2016)
# M2S. Localized normal statistics (based on Silva and Arellano, 2017)
# M3. X-STILT overpass-specific (based on Wu et al., GMDD)
method <- c('M1', 'M2H', 'M2S', 'M3')[4]

## output and input paths, txtfile name for storing background values
input.path  <- '/uufs/chpc.utah.edu/common/home/lin-group5/wde/input_data'
output.path <- file.path(homedir, 'lin-group5/wde/github/stilt', site)
txtfile     <- file.path(output.path, 
                         paste0(method, '_bg_', site, '_', oco2.ver, '.txt'))

# ----------------------------- M1 Trajec endpoint  ------------------------- #
# need to convolve footprint with diff fluxes and get endpoint CO2 as well
if (method == 'M1') {

  # get CT fluxes and mole fraction paths
  ct.ver    <- ifelse(substr(all.timestr, 1, 4) >= '2016', 'v2017', 'v2016-1')
  flux.path <- file.path(input.path, 'CT-NRT', ct.ver, 'fluxes/optimized')
  mf.path   <- file.path(input.path, 'CT-NRT', ct.ver, 'molefractions/co2_total')
  foot.path <- file.path(output.path, paste0('out_', all.timestr,'_v7_zisf1p0'), 
                         'by-id')
  traj.path <- foot.path 
  
  # call function to get background
  source('r/dependencies.r') # source all functions
  bg.info <- calc.M1.bg(all.timestr, foot.path, traj.path, ct.ver, flux.path, 
                        mf.path, output.path = file.path(output.path, 'bg'), 
                        oco2.ver, txtfile, writeTF = T, nhrs = -72)
} # end if M1

# ------------------------ M2H. Regional daily median  ---------------------- #
if (method == 'M2H')
  bg.info <- calc.M2H.bg(lon.lat, all.timestr, output.path, oco2.ver, oco2.path, 
                         txtfile, plotTF = F)

# ------------------------ M2S. Normal statistics  -------------------------- #
if (method == 'M2S') {
  library(MASS)

  # 4x4 deg box around hotsplot, from Silva and Arellano, 2017
  lon.lat <- get.lon.lat(site, dlat = 2, dlon = 2)

  silva.bg <- NULL
  for (t in 1:length(all.timestr)) {

    # grab observations and calculate the mean and SD
    obs    <- grab.oco2(oco2.path, all.timestr[t], lon.lat) %>% filter(qf == 0)
    tmp.bg <- as.numeric(fitdistr(obs$xco2, 'normal')$estimate[1]) -
              as.numeric(fitdistr(obs$xco2, 'normal')$estimate[2])
    silva.bg <- c(silva.bg, tmp.bg)

  }  # end for t
  
  silva.bg <- data.frame(timestr = all.timestr, silva.bg = silva.bg)
  write.table(silva.bg, quote = F, row.names = F, sep =',', file = txtfile)
} # end if M2S


# -------------------------- M3. Overpass-specific --------------------------- #
# need to run forward runs from a box around the city
if (method == 'M3') {

  #------------------------------ STEP 1 --------------------------------- #
  #### Whether forward/backward, release from a column or a box
  run_trajec <- F  # whether to run forward traj, if T, will overwrite existing
  plotTF     <- T  # whether to calculate background and plot 2D density
  delt       <- 2  # fixed timestep [min]; set = 0 for dynamic timestep
  
  ### MUST-HAVE parameters about errors
  run_hor_err <- T  # run trajec with hor wind errors/calc trans error
  run_ver_err <- F  # run trajec with mixed layer height scaling

  # release particles from a box (dxyp x dxyp) around the city center
  dxyp  <- 0.4               # degree around city center
  nhrs  <- 12                # forward run, nhrs should always be positive

  # release FROM # hrs (e.g., -10) before overpass time, 
  # TO overpass time (e.g., 0), with every 0.5 hour release
  dtime <- seq(-10, 0, 0.5)  # time windows to release particles continuously

  ### release particles from fixed levels
  agl    <- 10         # in mAGL
  numpar <- 1000       # particle # per each time window
  cat(paste('\n\nDone with receptor setup...\n'))

  #------------------------------ STEP 2 --------------------------------- #
  # path for the ARL format of WRF and GDAS
  # simulation_step() will find corresponding met files
  met.indx   <- 4
  met        <- c('hrrr', '1km', 'gdas0p5', 'edas40')[met.indx]
  met.path   <- file.path(homedir, 'u0947337', met)
  met.num    <- 1     # min number of met files needed

  # met file name convention
  met.format <- c('%Y%m%d.%Hz.hrrra', 'wrfout_', '%Y%m%d_gdas0p5', 'edas.%Y%m')[met.indx]
  
  # plot will be stored in 'outpath' as well
  if (run_trajec == T) outpath <- NULL  # will be generated in copies
  if (run_trajec == F) outpath <- file.path(output.path, 'out_forward/') 
  data.filter <- c('QF', 0)      # data filtering on observations

  # **** NEED CHANGES: which side for background, north, south, or both
  # can be interpolated from forward-plume, after being generated
  clean.side <- rep('south', length(all.timestr))

  # path for radiosonde data, 
  # needed if adding wind error component, run_hor_err = T
  raob.path  <- file.path(input.path, 'RAOB/middle.east/riyadh') 

  ### loop over each overpass
  bg.info <- NULL
  for (t in 1 : length(all.timestr)) {

    timestr <- all.timestr[t]
    cat(paste('\n\n## ----- working on overpass on', timestr, '----- ##\n'))

   #------------------------------ STEP 3 --------------------------------- #
    ## get horizontal transport error component if run_hor_err = T
    # calculate from model-data wind comparisons, if available
    hor.err <- get.uverr(run_hor_err, site, timestr, workdir, overwrite = F,
                         raob.path, raob.format = 'fsl', nhrs, met, met.path, 
                         met.format, met.files, lon.lat, agl = c(0, 100), 
                         nfTF = T, forwardTF = T, 
                         err.path = file.path(output.path, 'wind_err'))

    ## get vertical transport error component if run_ver_err = T
    # set zisf = 1 if run_ver_err = F
    zisf    <- c(0.6, 0.8, 1.0, 1.2, 1.4)[3]   
    pbl.err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = zisf)
    cat('Done with choosing met & inputting wind errors...\n')

    #------------------------------ STEP 4 -------------------------------- #
    # !!! need to add makefile for AER_NOAA_branch in fortran ;
    # link two hymodelcs in exe directory
    tmp.info <- run.forward.trajec(site = site, timestr = timestr,
                                   overwrite = run_trajec, nummodel = t,
                                   lon.lat = lon.lat, delt = delt, dxyp = dxyp,
                                   dzp = 0, dtime = dtime, agl = agl,
                                   numpar = numpar, nhrs = nhrs,
                                   workdir = workdir, outpath = outpath,
                                   hor.err = hor.err, pbl.err = pbl.err,
                                   met.format = met.format, met.path = met.path,
                                   met.num = 1, plotTF = plotTF,
                                   oco2.path = oco2.path, oco2.ver = oco2.ver,
                                   zoom = 7, td = 0.05, bg.dlat = 1,
                                   perc = 0.1, clean.side = clean.side[t],
                                   data.filter = data.filter)
    if (is.null(tmp.info)) next
    bg.info <- rbind(bg.info, tmp.info)
  } # end for t

  write.table(bg.info, quote = F, row.names = F, sep = ',', file = txtfile)
} # end if method == 'M3'
