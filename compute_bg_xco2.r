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
#' update code according to v9 changes, DW, 10/19/2018
#' remove data filtering, always use QF = 0 for background estimates, DW, 10/31/2018
#' refaectoring based on new stilt-submodule, DW, 01/25/2019 

## source all functions and load all libraries
# CHANGE working directory ***
homedir <- '/uufs/chpc.utah.edu/common/home'
xstilt_wd <- file.path(homedir, 'lin-group7/wde/X-STILT') #current dir
setwd(xstilt_wd)   # move to working directory
source('r/dependencies.r') # source all functions

# insert your google API for the use of ggplot and ggmap
api.key <- ''
register_google(key = api.key)


#------------------------------ STEP 1 --------------------------------------- #
site     <- 'Riyadh'  # insert target city
lon.lat  <- get.lon.lat(site = site, dlon = 1, dlat = 2)
site.lon <- lon.lat$citylon; site.lat <- lon.lat$citylat
oco2.ver <- c('b7rb', 'b8r', 'b9r')[1]  # OCO-2 version

## choose which background method:
# M1. Trajec-endpoint (using CarbonTracker)
# M2H. Regional daily median (based on Hakkareinen et al., 2016)
# M2S. Localized normal statistics (based on Silva and Arellano, 2017)
# M3. X-STILT overpass-specific (based on Wu et al., GMDD)
method <- c('M1', 'M2H', 'M2S', 'M3')[4]

## required output and input paths, txtfile name for storing background values
input.path  <- file.path(homedir, 'lin-group7/wde/input_data')
output.path <- file.path(homedir, 'lin-group7/wde/output', site)
bg.txtfile  <- file.path(output.path, paste0(method, '_bg_', site, '_', 
                                             oco2.ver, '.txt'))
oco2.path   <- file.path(input.path, 'OCO-2/L2', paste0('OCO2_lite_', oco2.ver))

# path for storing overpass sampling info by 'get,site.track'
txt.path   <- file.path(input.path, 'OCO-2/overpass_city') 
oco2.track <- get.site.track(site, oco2.ver, oco2.path, searchTF = F,
                             date.range = c('20140901', '20181231'), 
                             thred.count.per.deg = 100, lon.lat = lon.lat, 
                             urbanTF = T, dlon.urban = 0.5, dlat.urban = 0.5,
                             thred.count.per.deg.urban = 100, txt.path) %>% 
              filter(qf.urban.count > 80)
all.timestr <- oco2.track$timestr; print(all.timestr)


# ----------------------------- M1 Trajec endpoint  ------------------------- #
# need to convolve footprint with diff fluxes and get endpoint CO2 as well
if (method == 'M1') {

  # get CT fluxes and mole fraction paths
  ct.ver    <- ifelse(substr(all.timestr, 1, 4) >= '2016', 'v2017', 'v2016-1')
  flux.path <- file.path(input.path, 'CT-NRT', ct.ver, 'fluxes/optimized')
  mf.path   <- file.path(input.path, 'CT-NRT', ct.ver, 'molefractions/co2_total')

  # trajec and footprint paths, pointing to 'by-id', change them if needed
  foot.path <- file.path(output.path, paste('out', timestr, site, sep = '_'), 
                         'by-id')
  traj.path <- foot.path 
  
  # call function to get background
  bg.info <- calc.bg.M1(all.timestr, foot.path, traj.path, ct.ver, flux.path, 
                        mf.path, output.path, oco2.ver, txtfile = bg.txtfile, 
                        writeTF = T, nhrs = -72)
} # end if M1

# ------------------------ M2H. Regional daily median  ---------------------- #
if (method == 'M2H') bg.info <- calc.bg.M2H(lon.lat, all.timestr, output.path, 
                                            oco2.ver, oco2.path, 
                                            txtfile = bg.txtfile, plotTF = F)

# ------------------------ M2S. Normal statistics  -------------------------- #
if (method == 'M2S') bg.info <- calc.bg.M2S(site, all.timestr, oco2.path, 
                                            oco2.ver, txtfile = bg.txtfile)

# -------------------------- M3. Overpass-specific --------------------------- #
# need to run forward runs from a box around the city
if (method == 'M3') {

  #------------------------------ STEP 1 --------------------------------- #
  #### Whether forward/backward, release from a column or a box
  run_trajec  <- T  # whether to run forward traj, if T, will overwrite existing
  run_bg      <- F  # whether to calculate background and plot 2D density 
                    # this requires forward trajec ready
  run_hor_err <- T  # run trajec with hor wind errors/calc trans error
  run_ver_err <- F  # run trajec with mixed layer height scaling
  traj.path <- file.path(output.path, 'out_forward')

  # path for radiosonde data, FSL format from NOAA
  # needed if adding wind error component, run_hor_err = T
  raob.path <- file.path(input.path, 'RAOB', site)  # radiosonde
  err.path  <- file.path(output.path, 'wind_err') # needed if run_hor_err = T

  # if one does not want to run the whole wind error interpolation 
  # (it needs RAOB data and takes time), thus set overwrite == F 
  # and give a reasonable wind RMSE value, m/s
  overwrite <- F
  siguverr  <- 3     # wind RMSE, m/s

  # set zisf = 1 if run_ver_err = F
  zisf <- c(0.6, 0.8, 1.0, 1.2, 1.4)[3]; if (!run_ver_err) zisf <- 1.0       
  
  # release particles from a box (dxyp x dxyp) around the city center
  dxyp  <- 0.4               # degree around city center
  nhrs  <- 12                # forward run, nhrs should always be positive

  # release FROM # hrs (e.g., -10) before overpass time, 
  # TO overpass time (e.g., 0), with every 0.5 hour release
  dtime <- seq(-10, 0, 0.5)  # time windows to release particles continuously
  
  ### release particles from fixed levels, recort particles for every 2mins
  delt   <- 2          # in mins
  agl    <- 10         # in mAGL
  numpar <- 1000       # particle # per each time window
  cat(paste('\n\nDone with receptor setup...\n'))

  #------------------------------ STEP 2 --------------------------------- #
  # path for the ARL format of WRF and GDAS
  # simulation_step() will find corresponding met files
  met.indx   <- 2
  met        <- c('hrrr', 'gdas0p5', 'edas40', 'wrf')[met.indx]
  met.path   <- file.path(homedir, 'u0947337', met)
  met.num    <- 1     # min number of met files needed
  met.format <- c('%Y%m%d.%Hz.hrrra', '%Y%m%d_gdas0p5', 'edas.%Y%m', '%Y%m%d')[met.indx]
  
  #------------------------------ STEP 4 -------------------------------- #
  # !!! need to add makefile for AER_NOAA_branch in fortran ;
  # link two hymodelcs in exe directory
  # parallel running, DW, 11/06/2019
  n_nodes <- 1
  n_cores <- ceiling(length(all.timestr) / n_nodes)
  slurm_options <- list(time = '04:00:00', account = 'lin-kp', partition = 'lin-kp')
  
  if (run_trajec) {
    dtime <- list(dtime)
    xstilt_apply(FUN = run.forward.trajec, slurm = T, slurm_options, n_nodes, 
                 n_cores, jobname = paste0('XSTILT_forward_', site), 
                 site, site.lon, site.lat, timestr = all.timestr, 
                 run_trajec, run_hor_err, run_ver_err, agl, delt, dtime, 
                 dxyp, dzp = 0, numpar, nhrs, xstilt_wd, output.path = traj.path, 
                 met, met.format, met.path, met.num, raob.path, 
                 raob.format = 'fsl', siguverr, err.path, nhrs.zisf = 24, zisf)
    q('no')
  }


  # ----- calculating the forward background and plot kernel density map ----- #
  # **** NEED CHANGES: which side for background: 'north', 'south', or 'both'
  # can be interpolated from forward-plume, after being generated
  clean.side <- rep('both', length(all.timestr))  # just an example

  if (run_trajec == F & run_bg) { # need forward trajec ready
    
    bg.info <- NULL
    for (tt in 1: length(all.timestr)) {  ### loop over each overpass
      tmp.bg <- calc.bg.forward.trajec(trajpath = traj.path, site, 
                                       timestr = all.timestr[tt], agl, dtime, 
                                       numpar, dxyp, oco2.path, oco2.ver, qfTF = T, 
                                       met, zoom = 8, td = 0.05, bg.dlat = 0.5, 
                                       perc = 0.2, clean.side[tt], api.key)
      bg.info <- rbind(bg.info, tmp.bg)
    }
   
    write.table(bg.info, quote = F, row.names = F, sep = ',', file = bg.txtfile)
  } # end if run_bg
    
} # end if method == 'M3'
