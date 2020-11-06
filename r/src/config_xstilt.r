# STILT R Executable
# For documentation, see https://uataq.github.io/stilt/
# Ben Fasoli
# edited by Dien Wu for X-STILT


# preserve original 'run_stilt.r' and modify 'run_stilt.r' as a subroutine
# add ziscale, DW, 07/25/2018
# add two before_* functions for getting ground height and AK PW weighting, DW, 01/24/2019
# allows for generate footprints with different resolutions, DW, 02/11/2019 
# add additional variables due to changes in simulation_step()'s argument, DW, 07/02/2020
# add/remove variables as migrating to STILT-HYSPLIT, DW, 07/03/2020 

# refactoring: move some of the calculation in the main script to this function
# DW, 10/26/2020

config_xstilt <- function(namelist){

  # User inputs ----------------------------------------------------------------
  xstilt_wd <- namelist$xstilt_wd
  lib.loc   <- .libPaths()[1]
  site      <- namelist$site
  timestr   <- namelist$timestr
  lon.lat   <- namelist$lon.lat[[1]]
  oco.ver   <- namelist$oco.ver
  oco.path  <- namelist$oco.path
  oco.sensor <- namelist$oco.sensor
  store.path <- namelist$store.path

  # Model control -------------------------------------------------------------
  rm_dat      <- T
  timeout     <- namelist$timeout  # in sec
  run_trajec  <- namelist$run_trajec
  run_foot    <- namelist$run_foot
  run_sim     <- namelist$run_sim
  run_hor_err <- namelist$run_hor_err
  run_ver_err <- namelist$run_ver_err

  if (run_trajec) cat('Need to generate trajec...\n')
  if (run_foot)   cat('Need to generate footprint...\n')
  if (!run_trajec & !run_foot & !run_sim) 
    cat('NO calculations will be performed, please check run_* flags..\n')

  # line source for agl
  agl     <- c(namelist$minagl, namelist$maxagl) 
  n_hours <- namelist$nhrs
  numpar  <- namelist$numpar

  varsiwant <- namelist$varstrajec
  if (length(varsiwant) == 0)
    varsiwant  <- c('time', 'indx', 'long', 'lati', 'zagl', 'zsfc', 'foot', 
                    'mlht', 'dens', 'samt', 'sigw', 'tlgr', 'temp', 'pres')

  # met fields ----------------------------------------------------------------
  met       <- namelist$met
  met_path  <- namelist$met_path
  n_met_min <- namelist$n_met_min
  met_file_format <- namelist$met_file_format
  met_subgrid_buffer <- namelist$met_subgrid_buffer
  met_subgrid_enable <- namelist$met_subgrid_enable
  met_subgrid_levels <- namelist$met_subgrid_levels


  # Footprint params ----------------------------------------------------------
  xmn  <- round(lon.lat$citylon) - namelist$foot_dlon
  xmx  <- round(lon.lat$citylon) + namelist$foot_dlon
  ymn  <- round(lon.lat$citylat) - namelist$foot_dlat
  ymx  <- round(lon.lat$citylat) + namelist$foot_dlat
  xres <- namelist$foot_res
  yres <- namelist$foot_res
  hnf_plume      <- namelist$hnf_plume
  smooth_factor  <- namelist$smooth_factor
  time_integrate <- namelist$time_integrate
  projection     <- namelist$projection 


  # TROPOMI info if needed ------------------------------------------------
  tropomi.speci <- unlist(namelist$tropomi_speci)
  if (!NA %in% tropomi_speci) {

    if (!is.null(namelist$tropomi_hir_path))
      tropomi.path <- unlist(ifelse(substr(timestr, 1, 8) >= 20190806, 
                                    namelist$tropomi_hir_path, 
                                    namelist$tropomi_path))

    tropomi.config <- config_tropomi(timestr, tropomi.speci, tropomi.path, lon.lat)
    tropomiTF  <- tropomi.config$tropomiTF  # if true, separate OCO from TROPOMI runs
    tropomi.fn <- tropomi.config$tropomi.fn
    if (tropomiTF) cat('Treat OCO and TROPOMI as two separate runs given different overpass hours...\n')
    if (!tropomiTF) cat('Treat OCO and TROPOMI as one single run...\n')

    # if tropomiTF == T, trajec for OCO and TROPOMI will be treated as two separate runs; 
    # however, if run_tropomi == FALSE, still perform the OCO runs and turn tropomiTF off
    if ( !namelist$run_tropomi & tropomiTF ) {
      cat('\nWork on OCO runs since run_tropomi is turned off...\n')
      
      # if we focus on the OCO part of the run, turn tropomiTF off and set 
      # tropomi.species as NA, so that X-STILT will not work on TROPOMI weighting
      # in before_footprint_xstilt.r, DW, 10/28/2020
      tropomiTF = FALSE; tropomi.speci = tropomi.fn = NA

    } else cat('Working on the TROPOMI runs...\n')

  } else { tropomiTF = FALSE; tropomi.fn = NA}


  # Namimg output dir ---------------------------------------------------------
  # path to grab or store trajec, foot and potential trans err stat DW, 07/31/2018
  output_wd <- file.path(store.path, paste('out', timestr, met, oco.sensor, sep = '_'))
  if (tropomiTF) 
    output_wd <- file.path(store.path, paste('out', timestr, met, 'TROPOMI', sep = '_'))
  if (run_hor_err) output_wd <- gsub('out', 'outerr', output_wd)
  print(output_wd)


  # Compute XCO2 using existing traj & foot ------------------------------------
  # calculate XCO2 concentration and its error (need trajec and footprint ready)
  if (!run_trajec & !run_foot & run_sim) {

    # ------------------------  Horizontal trans error ---------------------- #
    ### simulate transport error in XCO2 due to met errors, DW, 07/25/2018
    # requires two sets of trajectories before running the following section:
    if (run_hor_err) { # this does not need footprint

      ## call function cal.trans.err() to estimate receptor-level trans err [ppm]
      # get actual ppm error, need to have error statistics ready
      # see cal.trajfoot.stat() in called before_footprint_xstilt.r for err stat
      cat('Start simulations of XCO2 error due to horizontal trans err...\n')
      result <- cal.trans.err(site, timestr, workdir = xstilt_wd, 
                              outdir = output_wd, store.path, met)
      if (is.null(result)) stop('No results calculated, check cal.trans.err()\n')

    } else {

      # ------------------- XCO2 or emiss error simulations ----------------- #
      # requires trajec and footprints ready for running things below, DW, 06/04/2018  
      # call func to match ODIAC emissions with xfoot & sum up to get 'dxco2.ff'
      cat('Start simulations of XCO2.ff or its error due to emiss err...\n')
      result <- run.xco2ff.sim(site, timestr, vname = namelist$odiac.ver, 
                               tiff.path = namelist$odiac.path, 
                               outdir = output_wd, foot.res = xres, 
                               workdir = xstilt_wd, store.path, nhrs = n_hours, 
                               oco.sensor, oco.ver, met, lon.lat,
                               run_emiss_err = namelist$run_emiss_err, 
                               edgar.file = namelist$edgar.file, 
                               ffdas.file = namelist$ffdas.file)
      if (is.null(result)) stop('No results calculated, check run.xco2ff.sim()\n')
      return(txt.file)
    } # end if run_hor_err
  } # end if run_sim



  # Receptor locations & time -------------------------------------------------
  # place denser receptors within lat range with high XCO2
  # whether to select receptors; or simulate all soundings
  selTF <- namelist$selTF  
  if (selTF) {  # recp.indx: how to subset from all screened soundings (QF = 0)
    peak.lat <- namelist$peak.lat 
    recp.indx <- c(seq(lon.lat$minlat,  peak.lat[1],     1 / namelist$num.bg),
                   seq(peak.lat[1],     peak.lat[2],     1 / namelist$num.peak),
                   seq(peak.lat[1],     lon.lat$maxlat,  1 / namelist$num.bg))
  } else recp.indx <- NULL

  # select satellite soundings, plotTF for whether plotting OCO-2 observed XCO2
  receptors <- get.recp.info(timestr, oco.ver, oco.path, lon.lat, selTF, 
                             recp.indx, recp.num = NULL, find.lat = NULL, 
                             agl, run_trajec, outdir = output_wd, 
                             data.filter = unlist(namelist$data.filter), 
                             tropomiTF, tropomi.fn = tropomi.fn[1])
  cat(paste('Done with receptor setup...total', nrow(receptors), 'receptor(s)..\n'))


  # Error stats if needed -----------------------------------------------------
  # Calculating error stats for horizontal and vertical transport errors
  hor_err <- get.uverr(run_hor_err, site, timestr, workdir = xstilt_wd, 
                       overwrite = F, raob.path = namelist$raob.path, 
                       raob.format = 'fsl', nhrs, met, met.path = met_path, 
                       met.format = met_file_format, lon.lat, agl, 
                       err.path = file.path(store.path, 'wind_err'))
  pbl_err <- get.zierr(run_ver_err, nhrs.zisf = 24, const.zisf = namelist$zisf)

  horcoruverr <- hor_err$horcoruverr
  siguverr    <- hor_err$siguverr
  tluverr     <- hor_err$TLuverr
  zcoruverr   <- hor_err$zcoruverr
  horcorzierr <- pbl_err$horcorzierr
  sigzierr    <- pbl_err$sigzierr
  tlzierr     <- pbl_err$TLzierr
  zicontroltf <- pbl_err$zicontroltf
  ziscale     <- pbl_err$ziscale   # prescribe PBL scaling, vector form

  # and prepare ODIAC based on footprint domain 
  if (run_hor_err) {
    emiss.file <- tif2nc.odiacv3(site, timestr, vname = namelist$odiac.ver, 
                                 workdir = xstilt_wd, 
                                 foot.ext = extent(xmn, xmx, ymn, ymx), 
                                 tiff.path = namelist$odiac.path, gzTF = F)
  } else emiss.file = NA


  # Transport and dispersion settings, use default setting in STILT-R v2
  capemin     <- -1
  cmass       <- 0
  conage      <- 48
  cpack       <- 1
  delt        <- 2
  dxf         <- 1
  dyf         <- 1
  dzf         <- 0.01
  efile       <- ''
  emisshrs    <- 0.01
  frhmax      <- 3
  frhs        <- 1
  frme        <- 0.1
  frmr        <- 0
  frts        <- 0.1
  frvs        <- 0.01
  hscale      <- 10800
  ichem       <- 8
  idsp        <- 2
  initd       <- 0
  k10m        <- 1
  kagl        <- 1
  kbls        <- 1
  kblt        <- 5
  kdef        <- 0
  khinp       <- 0
  khmax       <- 9999
  kmix0       <- 250
  kmixd       <- 3
  kmsl        <- 0
  kpuff       <- 0
  krand       <- 4
  krnd        <- 6
  kspl        <- 1
  kwet        <- 1
  kzmix       <- 0
  maxdim      <- 1
  maxpar      <- numpar
  mgmin       <- 10
  mhrs        <- 9999
  nbptyp      <- 1
  ncycl       <- 0
  ndump       <- 0
  ninit       <- 1
  nstr        <- 0
  nturb       <- 0
  nver        <- 0
  outdt       <- 0
  p10f        <- 1
  pinbc       <- ''
  pinpf       <- ''
  poutf       <- ''
  qcycle      <- 0
  rhb         <- 80
  rht         <- 60
  splitf      <- 1
  tkerd       <- 0.18
  tkern       <- 0.18
  tlfrac      <- 0.1
  tout        <- 0
  tratio      <- 0.75
  tvmix       <- 1
  veght       <- 0.5
  vscale      <- 200
  vscaleu     <- 200
  vscales     <- -1
  wbbh        <- 0
  wbwf        <- 0
  wbwr        <- 0
  wvert       <- FALSE
  w_option    <- 0
  zicontroltf <- 0
  ziscale     <- rep(list(rep(1, 24)), nrow(receptors))
  z_top       <- 25000

  # customized functions
  # before_*_xstilt() are two customized functions for OCO-2/XSTILT
  before_footprint <- before_footprint_xstilt
  before_trajec    <- before_trajec_xstilt

  # Parallel simulation settings
  slurm   <- namelist$slurm
  n_cores <- namelist$n_cores
  n_nodes <- namelist$n_nodes

  # OR specify the ammount of memory per cpu your job requires
  #mem_per_cpu <- 10 * 1024    
  slurm_options <- list(time = namelist$job.time, 
                        account = namelist$slurm_account, 
                        partition = namelist$slurm_partition, 
                        mem = namelist$mem_per_node)
  jobname <- paste0('XSTILT_', site, '_', timestr, 
                    '_', ifelse(tropomiTF, 'TROPOMI', oco.sensor))


  # Startup messages -----------------------------------------------------------
  message('\n\nInitializing X-STILT')
  message('Number of receptors: ', nrow(receptors))
  message('Number of parallel threads: ', n_nodes * n_cores)


  # Source dependencies --------------------------------------------------------
  setwd(xstilt_wd)
  source('r/dependencies.r')


  # Structure out directory ------------------------------------------------------
  # Outputs are organized in three formats. by-id contains simulation files by
  # unique simulation identifier. particles and footprints contain symbolic links
  # to the particle trajectory and footprint files in by-id
  system(paste0('rm -r ', output_wd, '/footprints'), ignore.stderr = T)
  if (run_trajec) {
    system(paste0('rm -r ', output_wd, '/by-id'), ignore.stderr = T)
    system(paste0('rm -r ', output_wd, '/met'), ignore.stderr = T)
    system(paste0('rm -r ', output_wd, '/particles'), ignore.stderr = T)
  }
  for (d in c('by-id', 'particles', 'footprints')) {
    d <- file.path(output_wd, d)
    if (!file.exists(d))
      dir.create(d, recursive = T)
  }


  # removed 'isot', 'iconvert', 'outfrac', 'random' -- 
  # not used by the latest STILT
  output <- xstilt_apply(FUN = simulation_step,
                         slurm = slurm, 
                         slurm_options = slurm_options,
                         n_cores = n_cores,
                         n_nodes = n_nodes,
                         jobname = jobname, 
                         before_footprint = list(before_footprint),
                         before_trajec = list(before_trajec),
                         lib.loc = lib.loc,
                         capemin = capemin,
                         cmass = cmass,
                         conage = conage,
                         cpack = cpack,
                         delt = delt,
                         dxf = dxf,
                         dyf = dyf,
                         dzf = dzf,
                         efile = efile,
                         emisshrs = emisshrs,
                         frhmax = frhmax,
                         frhs = frhs,
                         frme = frme,
                         frmr = frmr,
                         frts = frts,
                         frvs = frvs,
                         hnf_plume = hnf_plume,
                         horcoruverr = horcoruverr,
                         horcorzierr = horcorzierr,
                         hscale = hscale,
                         ichem = ichem,
                         idsp = idsp,
                         initd = initd,
                         k10m = k10m,
                         kagl = kagl,
                         kbls = kbls,
                         kblt = kblt,
                         kdef = kdef,
                         khinp = khinp,
                         khmax = khmax,
                         kmix0 = kmix0,
                         kmixd = kmixd,
                         kmsl = kmsl,
                         kpuff = kpuff,
                         krand = krand,
                         krnd = krnd,
                         kspl = kspl,
                         kwet = kwet,
                         kzmix = kzmix,
                         maxdim = maxdim,
                         maxpar = maxpar,
                         met_file_format = met_file_format,
                         met_path = met_path,
                         met_subgrid_buffer = met_subgrid_buffer,
                         met_subgrid_enable = met_subgrid_enable,
                         met_subgrid_levels = met_subgrid_levels,
                         mgmin = mgmin,
                         n_hours = n_hours,
                         n_met_min = n_met_min,
                         ncycl = ncycl,
                         ndump = ndump,
                         ninit = ninit,
                         nstr = nstr,
                         nturb = nturb,
                         numpar = numpar,
                         nver = nver,
                         outdt = outdt,
                         output_wd = output_wd,
                         p10f = p10f,
                         pinbc = pinbc,
                         pinpf = pinpf,
                         poutf = poutf,
                         wion = projection,
                         qcycle = qcycle,
                         r_run_time = receptors$run_time,
                         r_lati = receptors$lati,
                         r_long = receptors$long,
                         r_zagl = receptors$zagl,
                         rhb = rhb,
                         rht = rht,
                         rm_dat = rm_dat,
                         run_foot = run_foot,
                         run_trajec = run_trajec,
                         siguverr = siguverr,
                         sigzierr = sigzierr,
                         smooth_factor = smooth_factor,
                         splitf = splitf,
                         stilt_wd = xstilt_wd,
                         time_integrate = time_integrate,
                         timeout = timeout,
                         tkerd = tkerd,
                         tkern = tkern,
                         tlfrac = tlfrac,
                         tluverr = tluverr,
                         tlzierr = tlzierr,
                         tout = tout,
                         tratio = tratio,
                         tvmix = tvmix,
                         varsiwant = list(varsiwant),
                         veght = veght,
                         vscale = vscale,
                         vscaleu = vscaleu,
                         vscales = vscales,
                         w_option = w_option,
                         wbbh = wbbh,
                         wbwf = wbwf,
                         wbwr = wbwr,
                         wvert = wvert,
                         xmn = xmn,
                         xmx = xmx,
                         xres = xres,
                         ymn = ymn,
                         ymx = ymx,
                         yres = yres,
                         zicontroltf = zicontroltf,
                         ziscale = ziscale,
                         z_top = z_top,
                         zcoruverr = zcoruverr, 

  # pass additional variables to stilt_apply and then to simulation_step() 
  # needed for before_*_xstilt() for X-STILT, DW, 02/11/2019
                         oco.path = oco.path, 
                         tropomiTF = tropomiTF,
                         tropomi.speci = list(tropomi.speci),
                         tropomi.fn = list(tropomi.fn), 
                         ak.wgt = namelist$ak.wgt, 
                         pwf.wgt = namelist$pwf.wgt, 
                         ctm_name = namelist$ctm_name, 
                         ctm_path = namelist$ctm_path, 
                         ctm_file_format = namelist$ctm_file_format, 
                         run_hor_err = run_hor_err,
                         emiss.file = emiss.file, 
                         ct.ver = namelist$ct.ver, 
                         ctflux.path = namelist$ctflux.path, 
                         ctmole.path = namelist$ctmole.path, 
                         xres2 = list(namelist$foot_res2), 
                         yres2 = list(namelist$foot_res2), 
                         foot_nhrs = namelist$foot_nhrs )
}


# for debugging
if (F) {
  
  X = 1
  r_run_time = receptors$run_time[X]
  r_lati = receptors$lati[X]
  r_long = receptors$long[X]
  r_zagl = receptors$zagl[X]

  args <- list(oco.path = oco.path, 
              tropomiTF = tropomiTF,
              tropomi.speci = list(tropomi_speci),
              tropomi.fn = list(tropomi.fn), 
              ak.wgt = namelist$ak.wgt, 
              pwf.wgt = namelist$pwf.wgt, 
              ctm_name = namelist$ctm_name, 
              ctm_path = namelist$ctm_path, 
              ctm_file_format = namelist$ctm_file_format, 
              run_hor_err = run_hor_err,
              emiss.file = emiss.file, 
              ct.ver = namelist$ct.ver, 
              ctflux.path = namelist$ctflux.path, 
              ctmole.path = namelist$ctmole.path, 
              xres2 = list(namelist$foot_res2), 
              yres2 = list(namelist$foot_res2), 
              foot_nhrs = foot_nhrs)

}

